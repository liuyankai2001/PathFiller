"""查询 KEGG 通路并将图解析为化合物网络。"""

from functools import lru_cache
from http.client import IncompleteRead
from pathlib import Path
import random
import re
from threading import Lock
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

import networkx as nx

from .map_parser import parse_kgml
from .logging_config import get_logger
from .schemas import Entity


logger = get_logger("kegg")
DEFAULT_CACHE_DIR = Path(__file__).resolve().parents[1] / "cache"
_request_lock = Lock()
_last_request_started = 0.0
_MAX_ATTEMPTS = 5
_BACKOFF_BASE = 1.0
_BACKOFF_CAP = 30.0


class KeggRequestError(RuntimeError):
    """请求最终失败；status_code 为 None 表示网络错误或超时。"""

    def __init__(
        self, endpoint: str, *, status_code: int | None = None,
        attempts: int = 1, reason: str | None = None,
    ):
        self.endpoint = endpoint
        self.status_code = status_code
        self.attempts = attempts
        detail = f"HTTP {status_code}" if status_code is not None else (reason or "网络错误或超时")
        super().__init__(f"KEGG 请求 {endpoint} 失败（{detail}；已尝试 {attempts} 次）")


def _compound_name(value: str) -> str:
    """选取首个名称，拒绝空值、控制字符及被误当作名称的 C 编号。"""
    name = value.split(";", 1)[0].strip()
    if (not name or any(ord(char) < 32 for char in name)
            or re.fullmatch(r"(?:cpd:)?C[0-9]{5}", name)):
        raise ValueError("化合物名称为空、含控制字符或只有编号")
    return name


def get_compound_names(
    compound_ids, *, cache_dir: str | Path | None = None, refresh: bool = False,
) -> dict[str, str]:
    """批量获取 KEGG 首个化合物名称；单批失败只记录提示并继续。

    名称缓存为 Cxxxxx.name.txt。仅查询缺失的编号，每批最多十个；请求
    复用统一限速与重试。返回成功项，缺失项由调用方保留编号或已有名称。
    """
    values = [compound_ids] if isinstance(compound_ids, str) else compound_ids
    identifiers = set()
    for value in values:
        if not isinstance(value, str) or re.fullmatch(r"C[0-9]{5}", value.strip().upper()) is None:
            raise ValueError("化合物编号应为 C 加五位数字，例如 C05432")
        identifiers.add(value.strip().upper())
    cache = Path(cache_dir) if cache_dir is not None else DEFAULT_CACHE_DIR
    names, pending = {}, []
    for cid in sorted(identifiers):
        path = cache / f"{cid}.name.txt"
        if path.exists() and not refresh:
            try:
                names[cid] = _compound_name(path.read_text(encoding="utf-8").strip())
                continue
            except (OSError, UnicodeError, ValueError) as error:
                logger.debug("%s 名称缓存无效，将重新查询：%s", cid, error)
        pending.append(cid)
    if identifiers:
        logger.info("化合物名称：共 %d 个，命中缓存 %d 个，待查询 %d 个",
                    len(identifiers), len(names), len(pending))
    last_progress = time.monotonic()
    cache_write_failures = 0
    for offset in range(0, len(pending), 10):
        batch = pending[offset:offset + 10]
        logger.debug("查询化合物名称批次 %d/%d：%s", offset // 10 + 1,
                     (len(pending) + 9) // 10, ", ".join(batch))
        try:
            text = _request_text("list/" + "+".join(batch))
            fetched = {}
            for line in text.splitlines():
                if not line.strip():
                    continue
                fields = line.split("\t")
                cid = fields[0].removeprefix("cpd:")
                if len(fields) != 2 or cid not in batch or cid in fetched:
                    raise ValueError("KEGG 化合物名称响应包含无效、重复或未请求的编号")
                fetched[cid] = _compound_name(fields[1])
        except (RuntimeError, OSError, UnicodeError, ValueError) as error:
            logger.debug("化合物名称查询失败（%s），保留编号或已有名称：%s", ", ".join(batch), error)
            fetched = {}
        names.update(fetched)
        for cid, name in fetched.items():
            try:
                cache.mkdir(parents=True, exist_ok=True)
                (cache / f"{cid}.name.txt").write_text(name + "\n", encoding="utf-8")
            except OSError as error:
                cache_write_failures += 1
                logger.debug("%s 名称已获取，但缓存写入失败：%s", cid, error)
        for cid in sorted(set(batch) - fetched.keys()):
            logger.debug("KEGG 未返回 %s 的名称，保留编号或已有名称", cid)
        now = time.monotonic()
        processed = offset + len(batch)
        if now - last_progress >= 10 or processed == len(pending):
            logger.info("名称查询进度：已处理 %d/%d 个待查询化合物，累计获取 %d/%d 个名称",
                        processed, len(pending), len(names), len(identifiers))
            last_progress = now
    if len(names) < len(identifiers):
        logger.warning("%d 个化合物未获取到名称，保留已有名称或编号；详情见报告或调试日志",
                       len(identifiers) - len(names))
    if cache_write_failures:
        logger.warning("%d 个化合物名称已获取但未能写入缓存；详情见调试日志", cache_write_failures)
    return dict(sorted(names.items()))


def get_maps_by_compound(compound_id: str, *, refresh: bool = False) -> list[str]:
    """通过 C 编号查询关联的全部 map，保留接口返回顺序。

    例如 get_maps_by_compound("C05432") 返回 ["map00906", ...]。
    编号错误抛出 ValueError；请求失败或响应异常抛出 RuntimeError。
    成功查询但无通路关联时返回空列表。
    查询结果在当前进程中缓存，refresh=True 时重新查询。
    选择哪张图由搜索算法决定。
    """
    normalized = compound_id.strip().upper()
    if re.fullmatch(r"C[0-9]{5}", normalized) is None:
        raise ValueError("化合物编号应为 C 加五位数字，例如 C05432")
    if refresh:
        _fetch_maps.cache_clear()
    logger.debug("查询 %s 关联的 map%s", normalized, "（强制刷新）" if refresh else "")
    # 返回副本，避免调用方修改列表后影响缓存中的数据。
    maps = list(_fetch_maps(normalized))
    logger.debug("%s 关联查询完成：找到 %d 张 map", normalized, len(maps))
    return maps


def _request_text(endpoint: str) -> str:
    """重试临时请求错误；404 和其他不可重试 HTTP 错误立即抛出。"""
    global _last_request_started

    url = f"https://rest.kegg.jp/{endpoint}"
    request = Request(url, headers={"User-Agent": "PathFiller/0.1"})
    for attempt in range(1, _MAX_ATTEMPTS + 1):
        # 每次尝试都限速；退避等待在锁外执行。
        with _request_lock:
            delay = 0.4 - (time.monotonic() - _last_request_started)
            if delay > 0:
                time.sleep(delay)
            _last_request_started = time.monotonic()
        try:
            logger.debug("请求 KEGG：%s（第 %d/%d 次，超时 30 秒）", endpoint, attempt, _MAX_ATTEMPTS)
            with urlopen(request, timeout=30) as response:
                text = response.read().decode("utf-8")
            logger.debug("KEGG 请求成功：%s", endpoint)
            return text
        except HTTPError as error:
            error.close()
            retryable = error.code in {408, 429} or 500 <= error.code <= 599
            if not retryable or attempt == _MAX_ATTEMPTS:
                if error.code == 404:
                    logger.debug("KEGG 资源不存在：%s（HTTP 404，不重试）", endpoint)
                else:
                    logger.error("KEGG 请求失败：%s（HTTP %d，已尝试 %d 次，不再重试）",
                                 endpoint, error.code, attempt)
                raise KeggRequestError(
                    endpoint, status_code=error.code, attempts=attempt,
                ) from error
            reason = f"HTTP {error.code}"
        except (URLError, OSError, IncompleteRead) as error:
            # OSError 也包含套接字超时、连接中断和响应读取失败。
            if attempt == _MAX_ATTEMPTS:
                logger.error("KEGG 网络请求失败：%s（已尝试 %d 次，重试耗尽）：%s",
                             endpoint, attempt, str(error) or type(error).__name__)
                raise KeggRequestError(
                    endpoint, attempts=attempt, reason=str(error),
                ) from error
            reason = str(error) or type(error).__name__
        # Full Jitter：第 1～4 次重试前分别从 [0, 1/2/4/8] 秒采样。
        ceiling = min(_BACKOFF_CAP, _BACKOFF_BASE * 2 ** (attempt - 1))
        delay = random.uniform(0, ceiling)
        if attempt == 1:
            logger.warning("KEGG 请求暂时失败，正在自动重试：%s", endpoint)
        logger.debug("KEGG 请求暂时失败：%s（%s）；等待 %.2f 秒后重试（第 %d/%d 次）",
                     endpoint, reason, delay, attempt + 1, _MAX_ATTEMPTS)
        time.sleep(delay)


@lru_cache(maxsize=4096)
def _fetch_maps(compound_id: str) -> tuple[str, ...]:
    text = _request_text(f"link/pathway/{compound_id}")
    maps = []
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if (
            len(fields) != 2
            or fields[0] != f"cpd:{compound_id}"
            or re.fullmatch(r"path:map[0-9]{5}", fields[1]) is None
        ):
            raise RuntimeError(f"KEGG 返回了无法识别的通路列表：{compound_id}")
        map_id = fields[1].removeprefix("path:")
        if map_id not in maps:
            maps.append(map_id)
    return tuple(maps)


def _reaction_equation(text: str) -> str:
    parts = []
    reading = False
    for line in text.splitlines():
        field = line[:12].strip()
        if field == "EQUATION":
            reading = True
            parts.append(line[12:].strip())
        elif reading:
            if field:
                break
            parts.append(line[12:].strip())
    if not parts:
        raise ValueError("R 条目缺少 EQUATION 字段")
    return " ".join(parts)


def parse_map(
    map_id: str,
    *,
    entities: dict[str, Entity] | None = None,
    cache_dir: str | Path | None = None,
    refresh: bool = False,
) -> nx.MultiDiGraph:
    """输入 map 编号，返回保存 Entity、Connection 的 MultiDiGraph。

    RN 视图提供连接和方向，EC、KO 视图分别提供酶编号和 KO 注释。
    默认复用项目根目录 cache 中的原始数据；refresh=True 时重新获取。
    每次解析都创建新网络和连接，避免内源状态等修改污染后续调用。
    entities 可用于跨图复用实体；不会给已有对象增加额外属性。
    """
    normalized = map_id.strip().lower()
    if re.fullmatch(r"map[0-9]{5}", normalized) is None:
        raise ValueError("通路编号应为 map 加五位数字，例如 map00906")
    cache = Path(cache_dir) if cache_dir is not None else DEFAULT_CACHE_DIR
    logger.debug("开始解析 %s：联合读取 RN、EC、KO 视图", normalized)
    documents = {}
    for prefix in ("rn", "ec", "ko"):
        view = prefix + normalized[3:]
        kgml_file = cache / f"{view}.kgml"
        use_cache = kgml_file.exists() and not refresh
        logger.debug("%s %s 视图：%s", normalized, prefix.upper(),
                    f"读取缓存 {kgml_file.name}" if use_cache else "从 KEGG 下载")
        xml = kgml_file.read_text(encoding="utf-8") if use_cache else _request_text(f"get/{view}/kgml")
        try:
            root = ET.fromstring(xml)
        except ET.ParseError as error:
            raise ValueError(f"无法解析 KEGG 返回的 {prefix.upper()} KGML XML") from error
        if (
            root.tag != "pathway" or root.get("number") != normalized[3:]
            or root.get("name") != f"path:{view}" or root.get("org") != prefix
        ):
            raise ValueError(f"KGML 内容与请求的 {normalized} {prefix.upper()} 视图不一致")
        if not use_cache:
            cache.mkdir(parents=True, exist_ok=True)
            kgml_file.write_text(xml, encoding="utf-8")
            logger.debug("已缓存 %s：%s", view, kgml_file)
        documents[prefix] = xml

    @lru_cache(maxsize=None)
    def lookup(reaction_id: str) -> str:
        equation_file = cache / f"{reaction_id}.equation.txt"
        if equation_file.exists() and not refresh:
            logger.debug("读取反应方程缓存：%s", reaction_id)
            return equation_file.read_text(encoding="utf-8").strip()
        logger.debug("获取反应方程 %s，用于确认原图参与物", reaction_id)
        equation = _reaction_equation(_request_text(f"get/{reaction_id}"))
        cache.mkdir(parents=True, exist_ok=True)
        equation_file.write_text(equation + "\n", encoding="utf-8")
        return equation

    graph = parse_kgml(
        documents["rn"], ec_xml=documents["ec"], ko_xml=documents["ko"],
        entities=entities, reaction_lookup=lookup,
    )
    logger.debug("%s 解析完成：%d 个化合物、%d 条有向反应", normalized, len(graph), graph.number_of_edges())
    return graph


def get_tca_graph(
    *,
    entities: dict[str, Entity] | None = None,
    cache_dir: str | Path | None = None,
    refresh: bool = False,
) -> nx.MultiDiGraph:
    """获取 KEGG TCA 循环 map00020 的完整化合物网络。

    复用 parse_map 的下载、缓存和 KGML 方向处理。
    图中也包含辅因子等参与物。Entity.is_tca 保留为预览标记，
    跨图搜索通过通路链接指向 map00020 判断是否到达 TCA。
    传入共享 entities 字典，可以与目标通路复用相同的 Entity 对象。
    """
    return parse_map(
        "map00020", entities=entities, cache_dir=cache_dir, refresh=refresh,
    )
