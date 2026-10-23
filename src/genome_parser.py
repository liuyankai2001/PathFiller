"""从 GBFF 注释提取宿主 EC，生成 Genome 对象。"""

from pathlib import Path

from Bio import SeqIO

from .logging_config import get_logger
from .schemas import Genome


logger = get_logger("genome")


def _is_pseudo(feature) -> bool:
    return "pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers


def _organism(record) -> str | None:
    candidates = [record.annotations.get("organism", "")]
    candidates.extend(name for feature in record.features if feature.type == "source"
                      for name in feature.qualifiers.get("organism", []))
    return next((name.strip() for name in candidates if name and name.strip() not in {"", "."}), None)


def parse_genome(file_path: str | Path) -> Genome:
    """解析 GBFF 中的全部记录，返回宿主名称、绝对路径和去重 EC 集合。

    仅读取非伪基因 CDS 的 EC_number 注释；同时排除 CDS 自身和通过
    同一记录内 locus_tag 关联到伪基因 gene 的记录。多个 EC 全部保留，
    部分 EC 编号不扩展为通配符，也不从序列推断缺失注释。
    多条染色体/质粒记录合并到同一个 Genome，物种名称按出现顺序去重。
    有效文件没有 EC 时返回空集合；无记录或格式错误时抛出 ValueError。
    """
    path = Path(file_path).expanduser().resolve()
    genome = Genome(file_path=str(path))
    names = []
    record_count = cds_count = excluded_count = annotated_count = 0
    logger.info("开始解析宿主基因组：%s", path)
    try:
        if path.is_dir():
            raise IsADirectoryError(f"GBFF 路径指向目录，应提供文件：{path}")
        with path.open(encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                record_count += 1
                name = _organism(record)
                if name and name not in names:
                    names.append(name)
                pseudo_tags = {
                    tag for feature in record.features if feature.type == "gene" and _is_pseudo(feature)
                    for tag in feature.qualifiers.get("locus_tag", [])
                }
                for feature in record.features:
                    if feature.type != "CDS":
                        continue
                    cds_count += 1
                    if _is_pseudo(feature) or pseudo_tags.intersection(feature.qualifiers.get("locus_tag", [])):
                        excluded_count += 1
                        continue
                    ec_numbers = {ec.strip() for ec in feature.qualifiers.get("EC_number", []) if ec.strip()}
                    annotated_count += bool(ec_numbers)
                    genome.ec_numbers.update(ec_numbers)
                logger.debug("已解析基因组记录 %s；当前累计 %d 个 CDS、%d 种不同 EC",
                            record.id, cds_count, len(genome.ec_numbers))
    except OSError as error:
        logger.error("读取基因组文件失败：%s；%s", path, error)
        raise
    except ValueError as error:
        logger.error("GBFF 格式解析失败：%s；%s", path, error)
        raise ValueError(f"无法解析 GBFF 文件 {path}：{error}") from error
    if not record_count:
        logger.error("GBFF 文件中没有可解析的序列记录：%s", path)
        raise ValueError(f"GBFF 文件中没有可解析的序列记录：{path}")
    genome.name = "; ".join(names) or None
    logger.info("基因组解析完成：%d 条记录、%d 个 CDS，排除 %d 个伪基因 CDS；"
                "%d 个 CDS 含 EC 注释，得到 %d 种不同 EC",
                record_count, cds_count, excluded_count, annotated_count, len(genome.ec_numbers))
    if not genome.ec_numbers:
        logger.warning("该基因组未提取到非伪基因 CDS 的 EC 注释：%s", path)
    return genome
