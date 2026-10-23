"""PathFiller 命令行：宿主解析、组网、前 N 套不同补酶方案及统一输出。"""

import argparse
from pathlib import Path
import re

from src.genome_parser import parse_genome
from src.kegg import DEFAULT_CACHE_DIR
from src.logging_config import configure_logging, get_logger
from src.output_paths import default_output_directory, resolve_output_directory
from src.result_export import check_output_directory, export_results
from src.search import build_pathway_network, find_minimum_ec_paths

ROOT = Path(__file__).resolve().parent
logger = get_logger("cli")


def positive_int(value):
    try:
        number = int(value)
    except (TypeError, ValueError) as error:
        raise argparse.ArgumentTypeError("方案数量必须是正整数") from error
    if number < 1:
        raise argparse.ArgumentTypeError("方案数量必须是正整数")
    return number


def print_summary(network, result, output):
    print(f"目标：{network.target_id}；请求 {result.requested_count} 套方案，返回 {len(result.candidates)} 套。")
    print(f"数据完整：{result.data_is_complete}；前 N 套判定完成：{result.search_complete}；全部方案穷尽：{result.exhausted}。")
    for rank, candidate in enumerate(result.candidates, 1):
        start = candidate.pathway.connections[0].source.id if candidate.pathway.connections else network.target_id
        print(f"方案 {rank}：{start} → {network.target_id}；{len(candidate.pathway.connections)} 步；"
              f"需补 {len(candidate.extra_ecs)} 种 EC：{', '.join(sorted(candidate.extra_ecs)) or '无需补充'}")
    if not result.data_is_complete:
        print("数据不完整，以上是当前网络内的候选；其他原料供应和碳流尚未验证。")
    print("报告：", (Path(output) / "report.md").resolve())


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", nargs="?", default="C05432", help="目标 KEGG 化合物编号，默认 C05432")
    parser.add_argument("--gbff", type=Path, default=ROOT / "data" / "genomic.gbff", help="宿主 GBFF 文件")
    parser.add_argument("--top-n", type=positive_int, default=3, help="输出前几套不同补酶方案，默认 3；相同 EC 集合只计一套")
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--verbose", action="store_true", help="显示请求、缓存和逐图处理的调试日志")
    parser.add_argument("--export-network-svg", action="store_true", help="额外导出按候选着色的完整网络 SVG")
    parser.add_argument("--output-dir", type=Path, help="outputs 内的输出目录；已有同名业务结果时拒绝覆盖")
    args = parser.parse_args(argv)
    configure_logging("DEBUG" if args.verbose else "INFO")
    target = args.target.strip().upper()
    if re.fullmatch(r"C[0-9]{5}", target) is None:
        parser.error("目标编号应为 C 加五位数字，例如 C05432")
    try:
        output = check_output_directory(resolve_output_directory(args.output_dir) if args.output_dir
                                        else default_output_directory(target))
        logger.info("开始分析 %s：请求前 %d 套不同补酶方案", target, args.top_n)
        logger.debug("输出目录：%s", output)
        genome = parse_genome(args.gbff)
        network = build_pathway_network(target, cache_dir=args.cache_dir, refresh=args.refresh)
        result = find_minimum_ec_paths(network, genome, top_n=args.top_n)
        export_results(network, genome, result, output, cache_dir=args.cache_dir, refresh=args.refresh,
                       include_network_svg=args.export_network_svg)
        print_summary(network, result, output)
    except (OSError, ValueError, RuntimeError) as error:
        logger.error("运行失败：%s", error)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
