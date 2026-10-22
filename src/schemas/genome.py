"""宿主基因组及其已有的 EC 注释。"""

from dataclasses import dataclass, field


@dataclass
class Genome:
    """保存宿主名称、来源 GBFF 路径和用于匹配的 EC 集合。

    ec_numbers 由解析器从非伪基因 CDS 中提取并去重。
    本对象只保存数据，不负责读取或解析基因组文件。
    """

    file_path: str
    name: str | None = None
    ec_numbers: set[str] = field(default_factory=set)
