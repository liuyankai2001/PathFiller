"""化合物之间的反应连接。"""

from dataclasses import dataclass, field
from typing import Literal

from .entity import Entity


@dataclass
class Connection:
    """保存一步图示转化及其关联的反应和 EC。

    source、target 按 map 中的底物、产物确定，直接引用 Entity 对象。
    direction 取自当前 map 的 KGML reaction.type，不用 R 条目的方程覆盖。
    substrate_ids、product_ids 保留该反应在图中列出的全部底物、产物。
    ec_numbers 保存关联 EC，不表示必须同时补充全部 EC。
    is_endogenous 为布尔值，默认 False；与 Genome 匹配后，
    True 表示至少一个 EC 有宿主注释支持，False 表示未匹配。
    """

    source: Entity
    target: Entity
    reaction_ids: list[str] = field(default_factory=list)
    ec_numbers: set[str] = field(default_factory=set)
    direction: Literal["irreversible", "reversible", "unknown"] = "unknown"
    substrate_ids: list[str] = field(default_factory=list)
    product_ids: list[str] = field(default_factory=list)
    is_endogenous: bool = False
