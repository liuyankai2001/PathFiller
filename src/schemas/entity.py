"""通路中的化合物实体。"""

from dataclasses import dataclass


@dataclass
class Entity:
    """一个 C 编号对应的化合物，可被多条连接共同引用。"""

    id: str
    name: str | None = None
    is_tca: bool = False
