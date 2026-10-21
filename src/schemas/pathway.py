"""由有序反应连接组成的一条通路。"""

from dataclasses import dataclass, field

from .connection import Connection


@dataclass
class PathWay:
    """connections 按照上游到目标化合物的顺序保存。"""

    name: str | None = None
    connections: list[Connection] = field(default_factory=list)
