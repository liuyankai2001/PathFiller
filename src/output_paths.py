"""主程序与检查脚本共用的输出目录约定。"""

from datetime import datetime
from pathlib import Path
from uuid import uuid4


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "outputs"


def resolve_output_directory(value):
    """相对目录按项目根目录解释，只允许写入 outputs 及其子目录。"""
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    path = path.resolve()
    root = DEFAULT_OUTPUT_DIR.resolve()
    if not path.is_relative_to(root):
        raise ValueError(f"输出目录必须位于项目 outputs 目录内：{root}")
    return path


def default_output_directory(target):
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S_%f") + "_" + uuid4().hex[:8]
    return resolve_output_directory(DEFAULT_OUTPUT_DIR / target / run_id)
