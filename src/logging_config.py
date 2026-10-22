"""项目独立的中文控制台日志，不改变调用方的根日志配置。"""

from copy import copy
import logging


class _ChineseFormatter(logging.Formatter):
    def format(self, record):
        localized = copy(record)
        localized.levelname = {
            logging.DEBUG: "调试", logging.INFO: "信息", logging.WARNING: "警告",
            logging.ERROR: "错误", logging.CRITICAL: "严重错误",
        }.get(record.levelno, record.levelname)
        return super().format(localized)


def configure_logging(level: int | str = logging.INFO) -> logging.Logger:
    """默认实时输出到控制台；重复调用不重复添加处理器，可调整日志级别。"""
    logger = logging.getLogger("pathfiller")
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(_ChineseFormatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S"))
        logger.addHandler(handler)
    logger.setLevel(level)
    logger.propagate = False
    return logger


def get_logger(name: str) -> logging.Logger:
    logger = logging.getLogger("pathfiller")
    if not logger.handlers:
        configure_logging()
    return logger.getChild(name)
