import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path

_FMT = logging.Formatter(
    "%(asctime)s | %(levelname)s | %(name)s | "
    "%(filename)s:%(lineno)d | %(message)s"
)

def setup_logging(
    *,
    level=logging.INFO,
    log_dir="logs",
    filename="pipeline.log",
    console=False,
    console_level=None,
    max_bytes=5_000_000,
    backup_count=5,
):
    """
    Simple logging setup:
    - always logs to file
    - optional console output
    - safe to call multiple times (no duplicate handlers)
    """

    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / filename

    root = logging.getLogger()
    root.setLevel(level)

    # Avoid duplicate handlers
    if not root.handlers:
        # File handler
        file_handler = RotatingFileHandler(
            log_file,
            maxBytes=max_bytes,
            backupCount=backup_count,
            encoding="utf-8"
        )
        file_handler.setFormatter(_FMT)
        root.addHandler(file_handler)

        # Optional console
        if console:
            console_handler = logging.StreamHandler()
            console_handler.setFormatter(_FMT)
            console_handler.setLevel(console_level or level)
            root.addHandler(console_handler)


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(name)
