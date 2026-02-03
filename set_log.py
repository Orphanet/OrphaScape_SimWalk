# set_log.py — Simplifié et optimisé
"""
Configuration du logging centralisée.
Console OFF par défaut, idempotent, pas de re-config dans get_logger.
"""
from pathlib import Path
import logging
from logging.handlers import RotatingFileHandler
from typing import Optional

# Sera initialisé par setup_logging ou au premier appel get_logger
_LOG_DIR: Optional[Path] = None
_LOG_FILE: Optional[Path] = None
_CONFIGURED = False

_FMT = logging.Formatter(
    "%(asctime)s | %(levelname)s | %(name)s | %(filename)s:%(lineno)d | %(message)s"
)


def setup_logging(
    *,
    level: int = logging.INFO,
    console: bool = False,
    console_level: Optional[int] = None,
    filename: Optional[str] = None,
    log_dir: Optional[Path] = None,
    max_bytes: int = 5_000_000,
    backup_count: int = 5,
) -> None:
    """
    Configure le logging: fichier toujours présent, console optionnelle.
    Ré-appeler met à jour la console sans dupliquer les handlers.
    
    Args:
        level: Niveau de log (default: INFO)
        console: Activer la sortie console (default: False)
        console_level: Niveau spécifique pour la console
        filename: Nom du fichier de log (sera placé dans log_dir)
        log_dir: Dossier des logs (default: PROJECT_ROOT/logs)
        max_bytes: Taille max avant rotation
        backup_count: Nombre de fichiers de backup
    """
    global _LOG_DIR, _LOG_FILE, _CONFIGURED
    
    # Import conditionnel pour éviter import circulaire
    try:
        import path_variable as PV
        default_log_dir = Path(PV.PROJECT_ROOT) / "logs"
    except ImportError:
        default_log_dir = Path.cwd() / "logs"
    
    _LOG_DIR = log_dir or default_log_dir
    _LOG_DIR.mkdir(parents=True, exist_ok=True)
    
    # Détermine le fichier de log
    if filename:
        p = Path(filename)
        if p.parent != Path("."):
            # Chemin absolu fourni
            p.parent.mkdir(parents=True, exist_ok=True)
            _LOG_FILE = p
        else:
            _LOG_FILE = _LOG_DIR / p.name
    else:
        _LOG_FILE = _LOG_DIR / "pipeline.log"
    
    _LOG_FILE.touch(exist_ok=True)
    
    root = logging.getLogger()
    root.setLevel(level)
    
    # FileHandler unique
    has_file_handler = any(
        isinstance(h, logging.FileHandler) 
        for h in root.handlers
    )
    if not has_file_handler:
        fh = RotatingFileHandler(
            _LOG_FILE,
            maxBytes=max_bytes,
            backupCount=backup_count,
            encoding="utf-8"
        )
        fh.setFormatter(_FMT)
        root.addHandler(fh)
    
    # Gestion console: enlève les anciennes, ajoute si demandé
    for h in list(root.handlers):
        if isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler):
            root.removeHandler(h)
    
    if console:
        ch = logging.StreamHandler()
        ch.setFormatter(_FMT)
        if console_level is not None:
            ch.setLevel(console_level)
        root.addHandler(ch)
    
    _CONFIGURED = True


def get_logger(name: str = "app") -> logging.Logger:
    """
    Retourne un logger configuré.
    Si setup_logging() n'a pas été appelé, configure avec les valeurs par défaut.
    """
    if not _CONFIGURED:
        setup_logging()
    return logging.getLogger(name)


def get_log_file_path() -> str:
    """Retourne le chemin absolu du fichier de log."""
    if _LOG_FILE is None:
        setup_logging()
    return str(_LOG_FILE.resolve()) if _LOG_FILE else ""
