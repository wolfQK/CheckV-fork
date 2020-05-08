from checkv.modules import (
    download_database,
    contamination,
    completeness,
    repeats,
    quality_summary,
)

try:
    from importlib import metadata
except ImportError:
    import importlib_metadata as metadata

__version__ = metadata.version("checkv")
