import importlib_metadata

try:
    __version__ = importlib_metadata.version(__package__)  # type: ignore
except Exception:  # pragma: no cover
    # package is not installed
    __version__ = "0.0.0"


from .logo import *  # noqa: F401, F403
from .logo_formatter import *  # noqa: F401, F403
from .seq import *  # noqa: F401, F403
from .seq_io import *  # noqa: F401, F403
