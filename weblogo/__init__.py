from pkg_resources import DistributionNotFound, get_distribution

from .logo import *  # noqa: F401, F403
from .logo_formatter import *  # noqa: F401, F403
from .seq import *  # noqa: F401, F403
from .seq_io import *  # noqa: F401, F403

__version__: str

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:  # pragma: no cover
    __version__ = "?.?.?"
