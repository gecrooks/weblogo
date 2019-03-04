

from setuptools_scm import get_version

try:
    __version__ = get_version(root='..', relative_to=__file__)
except:
    __version__ = '?.?.?'

from .logo import *								# noqa: F401, F403
from .seq import *								# noqa: F401, F403
from .seq_io import *							# noqa: F401, F403
from .logo_formatter import *					# noqa: F401, F403
