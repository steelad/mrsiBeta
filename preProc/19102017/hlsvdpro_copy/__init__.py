from __future__ import division

# 3rd party imports
import pkg_resources

# This removes a useless layer of naming indirection.
from hlsvdpro.hlsvd import *

__version__ = pkg_resources.get_distribution('hlsvdpro').version
