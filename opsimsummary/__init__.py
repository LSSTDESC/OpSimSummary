from __future__ import absolute_import
import os
from .summarize_opsim import *
from . import summarize_opsim
from .simlib import *
from .version import __VERSION__ as __version__
from .trig import *
# from .tessellations import *
from .opsim_out import *

here = __file__
basedir = os.path.split(here)[0]
example_data = os.path.join(basedir, 'example_data')

# from .healpix import *
# from .healpixTiles import *
# import os

# here = __file__
# basedir = os.path.split(here)[0]
# example_data = os.path.join(basedir, 'example_data')
