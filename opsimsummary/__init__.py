from __future__ import absolute_import
import os
from .summarize_opsim import *
from . import summarize_opsim
from .simlib import SimlibMixin
from .version import __VERSION__ as __version__
here = __file__
basedir = os.path.split(here)[0]
example_data = os.path.join(basedir, 'example_data')
from .trig import *
from .healpix import *
from .tessellations import *
from .opsim_out import *
from .healpixTiles import *
import os

here = __file__
basedir = os.path.split(here)[0]
example_data = os.path.join(basedir, 'example_data')
try:
    from .visualization import *
except ImportError as e:
    print('Some imports failed, which implies some dependencies are missing as described below')
    print(e)
    print('Visulization functions based on maps will not work')
    pass
