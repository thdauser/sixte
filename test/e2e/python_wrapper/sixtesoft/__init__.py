"""
Define functions for the user of the sixtesoft package
"""

# from .sixte import *

# Generic functions
from .sixte import enable_sixte_slurm_use

# Low-level SIXTE functions
from .sixte import phogen
from .sixte import phoimg
from .sixte import gendetsim

# High-level SIXTE functions
from .sixte import runsixt
from .sixte import exposure_map
from .sixte import makelc
from .sixte import imgev

# Mission specific SIXTE functions
from .sixte import erosim
from .sixte import ero_vis
from .sixte import ero_calevents
