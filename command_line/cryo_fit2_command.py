# this .py is essential for commandline cryo_fit2 running

# LIBTBX_SET_DISPATCHER_NAME phenix.cryo_fit2 # this is essential to recognize phenix.cryo_fit2 10/15/2018
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

from __future__ import division, print_function
from iotbx.cli_parser import run_program
import mmtbx.restraints
from mmtbx import monomer_library
import subprocess, sys

import cryo_fit2_program # ImportError: No module named cryo_fit2.programs

if __name__ == '__main__':
  run_program(program_class=cryo_fit2_program.Program)
