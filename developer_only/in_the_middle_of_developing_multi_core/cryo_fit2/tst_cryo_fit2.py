# LIBTBX_SET_DISPATCHER_NAME phenix.cryo_fit2.run_tests
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

from __future__ import division, print_function
from mmtbx.programs import cryo_fit2
from iotbx.file_reader import any_file
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.utils import null_out, Sorry
import libtbx.load_env
import warnings
import os.path
from iotbx.cli_parser import run_program

def exercise_cryo_fit2(): #Checks that cryo_fit2 runs well
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/mmtbx/cryo_fit2/tst_cryo_fit2_model.pdb",
    test=os.path.isfile)
  map_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/mmtbx/cryo_fit2/tst_cryo_fit2_map.ccp4",
    test=os.path.isfile)
  assert (not None in [pdb_file, map_file])
  cryo_fit2_results = run_program(program_class=cryo_fit2.Program, args=[pdb_file, map_file, 'quiet=True'])

if __name__=='__main__':
  keep_going=True
  try:
    import wx # special import
  except ImportError:
    print("Required cctbx irrelevant dependencies are missing, skipping test.")
    keep_going=False
  tstdir = libtbx.env.find_in_repositories("phenix_regression/mmtbx/cryo_fit2")
  if (tstdir is None) :
    warnings.warn("phenix_regression not available, skipping test")
  else :
    if(keep_going):
      exercise_cryo_fit2()      
      print("OK")
