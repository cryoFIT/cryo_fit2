from __future__ import division, print_function
import subprocess, sys

# this is needed to import all common functions
path = subprocess.check_output(["which", "phenix.cryo_fit2"])
splited = path.split("/")
command_path = ''
for i in range(len(splited)-3):
  command_path = command_path + splited[i] + "/"
command_path = command_path + "modules/cryo_fit2/"
command_line_path = command_path + "command_line/"
sys.path.insert(0, command_line_path)

import cryo_fit2_program
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
    relative_path="cryo_fit2/regression/input/tst_cryo_fit2_helix_1_8.pdb",
    test=os.path.isfile)
  map_file = libtbx.env.find_in_repositories(
    relative_path="cryo_fit2/regression/input/tst_cryo_fit2_helix_1_8_reso_10.ccp4",
    test=os.path.isfile)
  resolution = "resolution=10"
  strong_ss = "strong_ss=False"
  explore = "explore=False"
  total_steps = "total_steps=50"
  total_steps_for_exploration = "total_steps_for_exploration=50"
  start_temperature = "start_temperature=300"
  final_temperature = "final_temperature=280"
  MD_in_each_cycle = "MD_in_each_cycle=2"
  number_of_steps = "number_of_steps=50"
  assert (not None in [pdb_file, map_file])
  cryo_fit2_results = run_program(program_class=cryo_fit2_program.Program, \
                                  args=[pdb_file, map_file, resolution, start_temperature, final_temperature, \
                                        MD_in_each_cycle, number_of_steps, strong_ss, explore, \
                                        total_steps, total_steps_for_exploration])
############## end of exercise_cryo_fit2()


if __name__=='__main__':
  keep_going=True
  try:
    import wx # special import
  except ImportError:
    print("Required cctbx irrelevant dependencies are missing, skipping test.")
    keep_going=False
  tstdir = libtbx.env.find_in_repositories("cryo_fit2/regression")
  if (tstdir is None) :
    warnings.warn("phenix_regression not available, skipping test")
  else :
    if(keep_going):
      exercise_cryo_fit2()      
      print("OK")
