import subprocess, sys
import iotbx.phil, libtbx
from libtbx import group_args
from libtbx.utils import Sorry
from libtbx.utils import null_out

# this is needed to import util.py
path = subprocess.check_output(["which", "phenix.cryo_fit2"])
splited_path = path.split("/")
command_path = ''
for i in range(len(splited_path)-3):
  command_path = command_path + splited_path[i] + "/"
command_path = command_path + "modules/cryo_fit2/"
util_path = command_path + "util/"
sys.path.insert(0, util_path)
#print ("util_path:",util_path)
from util import *

def rewrite_custom_geometry(args):
  input_model_file_name = args[0]
  ######## produce pymol format secondary structure restraints #########
  # I heard that running phenix commandline directly is not ideal.
  # Therefore, I had used code directly rather than executing phenix executables at commandline such as calculating rmsd
  # However, I think that running phenix.secondary_structure_restraints is the best option here.
  # The reason is that I need to copy most of the codes in cctbx_project/mmtbx/command_line/secondary_structure_restraints.py
  #to use codes directly instead of running executables at commandline
  print("\nGenerate default secondary structure restraints for user input model file to enforce stronger secondary structure restraints\n")
  make_pymol_ss_restraints = "phenix.secondary_structure_restraints " + input_model_file_name + " format=pymol"
  libtbx.easy_run.fully_buffered(make_pymol_ss_restraints)
  
  splited_input_model_file_name = input_model_file_name.split("/")
  input_model_file_name_wo_path = splited_input_model_file_name[len(splited_input_model_file_name)-1]
  ss_restraints_file_name = input_model_file_name_wo_path + "_ss.pml"
  rewrite_pymol_ss_to_custom_geometry_ss(ss_restraints_file_name)
########### end of rewrite_custom_geometry(args)


if (__name__ == "__main__"):
  args = sys.argv[1:]
  if len(args) == 0:
    print "Provide xxx.pdb file"
    print "Example: phenix.python make_custom_geometry.py xxx.pdb"
    exit(1)
  rewrite_custom_geometry(args)
  print "OK"
