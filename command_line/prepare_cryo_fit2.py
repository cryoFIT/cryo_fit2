######to do list
# To run cryo_fit2 at windows as well, replace all unix command like libtbx.easy_run.fully_buffered

from __future__ import division, print_function

from cryo_fit2_util import *

from iotbx import file_reader
from iotbx.xplor import crystal_symmetry_from_map

from libtbx import easy_mp
from libtbx import easy_pickle
from libtbx import phil
import iotbx.pdb
import iotbx.phil
from libtbx.phil import change_default_phil_values
import libtbx.phil
import libtbx.phil.command_line
from libtbx.utils import date_and_time, multi_out, Sorry

import mmtbx
import cryo_fit2_run
import os, shutil, subprocess, sys, time

try:
  from phenix.program_template import ProgramTemplate # at Doonam's laptop, phenix.program_template works
except ImportError:
  from libtbx.program_template import ProgramTemplate

# this is needed to import util py files
path = subprocess.check_output(["which", "phenix.cryo_fit2"])
splited_path = path.split("/")
command_path = ''
for i in range(len(splited_path)-3):
  command_path = command_path + splited_path[i] + "/"
command_path = command_path + "modules/cryo_fit2/"
util_path = command_path + "util/"
sys.path.insert(0, util_path)
print ("util_path:",util_path)
from util import *

base_master_phil_str = '''
include scope libtbx.phil.interface.tracking_params
start_temperature = 300
  .type = int
  .short_caption = Starting temperature of annealing in Kelvin
final_temperature = 0
  .type = int
  .short_caption = Final temperature of annealing in Kelvin
cool_rate = 10
  .type = int
  .short_caption = cooling rate of annealing in Kelvin
number_of_steps = 1000
  .type = int
  .short_caption = number of steps in phenix.dynamics
map_weight = None
  .type = float
  .short_caption = cryo-EM map weight. \
                   A user is recommended NOT to specify this, so that it will be automatically optimized. \
                   If the map is derived from SAXS, map_weight < 0.3 is recommended so that base pairs of nucleic acids are intact.
resolution = None
  .type = float
  .short_caption = cryo-EM map resolution (angstrom) that needs to be specified by a user
output_dir = output
  .type = path
  .short_caption = Output folder PREFIX
progress_on_screen = True
    .type          = bool
    .help          = If True, temp=xx dist_moved=xx angles=xx bonds=xx is shown on screen rather than cryo_fit2.log \
                     If False, temp=xx dist_moved=xx angles=xx bonds=xx is NOT shown on screen, and saved into cryo_fit2.log
loose_ss_def = False
    .type   = bool
    .help   = If True, secondary structure definition for nucleic acid is loose. Use this with great caution.  \
              If False, use Oleg's original strict definition. 
    .short_caption = Keep origin of a resulted atomic model
keep_origin = True
    .type   = bool
    .help   = If True, write out model with origin in original location.  \
              If False, shift map origin to (0,0,0). 
    .short_caption = Keep origin of a resulted atomic model
devel = False
    .type   = bool
    .help   = If True, run quickly only to check sanity
include scope mmtbx.monomer_library.pdb_interpretation.grand_master_phil_str # to use secondary_structure.enabled
include scope mmtbx.monomer_library.pdb_interpretation.geometry_restraints_remove_str # to use nucleic_acid.base_pair.restrain_planarity but not works as expected

selection_fixed = None
  .type = str
  .short_caption = Selection for fixed model
  .input_size = 400
  .help = Selection of the target atoms to fit to (optional)

selection_moving = None
  .type = str
  .short_caption = Selection for moving model
  .input_size = 400
  .help = Selection of the atoms that will be fit to selection_fixed (optional)
  
selection_fixed_preset = * ca backbone all
  .type = choice
  .help = Selection preset for fixed model.
  
selection_moving_preset = * ca backbone all
  .type = choice
  .help = Selection preset for moving model.
  
''' ############## end of base_master_phil_str



new_default = 'pdb_interpretation.secondary_structure.enabled = True'
modified_master_phil_str = change_default_phil_values(
  base_master_phil_str, new_default, phil_parse=iotbx.phil.parse)

new_default = 'pdb_interpretation.secondary_structure.protein.remove_outliers = True'
modified_master_phil_str = change_default_phil_values(
  modified_master_phil_str, new_default, phil_parse=iotbx.phil.parse)


class Program(ProgramTemplate):

  description = '''
Program for running cryo_fit2.\n

Minimum required inputs:
  Model file (.pdb or .cif)
  Map file   (MRC/ccp4 format)
  Map resolution

Example running command:
  phenix.cryo_fit2 model.pdb map.ccp4 resolution=5

Options:
  resolution                   (cryo-EM map resolution in angstrom that needs to be entered by a user)
  map_weight                   (cryo-EM map weight.
                               A user is recommended NOT to specify this, so that it will be automatically optimized.
                               If the map is derived from SAXS, map_weight < 0.3 is recommended so that base pairs of nucleic acids are intact.)
  start_temperature            (default: 300)
  final_temperature            (default: 0)
  cool_rate                    (default: 10)
  number_of_steps              (default: 1000)
  secondary_structure.enabled  (default: True)
                               Most MD simulations tend to break secondary structure. 
                               Therefore, turning on this option is recommended. 
                               If HELIX/SHEET records are present in supplied .pdb file, 
                               automatic search of the existing secondary structures in the given 
                               input pdb file will not be executed.
  secondary_structure.protein.remove_outliers (default: True)
                               False may be useful for very poor low-resolution structures by
                               ignoring some hydrogen "bond" if it exceed certain distance threshold
  output_dir                   (output folder name prefix, default: output)
  keep_origin                  (default: True)
                               If True, write out model with origin in original location.
                               If False, shift origin to (0,0,0). 
  progress_on_screen           (default: False)
                               If True, temp= xx dist_moved= xx angles= xx bonds= xx is shown on screen rather than cryo_fit2.log 
                               If False, temp= xx dist_moved= xx angles= xx bonds= xx is NOT shown on screen, and saved into cryo_fit2.log
  devel                        (default: False)
                               If True, run quickly only to check sanity
'''

  datatypes = ['model', 'real_map', 'phil']
  
  master_phil_str = modified_master_phil_str # this is ESSENTIAL to avoid
  '''
     Sorry: Some PHIL parameters are not recognized by phenix.cryo_fit2.
     Please run this program with the --show-defaults option to see what parameters are available.
     PHIL parameters in files should be fully specified (e.g. "output.overwrite" instead of just "overwrite")
  '''
  
  # ---------------------------------------------------------------------------
  def validate(self): # this validate function runs as well (1/26/2019)
    print('Validating inputs', file=self.logger)
    if not (self.data_manager.has_models()):
      raise Sorry("Supply an atomic model file. Type \"phenix.cryo_fit2\" to know minimally required options")
    if not (self.data_manager.has_real_maps()):
      raise Sorry("Supply a map file. Type \"phenix.cryo_fit2\" to know minimally required options")
    if not (self.params.resolution):
      raise Sorry("Map resolution is required. A user can get the resolution either by EMDB reported value or by running phenix.mtriage. Type \"phenix.cryo_fit2\" to know minimally required options")

  # ---------------------------------------------------------------------------
  def run(self):
    args = sys.argv[1:]
    print ("args", args)
    checked_whether_args_has_eff = check_whether_args_has_eff(args)
    
    input_model_file_name = self.data_manager.get_default_model_name()
    print('User input model: %s' % input_model_file_name, file=self.logger)
    
    log = multi_out()
    out=sys.stdout
    log.register("stdout", out)
    
    splited = input_model_file_name.split("/")
    input_model_file_name_wo_path = splited [len(splited)-1]
    
    log_file_name = "cryo_fit2.log"
    
    logfile = open(log_file_name, "w") # since it is 'w', an existing file with the same name will be erased
    log.register("logfile", logfile)
    
    if (checked_whether_args_has_eff == True):
      write_this = "User entered custom geometry restraints already.\n"
      print (write_this)
      logfile.write(write_this)
    else:
      write_this = "User did not enter custom geometry restraints, make it now.\n"
      print (write_this)
      logfile.write(write_this)
      
      ######## produce pymol format secondary structure restraints #########
      # I heard that running phenix commandline directly is not ideal.
      # Therefore, I had used code directly rather than executing phenix executables at commandline such as calculating rmsd
      # However, I think that running phenix.secondary_structure_restraints is the best option here.
      # The reason is that I need to copy most of the codes in cctbx_project/mmtbx/command_line/secondary_structure_restraints.py
      #to use codes directly instead of running executables at commandline
      logfile.write("\nGenerate default secondary structure restraints for user input model file to enforce stronger secondary structure restraints\n")
      make_pymol_ss_restraints = "phenix.secondary_structure_restraints " + input_model_file_name + " format=pymol"
      logfile.write(make_pymol_ss_restraints)
      logfile.write("\n")
      libtbx.easy_run.fully_buffered(make_pymol_ss_restraints)
  
      splited_input_model_file_name = input_model_file_name.split("/")
      input_model_file_name_wo_path = splited_input_model_file_name[len(splited_input_model_file_name)-1]
      ss_restraints_file_name = input_model_file_name_wo_path + "_ss.pml"
      rewrite_to_custom_geometry(ss_restraints_file_name)
      custom_geom_file_name = ss_restraints_file_name[:-4] + "_custom_geom.eff"
    
    logfile.write("\n")
    logfile.close()
