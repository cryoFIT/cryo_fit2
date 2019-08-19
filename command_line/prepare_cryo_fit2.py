# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_SIGNALS_DEFAULT=1
#### (source) http://cci.lbl.gov/cctbx_sources/crys3d/command_line/hklview.py

######to do list
# To run cryo_fit2 at windows as well, replace all unix command like libtbx.easy_run.fully_buffered
from __future__ import division, print_function

import math, os, shutil, subprocess, sys, time

from iotbx import file_reader
from iotbx.xplor import crystal_symmetry_from_map

from libtbx import easy_mp
from libtbx import easy_pickle
from libtbx import phil
import iotbx.pdb
import iotbx.phil
from libtbx.introspection import number_of_processors # modules/dials/command_line/batch_analysis.py
from libtbx.phil import change_default_phil_values
import libtbx.phil
import libtbx.phil.command_line
from libtbx.utils import date_and_time, multi_out, Sorry

import mmtbx
import cryo_fit2_run

from mmtbx.refinement.real_space import weight

try:
  from phenix.program_template import ProgramTemplate # at Doonam's laptop, phenix.program_template works
except ImportError:
  from libtbx.program_template import ProgramTemplate

os.environ['BOOST_ADAPTBX_FPE_DEFAULT'] = "1"
os.environ['BOOST_ADAPTBX_SIGNALS_DEFAULT'] = "1"


########## <begin> import util py files
cryo_fit2_repository_dir = libtbx.env.dist_path("cryo_fit2") # Locate phenix.cryo_fit.run_tests executable
util_path = cryo_fit2_repository_dir + "/util/"
sys.path.insert(0, util_path)
from util import *
########## <end> import util py files


program_citations = libtbx.phil.parse('''
citation {
  article_id = Pavel_s_2018_realspace_refine mentioned and used this phenix.dynamics
}
''')


base_master_phil_str = '''
include scope libtbx.phil.interface.tracking_params
cool_rate        = None
  .type          = float
  .short_caption = Cooling rate of annealing in Kelvin. Will be automatically determined by cryo_fit2.
explore          = True
  .type          = bool
  .short_caption = If True, cryo_fit2 will use maximum number of multiple cores to explore the most optimal MD parameters.\
                   However, this exploration requires a lot of computing power (e.g. > 128 GB memory, > 20 cores).\
                   Exploring at a macbook pro (16 GB memory, 2 cores) crashed.
final_temperature = 0
  .type           = float
  .short_caption  = Final temperature of annealing in Kelvin
keep_origin      = True
  .type          = bool
  .help          = If True, write out model with origin in original location.  \
                   If False, shift map origin to (0,0,0).
  .short_caption = Keep origin of a resulted atomic model
loose_ss_def = False
  .type      = bool
  .help      = If True, secondary structure definition for nucleic acid is loose. Use this with great caution.  \
               If False, use Oleg's original strict definition.
map_weight       = None
  .type          = float
  .short_caption = cryo-EM map weight. \
                   A user is recommended NOT to specify this, so that it will be automatically optimized.
max_steps_for_exploration    = 10000
  .type                      = int
  .short_caption             = The total number of steps for MD parameter exploration. \
                               10k is enough to discern Mg Channel \
                               15k is not enough for tRNA
max_steps_for_final_MD    = None
  .type                   = int
  .short_caption          = The maximum number of steps in final running of phenix.dynamics.\
                            If specified, run up to this number of steps no matter what.
MD_in_each_cycle = None
  .type          = int
  .short_caption = An cycle here is different from the one in deep learning. \
                   Here, the cycle is each iteration of MD from start_temperature to final_temperature. \
                   If not specified, cryo_fit2 will use the optimized value by automatic exploration.
nproc            = None
  .type          = int
  .short_caption = Number of cores to use for MD parameter exploration.\
                   If not specified, cryo_fit2 will use available number of cores/3
number_of_steps  = None
  .type          = int
  .short_caption = The number of MD steps in each phenix.dynamics \
                   If not specified, cryo_fit2 will use the optimized value by automatic exploration.
output_dir       = output
  .type          = path
  .short_caption = Output folder PREFIX
progress_on_screen = True
  .type            = bool
  .help            = If True,  temp=x dist_moved=x angles=x bonds=x is shown on screen rather than cryo_fit2.log \
                     If False, temp=x dist_moved=x angles=x bonds=x is NOT shown on screen, but saved into cryo_fit2.log
record_states    = False
  .type          = bool
  .help          = If True, cryo_fit2 records all states and save it to all_states.pdb (only when cryo_fit2 is successfully completed)\
                   However, 3k atoms molecules (like L1 stalk in a ribosome) require more than 160 GB of memory. \
                   If False, cryo_fit2 doesn't record each state of molecular dynamics.
reoptimize_map_weight_after_each_cycle_during_final_MD = False
  .type                                = bool
  .help                                = If True, cryo_fit2 will reoptimize map_weight after 5~100 cycles. \
                                         It will lengthens cryo_fit2 running time significantly longer.\
                                         However, Doo Nam confirmed that it is effective to prevent \
                                         nan error during core cryo-EM map based core dynamics run for full-tRNA. \
                                         However, now map_weight is not multiplied in a crazy manner, so this is False by default.
resolution       = None
  .type          = float
  .short_caption = cryo-EM map resolution (angstrom) that needs to be specified by a user
short            = False
  .type          = bool
  .help          = If True, run quickly only to check sanity
sigma_for_custom_geom = None
  .type               = float
  .short_caption      = The lower this value, the stronger the custom made secondary structure restraints will be. \
                        Oleg recommended 0.021 which is the sigma value for covalent bond.
start_temperature = None
  .type           = float
  .short_caption  = Starting temperature of annealing in Kelvin. \
                   If not specified, cryo_fit2 will use the optimized value after automatic exploration between 300 and 900.
strong_ss = True
  .type   = bool
  .help   = If True, cryo_fit2 will use a stronger sigma_for_custom_geom for secondary structure restraints. \
            If False, it will not use custom geometry
weight_multiply  = None
  .type          = float
  .short_caption = Cryo_fit2 will multiply cryo-EM map weight by this much. \ 
                   If not specified, cryo_fit2 will use the default value (e.g. 1) \
                   Usually a small molecule (a helix) requires only 1 (not multiply). \
                   For a helix, 20 keeps geometry, 100 breaks it (w/o special sigma) \
                   However, a large molecule needs a larger value (e.g. 10~50).
include scope mmtbx.monomer_library.pdb_interpretation.grand_master_phil_str # to use secondary_structure.enabled
include scope mmtbx.monomer_library.pdb_interpretation.geometry_restraints_remove_str # to use nucleic_acid.base_pair.restrain_planarity but not works as expected
selection_fixed  = None
  .type          = str
  .short_caption = Selection for fixed model
  .input_size    = 400
  .help          = Selection of the target atoms to fit to (optional)
selection_moving = None
  .type          = str
  .short_caption = Selection for moving model
  .input_size    = 400
  .help          = Selection of the atoms that will be fit to selection_fixed (optional)
selection_fixed_preset = * ca backbone all
  .type                = choice
  .help                = Selection preset for fixed model.
selection_moving_preset = * ca backbone all
  .type                 = choice
  .help                 = Selection preset for moving model.
''' ############## end of base_master_phil_str  
  

new_default = 'pdb_interpretation.secondary_structure.enabled = True'
modified_master_phil_str = change_default_phil_values(
  base_master_phil_str, new_default, phil_parse=iotbx.phil.parse)

new_default = 'pdb_interpretation.secondary_structure.protein.remove_outliers = True'
modified_master_phil_str = change_default_phil_values(
  modified_master_phil_str, new_default, phil_parse=iotbx.phil.parse)

'''
#print ("modified_master_phil_str:",modified_master_phil_str)
new_default = 'pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_planarity = True'
modified_master_phil_str = change_default_phil_values(
  modified_master_phil_str, new_default, phil_parse=iotbx.phil.parse)
print ("modified_master_phil_str:",modified_master_phil_str)
#STOP()

new_default = 'pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_hbonds = True'
modified_master_phil_str = change_default_phil_values(
  modified_master_phil_str, new_default, phil_parse=iotbx.phil.parse)

print ("modified_master_phil_str:",modified_master_phil_str)
#STOP()
'''


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
  resolution                   cryo-EM map resolution in angstrom that needs to be entered by a user
  
  map_weight                   cryo-EM map weight.
                               A user is recommended NOT to specify this, so that it will be automatically optimized.
                               If the map is derived from SAXS, map_weight < 0.3 is recommended so that base pairs of nucleic acids are intact.
  
  start_temperature            If not specified, cryo_fit2 will use the optimized value after automatic exploration between 300 and 1000
  
  final_temperature            (default: 0)
  
  MD_in_each_cycle             Cycle is each iteration of MD from start_temperature to final_temperature.
                               If not specified, cryo_fit2 will use the optimized value after automatic exploration.
  
  nproc                        Number of cores to use for MD parameter exploration.
                               If not specified, cryo_fit2 will use available number of cores/3
  
  number_of_steps              The number of MD steps in each phenix.dynamics
                               If not specified, cryo_fit2 will use the optimized value after automatic exploration.
  
  max_steps_for_final_MD       (default: None)
                               The maximum number of steps in final running of phenix.dynamics.
                               If specified, run up to this number of steps no matter what.
  
  secondary_structure.enabled  (default: True)
                               Most MD simulations tend to break secondary structure. 
                               Therefore, turning on this option is recommended. 
                               If HELIX/SHEET records are present in supplied .pdb file, 
                               automatic search of the existing secondary structures in the given 
                               input pdb file will not be executed.
  
  secondary_structure.protein.remove_outliers
                               (default: True)
                               False may be useful for very poor low-resolution structures by
                               ignoring some hydrogen "bond" if it exceed certain distance threshold
  
  output_dir                   (default: output)
                               output folder name prefix 
  
  keep_origin                  (default: True)
                               If True, write out model with origin in original location.
                               If False, shift origin to (0,0,0). 
  
  progress_on_screen           (default: False)
                               If True, temp= xx dist_moved= xx angles= xx bonds= xx is shown on screen rather than cryo_fit2.log 
                               If False, temp= xx dist_moved= xx angles= xx bonds= xx is NOT shown on screen, and saved into cryo_fit2.log
  
  record_states                (default: False)
                               If True, cryo_fit2 records all states and save it to all_states.pdb (only when cryo_fit2 is successfully completed)
                               However, 3k atoms molecules (like L1 stalk in a ribosome) require more than 160 GB of memory.
                               If False, cryo_fit2 doesn't record states of molecular dynamics.
  
  short                        (default: False)
                               If True, run quickly only to check sanity
'''

  #secondary_structure.nucleic_acid.base_pair.restrain_planarity  (default: True)
  #secondary_structure.nucleic_acid.base_pair.restrain_hbonds  (default: True)
  
  datatypes = ['model', 'real_map', 'phil']
  citations = program_citations
  
  master_phil_str = modified_master_phil_str
  # this is ESSENTIAL to avoid
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
    time_total_start = time.time()
    args = sys.argv[1:]
    
    log = multi_out()
    out=sys.stdout
    log.register("stdout", out)
    
    log_file_name = "cryo_fit2.log"
    logfile = open(log_file_name, "w") # since it is 'w', an existing file will be overwritten. (if this is "a", new info will be appended to an existing file)
    log.register("logfile", logfile)
    logfile.write(str(date_and_time()))

    write_this = "\nPreparing cryo_fit2...\n"
    print (write_this)
    logfile.write(write_this)
    
    print('A user input model: %s' % self.data_manager.get_default_model_name(), file=self.logger)
    model_inp = self.data_manager.get_model()
    
    old_style_RNA, removed_R_prefix_in_RNA_pdb_file_name = remove_R_prefix_in_RNA(self.data_manager.get_default_model_name())
    if (old_style_RNA == True):
      write_this ='''Archaic style of nucleic acids (e.g. RA, RU, RT, RG, RC) were detected in user's pdb file.
phenix can't run with this type of naming.
cryo_fit2 replaced these with A,U,T,G,C and rewrote into ''' + removed_R_prefix_in_RNA_pdb_file_name + '''
Please rerun cryo_fit2 with this re-written pdb file\n'''
      print (write_this)
      logfile.write(write_this)
      logfile.close()
      exit(1)
    
    # if (self.params.sigma_for_custom_geom != None):
    #   user_sigma_for_custom_geom = self.params.sigma_for_custom_geom
    # else:
    #   self.params.sigma_for_custom_geom = 0.04
      
    user_eff_file_provided, user_eff_file_name = check_whether_args_has_eff(args, logfile, "prepare_cryo_fit2", "NA")
    if ((user_eff_file_provided == False) and (self.params.strong_ss == True)):
      
      if (self.params.sigma_for_custom_geom != None):
        user_sigma_for_custom_geom = self.params.sigma_for_custom_geom
      else:
        self.params.sigma_for_custom_geom = 0.04
      
      write_this = "A user didn't provide an .eff file. Therefore, cryo_fit2 will make it automatically to enforce stronger secondary structure restraints.\n"
      print (write_this)
      logfile.write(write_this)
      generated_eff_file_name = write_custom_geometry(logfile, self.data_manager.get_default_model_name(), \
                                                      self.params.sigma_for_custom_geom)
      sys.argv.append(generated_eff_file_name)
    
    logfile.close()
