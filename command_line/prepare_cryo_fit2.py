# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_SIGNALS_DEFAULT=1
#### (source) http://cci.lbl.gov/cctbx_sources/crys3d/command_line/hklview.py

######to do list
# To run cryo_fit2 at windows as well, replace all unix command like libtbx.easy_run.fully_buffered
from __future__ import division, print_function

import math, os, shutil, subprocess, sys, time

import cryo_fit2_run

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
cc_only           = False
  .type           = bool
  .help           = If True, cryo_fit2 will calculate cc (overall) only, not running cryo_fit2 itself.
cool_rate        = None
  .type          = float
  .short_caption = Cooling rate of annealing in Kelvin. Will be automatically determined by cryo_fit2.
explore          = False
  .type          = bool
  .short_caption = If True, cryo_fit2 will use maximum number of multiple cores to explore the most optimal MD parameters.\
                   However, this exploration requires a lot of computing power (e.g. > 128 GB memory, > 20 cores).\
                   Exploring at a macbook pro (16 GB memory, 2 cores) crashed.
final_temperature = 0
  .type           = float
  .short_caption  = Final temperature of annealing in Kelvin
HE_angle_sigma_scale = 1
  .type              = float
  .short_caption     = Multiply sigmas for h-bond angles by this value. Original sigmas range from 5 to 10.\
                       As Doo Nam decreased this value, cryo_fit2 takes longer. Therefore, lower value may strengthen restraint.
HE_sigma            = 0.05
  .type             = float
  .short_caption    = The lower this value, the stronger the custom made secondary structure restraints will be. \
                      Oleg once recommended 0.021 which is the sigma value for covalent bond. \
                      According to a small benchmark with a RNA molecule (e.g. L1 stalk), 0.05 best preserves the number of base-pairs.
HE_slack            = 0
  .type             = float
  .short_caption    = As Doo Nam understands /modules/cctbx_project/mmtbx/monomer_library/pdb_interpretation.py, \
                      default value is 0. Indeed, Oleg confirmed that slack should be always 0 for proper geometry restraints. (~Sep, 2019)\
                      It is true that 3.5 Angstrom is a usual width with Go-model. \
                      Therefore, Doo Nam tempted to try 1.7 slack to allow more flexible equilibrium. \
                      However, Doo Nam confirmed that 3.5 HE_slack not only broke helix geometry but also hinders cc improvement in both benchmark cases.\
                      Similarly, selectively disabling helix secondary structure restraints not only broke helix geometry but also hinders cc improvement.
HE_top_out          = False
  .type             = bool
  .help             = If True, top_out potential is used rather than harmonic potential for helix and sheets.
# ignore_symmetry_conflicts = True
#   .type                   = bool
#   .help                   = You can ignore the symmetry information (CRYST1) from \
#                             coordinate files. This may be necessary if your model has been\
#                             placed in a box with box_map for example.
#   .expert_level           = 2
keep_origin      = True
  .type          = bool
  .help          = If True, write out model with origin in original location.  \
                   If False, shift map origin to (0,0,0).
  .short_caption = Keep origin of a resulted atomistic model
loose_ss_def = False
  .type      = bool
  .help      = If True, secondary structure definition for nucleic acid is loose. Use this with great caution.  \
               If False, use Oleg's original strict definition.
map_weight       = None
  .type          = float
  .short_caption = cryo-EM map weight. \
                   A user is recommended NOT to specify this, so that it will be automatically optimized.
map_weight_multiply  = None
  .type              = float
  .short_caption     = Cryo_fit2 will multiply cryo-EM map weight by this much. \ 
                       If not specified, cryo_fit2 will use the default value (e.g. 1) \
                       Usually, a small molecule (a helix) requires only 1 (not multiply). \
                       For a helix, 20 keeps geometry, 100 breaks it (w/o special sigma) \
                       However, a large molecule needs a larger value (e.g. 10~50).
max_steps_for_exploration = 10000
  .type                   = int
  .short_caption          = The total number of steps for MD parameter exploration. \
                            10k is enough to discern Mg Channel \
                            15k is not enough for tRNA
max_steps_for_final_MD = None
  .type                = int
  .short_caption       = The maximum number of steps in final running of phenix.dynamics.\
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
parallelity_sigma = 0.0335
  .type           = float
  .help           = 0.0335 is default by Oleg
planarity_sigma   = 0.176
  .type           = float
  .help           = 0.176 is default by Oleg
progress_on_screen = True
  .type            = bool
  .help            = If True,  temp=x dist_moved=x angles=x bonds=x is shown on screen rather than cryo_fit2.log \
                     If False, temp=x dist_moved=x angles=x bonds=x is NOT shown on screen, but saved into cryo_fit2.log
record_states    = True
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
stacking_pair_sigma = 0.027
  .type             = float
  .help             = 0.027 is default by Oleg
start_temperature = None
  .type           = float
  .short_caption  = Starting temperature of annealing in Kelvin. \
                    If not specified, cryo_fit2 will 300.\
                    (ref.) modules/cctbx_project/mmtbx/dynamics/simulated_anealing.py uses 5k Kelvin starting temperature
stronger_ss = False
  .type     = bool
  .help     = If True, cryo_fit2 will use a stronger HE_sigma for secondary structure restraints. \
              If False, it will not use custom geometry
write_custom_geom_only   = False
  .type                  = bool
  .help                  = If True, cryo_fit2 will write custom geometry eff file only so that users can modify it for user's need, and provide to cryo_fit2 (cryo_fit2 itself will not run).
include scope mmtbx.maps.map_model_cc.master_params
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


class Program(ProgramTemplate):

  description = '''
Program for running cryo_fit2.\n

Minimum required inputs:
  Model file     (.pdb or .cif)
  Map file       (MRC/ccp4 format)
  Map resolution (in angstrom)

Example running command:
  phenix.cryo_fit2 model.pdb map.ccp4 resolution=5

Options:
  cc_only                      (default: False)
                               If True, calculates cc only without running cryo_fit2 itself.
  
  final_temperature            (default: 0)
                               final_temperature of phenix.dynamics.
                               
  HE_angle_sigma_scale         (default: 1)
                               Multiply sigmas for hydrogen bond angles by this value. Original sigmas range from 5 to 10.
                               
  HE_sigma                     (default: 0.05)
                               The lower this value, the stronger the custom made secondary structure restraints will be.
                               Oleg once recommended 0.021 which is the sigma value for covalent bond.
                               According to a small benchmark with a RNA molecule (e.g. L1 stalk), 0.05 best preserves number of base-pairs.
  
  HE_slack                     (default: 0)
                               Other than developers, it is recommended not to change this default value.
  
  keep_origin                  (default: True)
                               If True, write out model with origin in original location.
                               If False, shift origin to (0,0,0). 
  
  max_steps_for_final_MD       (default: None)
                               The maximum number of steps in final running of phenix.dynamics.
                               If specified, run up to this number of steps no matter what.
                               
  map_weight                   cryo-EM map weight.
                               A user is recommended NOT to specify this, so that it will be automatically optimized.
                               If the map is derived from SAXS, map_weight < 0.3 is recommended so that base pairs of nucleic acids are intact.
                
  MD_in_each_cycle             Cycle is each iteration of MD from start_temperature to final_temperature.
                               If not specified, cryo_fit2 will use the optimized value after automatic exploration.
  
  nproc                        Number of cores to use for MD parameter exploration.
                               If not specified, cryo_fit2 will use available number of cores/3
                               
  number_of_steps              The number of MD steps in each phenix.dynamics
                               If not specified, cryo_fit2 will use the optimized value after automatic exploration.
                               
  output_dir                   (default: output)
                               output folder name prefix
                               
  progress_on_screen           (default: False)
                               If True, temp= xx dist_moved= xx angles= xx bonds= xx is shown on screen rather than cryo_fit2.log 
                               If False, temp= xx dist_moved= xx angles= xx bonds= xx is NOT shown on screen, and saved into cryo_fit2.log
                               
  record_states                (default: False)
                               If True, cryo_fit2 records all states and save it to all_states.pdb (only when cryo_fit2 is successfully completed)
                               However, 3k atoms molecules (like L1 stalk in a ribosome) require more than 160 GB of memory.
                               If False, cryo_fit2 doesn't record states of molecular dynamics.
                               
  resolution                   cryo-EM map resolution in angstrom that needs to be entered by a user.
  
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
  
  short                        (default: False)
                               If True, run quickly only to check sanity.
                               
  start_temperature            (default: 300)
                               start_temperature of phenix.dynamics.
  
  write_custom_geom_only       (default: False)
                               If True, write custom geometry eff file only (not running cryo_fit2).
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


    print('A user input model name:    %s' % self.data_manager.get_default_model_name(), file=self.logger)
    model_inp = self.data_manager.get_model()

    print('A user input map file name: %s' % self.data_manager.get_default_real_map_name(), file=self.logger)
    map_inp = self.data_manager.get_real_map()
    
    
    ######### <begin> calculates CC_overall before cryo_fit2
    add_CRYST1_to_pdb_file(self, logfile, map_inp, self.data_manager.get_default_model_name())
    
    try:
      cc_before_cryo_fit2 = round(calculate_overall_cc(map_data=map_inp.map_data(), model=model_inp, resolution=self.params.resolution), 4)
      # Pavel thinks that cc_box should be pretty much similar as this cc_before_cryo_fit2
      #Doo Nam confirmed that even without CRYST1 in pdb file, cc value here does not match cc_box (but close)
      
      write_this = "\n\nCC_overall before cryo_fit2 (both exploration and final MD): " + str(cc_before_cryo_fit2) + "\n"
      print('%s' %(write_this))
      logfile.write(str(write_this))
    except:
      write_error_message_for_overall_cc(logfile)
      
      
    if (self.params.cc_only == True):
      crystal_symmetry = mmtbx.utils.check_and_set_crystal_symmetry(
        models = [model_inp], map_inps=[map_inp])
      report_map_model_cc(self, map_inp, model_inp, crystal_symmetry, logfile)
      logfile.close()
      exit(1)
    ######### <end> calculates CC_overall before cryo_fit2
    
    
    write_this = "\nPreparing cryo_fit2...\n"
    print (write_this)
    logfile.write(write_this)
    
    
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
    
    cleaned_pdb_file_name, cleaned_unusual_residue = clean_unusual_residue (self.data_manager.get_default_model_name())
    if (cleaned_unusual_residue == True):
      write_this ='''
       Unusual residue names like 34G that real_space_refine can't deal were detected in user's pdb file.\nCryo_fit2 removed these and rewrote into ''' + cleaned_pdb_file_name + '''\nPlease rerun cryo_fit2 with this re-written pdb file\n'''
      print (write_this)
      logfile.write(write_this)
      logfile.close()
      exit(1)
    
    leave_one_conformer(logfile, self.data_manager.get_default_model_name())
    
    
    ############# (begin) Assign sigma/slack for H/E
    if ((self.params.HE_sigma != 0.05) or (self.params.HE_slack != 0.0) or (self.params.HE_angle_sigma_scale != 1) or (self.params.HE_top_out == True)):
      generated_eff_file_name = assign_ss_params_to_H_E(logfile, self.data_manager.get_default_model_name(), \
                                                      self.params.HE_sigma, self.params.HE_slack, self.params.HE_angle_sigma_scale,\
                                                      self.params.HE_top_out)
      if (generated_eff_file_name != False):
        sys.argv.append(generated_eff_file_name)
    ############# (end) Assign sigma/slack for H/E
    
    
    ############# (begin) Assign sigmas for nucleic_acids
    if ((self.params.parallelity_sigma != 0.0335) or (self.params.planarity_sigma != 0.176) or (self.params.stacking_pair_sigma != 0.027)):
      generated_eff_file_name_w_nucleic_acid_sigmas = assign_nucleic_acid_sigmas(logfile, \
                                                  self.data_manager.get_default_model_name(), self.params.parallelity_sigma, self.params.planarity_sigma, self.params.stacking_pair_sigma)
      if (generated_eff_file_name_w_nucleic_acid_sigmas != False):
        sys.argv.append(generated_eff_file_name_w_nucleic_acid_sigmas)
    ############# (end) Assign sigmas for nucleic_acids
    
    
    if (self.params.write_custom_geom_only == True):
      generated_eff_file_name = write_custom_geometry(logfile, self.data_manager.get_default_model_name(), \
                                                      self.params.stronger_ss_sigma, self.params.stronger_ss_slack)
      exit(1)
      
      
    '''
    # seems wrong to assign sigma and slack ? keep for now
    ############# (begin) deal with Doonam's stronger_ss
    if (self.params.stronger_ss == True):
      
      generated_eff_file_name = write_custom_geometry(logfile, self.data_manager.get_default_model_name(), \
                                                      self.params.stronger_ss_sigma, self.params.stronger_ss_slack)
      
      sys.argv.append(generated_eff_file_name)
    ############# (end) deal with Doonam's stronger_ss
    '''
    
    logfile.close()
