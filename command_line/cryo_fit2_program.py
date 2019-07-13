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


########## <begin> import util py files
path = subprocess.check_output(["which", "phenix.cryo_fit2"])
splited_path = path.split("/")
command_path = ''
for i in range(len(splited_path)-3):
  command_path = command_path + splited_path[i] + "/"
command_path = command_path + "modules/cryo_fit2/"
util_path = command_path + "util/"
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
cores_from_user  = None
  .type          = int
  .short_caption = Number of cores to use for MD parameter exploration.
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
MD_in_each_epoch = None
  .type          = int
  .short_caption = An epoch here is different from the one in deep learning. \
                   Here, the epoch is each iteration of MD from start_temperature to final_temperature. \
                   If not specified, cryo_fit2 will use the optimized value by automatic exploration.
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
  .help          = If True, cryo_fit2 records all states and save it to all_states.pdb. \
                   However, 3k atoms molecules (like L1 stalk in a ribosome) require more than 160 GB of memory. \
                   If False, cryo_fit2 doesn't record each state of molecular dynamics.
reoptimize_map_weight_after_each_epoch = False
  .type                                = bool
  .help                                = If True, cryo_fit2 will reoptimize map_weight after each epoch.
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
  .help   = If True, cryo_fit2 will use a stronger sigma_for_custom_geom (e.g. 0.021) for secondary structure restraints. \
            If False, it will use the original sigma_for_custom_geom (e.g. 1)
total_steps      = None
  .type          = int
  .short_caption = The total number of steps in phenix.dynamics.\
                   If specified, run up to this number of steps no matter what.
total_steps_for_exploration  = 5000
  .type                      = int
  .short_caption             = The total number of steps for MD parameter exploration
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
  
  MD_in_each_epoch             An epoch here is different from the one in deep learning.
                               Here, the epoch is each iteration of MD from start_temperature to final_temperature.
                               If not specified, cryo_fit2 will use the optimized value after automatic exploration.
  
  number_of_steps              The number of MD steps in each phenix.dynamics
                               If not specified, cryo_fit2 will use the optimized value after automatic exploration.
  
  total_steps                  (default: None)
                               If specified, run up to this number of step no matter what.
  
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
                               If True, cryo_fit2 records all states and save it to all_states.pdb.
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

    # Importantly declared initial global variables
    user_cool_rate = None
    user_MD_in_each_epoch = None 
    user_number_of_steps = None 
    user_sigma_for_custom_geom = None
    user_start_temperature = None
    user_weight_multiply = None
    
    # save a user entered params.* now
    if (self.params.cool_rate != None):
      user_cool_rate = self.params.cool_rate
    if (self.params.MD_in_each_epoch != None):
      user_MD_in_each_epoch = self.params.MD_in_each_epoch
    if (self.params.number_of_steps != None):
      user_number_of_steps = self.params.number_of_steps
    if (self.params.sigma_for_custom_geom != None):
      user_sigma_for_custom_geom = self.params.sigma_for_custom_geom
    if (self.params.start_temperature != None):
      user_start_temperature = self.params.start_temperature
    if (self.params.weight_multiply != None):
      user_weight_multiply = self.params.weight_multiply

    
    print ("A user entered resolution:", str(self.params.resolution))
    
    print('A user input model: %s' % self.data_manager.get_default_model_name(), file=self.logger)
    model_inp = self.data_manager.get_model()
    
    print('A user input map: %s' % self.data_manager.get_default_real_map_name(), file=self.logger)
    map_inp = self.data_manager.get_real_map()
    
    
    ################# <begin> Doonam's playground ################
    
    ####### works
    #print ("dir(map_inp):",dir(map_inp)) # just shows list of what items are available
    #['__doc__', '__init__', '__module__', 'cannot_be_sharpened', 'crystal_symmetry', 'data', 'external_origin', 'get_additional_labels', 'get_labels', 'get_limitation', 'get_limitations', 'grid_unit_cell', 'header_max', 'header_mean', 'header_min', 'header_rms', 'is_in_limitations', 'is_similar_map', 'labels', 'map_data', 'nxstart_nystart_nzstart', 'origin', 'pixel_sizes', 'show_summary', 'space_group_number', 'statistics', 'unit_cell', 'unit_cell_crystal_symmetry', 'unit_cell_grid', 'unit_cell_parameters']
    
    print ("map_inp.show_summary():", map_inp.show_summary())
    #print ("map_inp.unit_cell_grid():", map_inp.unit_cell_grid()) #TypeError: 'tuple' object is not callable
    #print ("map_inp.unit_cell_grid:", map_inp.unit_cell_grid)
    #print ("map_inp.unit_cell_grid[0]:", map_inp.unit_cell_grid[0]) 

    # for map_boxed map
      # unit cell grid: (360, 360, 360)
      # map grid:   (99, 87, 85)
    # for original map
      # unit cell grid: (360, 360, 360)
      # map grid:   (360, 360, 360)
    # it shows many items from header_min to pixel size, show_summary() itself shows "None"

    '''
    print ("map_inp.unit_cell_crystal_symmetry().unit_cell():",map_inp.unit_cell_crystal_symmetry().unit_cell())
    print ("map_inp.unit_cell_parameters:",map_inp.unit_cell_parameters)
    print ("map_inp.space_group_number:",(map_inp.space_group_number))
    #print ("map_inp.unit_cell_crystal_symmetry():",map_inp.unit_cell_crystal_symmetry()) # just shows the address of the object
    #print ("map_inp.crystal_symmetry():",map_inp.crystal_symmetry()) # just shows the address of the object
    '''
    map_inp_data = map_inp.map_data()
    #print ("map_inp_data:", map_inp_data) #<scitbx_array_family_flex_ext.double object at 0x119e171d0>
    #print ("map origin:", map_inp_data.origin())
    #132, 94, 203 for DN map_L1 stalk same as in util.py's target_map_data.origin()
    
    #'''
    #print ("map accessor:", map_inp_data.accessor()) # shows an object address
    #print ("map_inp.crystal_symmetry():",map_inp.crystal_symmetry()) # shows an object address
    
    ########## test
    #map_inp.space_group_number()
    #print ("map_inp.space_group_number():",map_inp.space_group_number())
    #print ("str(map_inp.space_group_number()):",str(map_inp.space_group_number()))
    
    ########## not works
    '''
    #print ("map_inp.origin():", map_inp.origin()) # ==> TypeError: 'list' object is not callable
    #print ("map_inp.accessor():", map_inp.accessor())
    map_inp.unit_cell_crystal_symmetry()
    print ("map_inp.space_group_number():",map_inp.space_group_number())
    map_inp.unit_cell_parameters()
    print ("map_inp.unit_cell_parameters().unit_cell():",map_inp.unit_cell_parameters().unit_cell())
    '''
    ################# <end> Doonam's playground ################
    
    
    if (self.params.loose_ss_def == True):
      self.params.pdb_interpretation.secondary_structure.nucleic_acid.hbond_distance_cutoff=4
      # default is 3.4, the longer this value, the loose the restraint
      
      self.params.pdb_interpretation.secondary_structure.nucleic_acid.angle_between_bond_and_nucleobase_cutoff=30
      # default is 35, this angle doesn't have much effect on bvht RNA structure
      
    # when I use phenix.secondary_structure_restraints at commandline, use
    # secondary_structure.nucleic_acid.hbond_distance_cutoff
    # not
    # pdb_interpretation.secondary_structure.nucleic_acid.hbond_distance_cutoff
    
    print ("self.params.pdb_interpretation.secondary_structure.enabled:",self.params.pdb_interpretation.secondary_structure.enabled)
    print ("self.params.pdb_interpretation.secondary_structure.protein.remove_outliers:",self.params.pdb_interpretation.secondary_structure.protein.remove_outliers)
    print ("self.params.pdb_interpretation.secondary_structure.nucleic_acid.enabled:",self.params.pdb_interpretation.secondary_structure.nucleic_acid.enabled)
    #print ("self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.enabled:",self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.enabled)
    #print ("self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_planarity:",self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_planarity)
    #print ("self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_hbonds:",self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_hbonds)
    
    old_style_RNA, removed_R_prefix_in_RNA_pdb_file_name = remove_R_prefix_in_RNA(self.data_manager.get_default_model_name())
    if (old_style_RNA == True):
      write_this ='''Archaic style of nucleic acids (e.g. RA, RU, RT, RG, RC) were detected in user's pdb file.
phenix can't run with this type of naming.
cryo_fit2 replaced these with A,U,T,G,C and rewrote into ''' + removed_R_prefix_in_RNA_pdb_file_name + '''
please rerun cryo_fit2 with this re-written pdb file\n'''
      print (write_this)
      logfile.write(write_this)
      logfile.close()
      exit(1)
    
    splited = self.data_manager.get_default_model_name().split("/")
    input_model_file_name_wo_path = splited [len(splited)-1]

    if (self.params.short == True) :
      self.params.start_temperature = 300
      self.params.final_temperature = 280
      self.params.MD_in_each_epoch = 2
      self.params.number_of_steps = 1
      self.params.total_steps = 50

    elif (input_model_file_name_wo_path == "tutorial_cryo_fit2_model.pdb"): 
      self.params.start_temperature = 1000
      self.params.final_temperature = 0
      self.params.MD_in_each_epoch = 5
      self.params.number_of_steps = 1000
      self.params.total_steps = 2000


    ########## <begin> Automatic map weight determination
    user_map_weight = ''
    if (self.params.map_weight == None): # a user didn't specify map_weight
      self.params.map_weight = determine_optimal_weight_by_template(self, logfile, map_inp ,'')
      logfile.write("\nAn automatically optimized map_weight (before any multiplication): ")
    else:
      user_map_weight = self.params.map_weight # this user_map_weight will be used later
      logfile.write("\nA user specified map_weight: ")
    
    logfile.write(str(round(self.params.map_weight,1)))
    logfile.write("\n")
    ########## <end> Automatic map weight determination
    
    
    bp_in_a_user_pdb_file, H_in_a_user_pdb_file, E_in_a_user_pdb_file, ss_file = \
      know_bp_H_E_in_a_user_pdb_file(self.data_manager.get_default_model_name(), logfile)
    
    if ((bp_in_a_user_pdb_file == 0) and (H_in_a_user_pdb_file == 0) and (E_in_a_user_pdb_file == 0)):
        write_this = "A user input file has no base pair or Helix/Sheet.\nMaybe this is an intrinsic molecule. Therefore, cryo_fit2 will not explore MD parameters\n"
        print(write_this)
        logfile.write(write_this)
        self.params.explore == False

    if (self.params.strong_ss == True):
      write_this = "A user turned strong_ss=True\n"
      print (write_this)
      logfile.write(write_this)
      
    
    print ("args outside of make_argstuples fn:",args)
    
    ####################### <begin> Explore the optimal combination of parameters
    if ((self.params.short == False) and (self.params.explore == True)):

      ########  Based on preliminary benchmarks (~500 combinations with L1 stalk and tRNA), Doonam believes that finding an
      ######## optimum combination of different parameters is a better approach than individually finding each "optimal" parameter
      bp_cutoff = bp_in_a_user_pdb_file * 0.97
      write_this = "bp_cutoff from a user input pdb file: " + str(round(bp_cutoff,1)) 
      print(write_this)
      logfile.write(write_this)
      
      # stricter cutoff for H/E than bp
      H_cutoff = H_in_a_user_pdb_file * 0.99
      write_this = "H_cutoff from a user input pdb file: " + str(round(H_cutoff,1))
      print(write_this)
      logfile.write(write_this)
      
      E_cutoff = E_in_a_user_pdb_file * 0.99
      write_this = "E_cutoff from a user input pdb file: " + str(round(E_cutoff,1))
      print(write_this)
      logfile.write(write_this)
        
      if (os.path.isdir("parameters_exploration") == True):
        shutil.rmtree("parameters_exploration")
      os.mkdir("parameters_exploration")
      
      total_combi_num, argstuples = make_argstuples(self, args, logfile, user_map_weight, bp_cutoff, H_cutoff, E_cutoff) # user_map_weight should tag along for a later usage
      
      cores_to_use = ''
      if (self.params.cores_from_user != None):
        cores_to_use = self.params.cores_from_user
      else:
        returned_nproc = number_of_processors(return_value_if_unknown=-1)
        # kaguya resulted in 32
        # sparky resulted in 40 (I expected to see 34 since I was running 6 cores at that time. \
        #It seems that number_of_processors returned just all # of processors)
        
        cores_to_use = math.ceil(returned_nproc/4) # will use at least 1 core, since math.ceil rounds up to the next greater integer
        # just to avoid crash, it seems like sparky linux machine can't handle more than 40 cores (even 20 cores)
        # when I used 13 cores, the load average reached 20!
      
      write_this = "Cryo_fit2 will use " + str(int(cores_to_use)) + " core(s) to explore up to " + str(total_combi_num) + " MD parameters.\n"
      print(write_this)
      logfile.write(write_this)
      
      success_exploration_count = 0
      for arguments_for_explore, res, errstr in easy_mp.multi_core_run( explore_parameters_by_multi_core, argstuples, \
                                                       cores_to_use): # the last argument is nproc
          print ("explore_parameters_by_multi_core ran")
          
          #print ('arguments_for_explore: %s ' %(arguments_for_explore)) # "arguments_for_explore:   [<cryo_fit2_program.Program object at 0x11b54ffd0>, <libtbx.phil.scope_extract object at 0x11b54fd50>, <open file 'cryo_fit2.log', mode 'w' at 0x11b45b9c0>, '', 0.0, 0.99, 0.0, 2, 1, 0.001, 300, 1]"
          #print ('arguments_for_explore: ', str(arguments_for_explore)) # "arguments_for_explore:  (<cryo_fit2_program.Program object at 0x10dc66910>, <libtbx.phil.scope_extract object at 0x10dc66b90>, 900, <open file 'cryo_fit2.log', mode 'w' at 0x10dceba50>, '')"
          
          if (res != None):
            write_this = 'bp:' + str(res[0]) + ', H:' + str(res[1]) + ', E:' + str(res[2])
            print (write_this) # 1, 0, 0
            logfile.write(write_this)
          
          write_this = 'error string: %s ' %(errstr) + '\n'
          
          # -> this errstr will be either "None" or
          '''/Users/builder/slave/phenix-nightly-mac-intel-osx-x86_64/modules/cctbx_project/cctbx/xray/sampling_base.h: expone\
nt_table: excessive range.
Traceback (most recent call last):
  File "/Users/doonam/bin/phenix-1.15rc3-3442/modules/cctbx_project/libtbx/scheduling/job_scheduler.py", line 64, in job_\
cycle
    value = target( *args, **kwargs )
  File "/Users/doonam/bin/phenix-1.15rc3-3442/modules/cryo_fit2/util/util.py", line 237, in explore_parameters_by_multi_c\
ore
    output_dir_final = task_obj.run()
  File "/Users/doonam/bin/phenix-1.15rc3-3442/modules/cryo_fit2/command_line/cryo_fit2_run.py", line 170, in run
    cc_after_small_MD = calculate_cc(map_data=map_data, model=self.model, resolution=self.params.resolution)
  File "/Users/doonam/bin/phenix-1.15rc3-3442/modules/cryo_fit2/util/util.py", line 26, in calculate_cc
    fc = xrs.structure_factors(d_min = resolution).f_calc()
  File "/Users/doonam/bin/phenix-1.15rc3-3442/modules/cctbx_project/cctbx/xray/structure.py", line 1573, in structure_fac\
tors
    algorithm=algorithm)
  File "/Users/doonam/bin/phenix-1.15rc3-3442/modules/cctbx_project/cctbx/xray/structure_factors/from_scatterers.py", lin\
e 53, in __call__
    algorithm=algorithm) # passing algorithm allows f to decide on CPU/GPU implementation
  File "/Users/doonam/bin/phenix-1.15rc3-3442/modules/cctbx_project/cctbx/xray/structure_factors/from_scatterers_fft.py",\
 line 38, in __init__
    tolerance_positive_definite=manager.tolerance_positive_definite())""
    '''

          print (write_this)
          logfile.write(write_this)
          
          if (errstr == None):
            success_exploration_count = success_exploration_count + 1
      
      write_this = "\ncryo_fit2 explored " + str(success_exploration_count) + " combination(s) of MD parameters " + \
                   "out of " + str(total_combi_num) + " total combinations.\nIt will run fully with optimized parameters.\n"
      
      print (write_this)
      logfile.write(write_this)

      optimum_MD_in_each_epoch, optimum_sigma_for_custom_geom, optimum_start_temperature, optimum_steps,  \
      optimum_weight_multiply = extract_the_best_cc_parameters(logfile)
      
      self.params.MD_in_each_epoch      = int(optimum_MD_in_each_epoch)
      self.params.sigma_for_custom_geom = float(optimum_sigma_for_custom_geom)
      self.params.start_temperature     = float(optimum_start_temperature) # make it as float to format it consistent as in parameter exploration and user input
      self.params.number_of_steps       = int(optimum_steps)
      self.params.weight_multiply       = float(optimum_weight_multiply)

      # Override self.params.* with user entered values
      if (user_MD_in_each_epoch != None):
        self.params.MD_in_each_epoch = user_MD_in_each_epoch
      if (user_number_of_steps != None):
        self.params.number_of_steps = user_number_of_steps
      if (user_sigma_for_custom_geom != None):
        self.params.sigma_for_custom_geom = user_sigma_for_custom_geom
      if (user_start_temperature != None):
        self.params.start_temperature = user_start_temperature
      if (user_weight_multiply != None):
        self.params.weight_multiply = user_weight_multiply
    ####################### <end> explore the optimal combination of parameters
  
  
    ### Assign default values if not specified till now
    if (self.params.MD_in_each_epoch == None):
      self.params.MD_in_each_epoch = 4
    if (self.params.number_of_steps == None):
      self.params.number_of_steps = 100
    if (self.params.sigma_for_custom_geom == None):
      self.params.sigma_for_custom_geom = 0.021
    if (self.params.start_temperature == None):
      self.params.start_temperature = 300
    if (self.params.weight_multiply == None):
      self.params.weight_multiply = 1


    ###############  (begin) core cryo_fit2    
    print ("Final MD parameters after user input/automatic optimization")
    print ("final_temperature     :", str(self.params.final_temperature))
    print ("MD_in_each_epoch      :", str(self.params.MD_in_each_epoch))
    print ("number_of_steps       :", str(self.params.number_of_steps))
    print ("sigma_for_custom_geom :", str(self.params.sigma_for_custom_geom))
    print ("start_temperature     :", str(self.params.start_temperature))
        
    # override self.params.* with user entered values
    if (user_cool_rate != None):
        self.params.cool_rate = user_cool_rate
    else:
      self.params.cool_rate = (self.params.start_temperature-self.params.final_temperature)/(self.params.MD_in_each_epoch-1)
    print ("cool_rate             :", str(round(self.params.cool_rate,1)), "\n")
    
    if (self.params.sigma_for_custom_geom == None):
      self.params.sigma_for_custom_geom = 0.021
    
    eff_file_provided = check_whether_args_has_eff(args, logfile)
    
    
    ''' # old works
    eff_file_provided = False
    for i in range(len(args)):
      if ((args[i][(len(args[i])-4):len(args[i])]) == ".eff"):
        user_eff_file_name = str(args[i])
        write_this = "A user provided an .eff file (e.g. " + user_eff_file_name + "), cryo_fit2 will use it."
        print (write_this)
        logfile.write(write_this)
        eff_file_provided = True
    '''
    
    if ((eff_file_provided == False) and (self.params.strong_ss == True)): # If optimal sigma is not found (or exploration is not tried in the first place)
      write_this = "A user didn't provide an .eff file. Therefore, cryo_fit2 will make it automatically to enforce stronger secondary structure restraints.\n"
      print (write_this)
      logfile.write(write_this)
      eff_file_name = write_custom_geometry(logfile, self.data_manager.get_default_model_name(), self.params.sigma_for_custom_geom)
      args.append(eff_file_name)
    
    
    output_dir = get_output_dir_name(self)
    
    # All parameters are determined (either by a user or automatic optimization)    
    cryo_fit2_input_command = "phenix.cryo_fit2 " + self.data_manager.get_default_model_name() \
                            + " " + self.data_manager.get_default_real_map_name()  \
                            + " resolution=" + str(self.params.resolution)  \
                            + " strong_ss=" + str(self.params.strong_ss) \
                            + " sigma_for_custom_geom=" + str(self.params.sigma_for_custom_geom) \
                            + " start_temperature=" + str(self.params.start_temperature)  \
                            + " final_temperature=" + str(self.params.final_temperature) \
                            + " MD_in_each_epoch=" + str(self.params.MD_in_each_epoch) \
                            + " cool_rate=" + str(round(self.params.cool_rate,1)) \
                            + " number_of_steps=" + str(self.params.number_of_steps) \
                            + " weight_multiply=" + str(round(self.params.weight_multiply,1)) \
                            + " record_states=" + str(self.params.record_states) \
                            + " reoptimize_map_weight_after_each_epoch=" + str(self.params.reoptimize_map_weight_after_each_epoch) \
                            + " explore=False"
                            #+ "secondary_structure.enabled=" + str(self.params.pdb_interpretation.secondary_structure.enabled) + " " \
                            #+ "secondary_structure.protein.remove_outliers=" + str(self.params.pdb_interpretation.secondary_structure.protein.remove_outliers) + " " \
                            #+ "secondary_structure.nucleic_acid.enabled=" + str(self.params.pdb_interpretation.secondary_structure.nucleic_acid.enabled) + " " \
                            #+ "secondary_structure.nucleic_acid.hbond_distance_cutoff=" + str(self.params.pdb_interpretation.secondary_structure.nucleic_acid.hbond_distance_cutoff) + " " \
                            #+ "secondary_structure.nucleic_acid.angle_between_bond_and_nucleobase_cutoff=" + str(self.params.pdb_interpretation.secondary_structure.nucleic_acid.angle_between_bond_and_nucleobase_cutoff) + " " \
                            #+ "map_weight=" + str(round(self.params.map_weight,1)) + " " \
    if (eff_file_provided == True):
      cryo_fit2_input_command = cryo_fit2_input_command + " " + user_eff_file_name
    if (self.params.total_steps != None):
      cryo_fit2_input_command = cryo_fit2_input_command + " total_steps=" + str(self.params.total_steps)
    cryo_fit2_input_command = cryo_fit2_input_command + "\n"
    
    print ("\ncryo_fit2_input_command:",cryo_fit2_input_command)
    
    input_command_file = open("cryo_fit2.input_command.txt", "w")
    input_command_file.write(str(cryo_fit2_input_command))
    input_command_file.close()
    
    logfile.write("\nAn input command for final cryo_fit2 MD run:\n")
    logfile.write(str(cryo_fit2_input_command))
    
    task_obj = cryo_fit2_run.cryo_fit2_class(
      model             = model_inp,
      model_name        = self.data_manager.get_default_model_name(),
      map_inp           = map_inp,
      params            = self.params,
      out               = self.logger,
      map_name          = self.data_manager.get_default_real_map_name(),
      logfile           = logfile,
      output_dir        = output_dir,
      user_map_weight   = user_map_weight,
      weight_multiply   = self.params.weight_multiply)

    task_obj.validate()
    
    output_dir_final = task_obj.run()
    ############### (end) core cryo_fit2
    
    if (self.params.strong_ss == True):
      pymol_ss = input_model_file_name_wo_path + "_ss.pml"
      mv_command_string = "mv " + pymol_ss + " " + eff_file_name + " " + output_dir_final
      libtbx.easy_run.fully_buffered(mv_command_string)
    
    write_geo(self, model_inp, "used_geometry_restraints.geo")
  
    # clean up
    mv_command_string = "mv cryo_fit2.input_command.txt " + ss_file + " used_geometry_restraints.geo " + log_file_name + " " + output_dir_final
    libtbx.easy_run.fully_buffered(mv_command_string)
    
    time_total_end = time.time()
    time_took = show_time("cryo_fit2", time_total_start, time_total_end)
    print (time_took)
    logfile.write(str(time_took))
    logfile.write("\n")
    
    logfile.close()
