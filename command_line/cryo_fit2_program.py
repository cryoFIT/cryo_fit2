######to do list
# To run cryo_fit2 at windows as well, replace all unix command like libtbx.easy_run.fully_buffered

from __future__ import division, print_function

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
start_temperature = None
  .type = float
  .short_caption = Starting temperature of annealing in Kelvin. \
                   If not specified, cryo_fit2 will use the optimized value after automatic exploration between 300 and 1000.
final_temperature = 0
  .type = float
  .short_caption = Final temperature of annealing in Kelvin
number_of_steps = 100
  .type = int
  .short_caption = number of steps in phenix.dynamics
number_of_MD_in_each_epoch = 4
  .type = int
  .short_caption = An epoch here is different from the one in deep learning. \
                   Here, the epoch is each iteration of MD from start_temperature to final_temperature.
cool_rate = None
  .type = float
  .short_caption = Cooling rate of annealing in Kelvin. Will be automatically determined by cryo_fit2.
total_number_of_steps = None
  .type = int
  .short_caption = The total number of steps in phenix.dynamics.\
                   If specified, run up to this number of steps no matter what.
map_weight = None
  .type = float
  .short_caption = cryo-EM map weight. \
                   A user is recommended NOT to specify this, so that it will be automatically optimized.
weight_boost = 1
  .type = float
  .short_caption = boost cryo-EM map weight by this much. For a helix, 20 keeps geometry, 100 breaks it.
resolution = None
  .type = float
  .short_caption = cryo-EM map resolution (angstrom) that needs to be specified by a user
output_dir = output
  .type = path
  .short_caption = Output folder PREFIX
progress_on_screen = True
    .type          = bool
    .help          = If True,  temp=x dist_moved=x angles=x bonds=x is shown on screen rather than cryo_fit2.log \
                     If False, temp=x dist_moved=x angles=x bonds=x is NOT shown on screen, but saved into cryo_fit2.log
record_states = False
    .type     = bool
    .help     = If True, cryo_fit2 records all states and save it to all_states.pdb. \
                However, 3k atoms molecules (like L1 stalk in a ribosome) require more than 160 GB of memory. \
                If False, cryo_fit2 doesn't record each state of molecular dynamics.
strong_ss = True
    .type   = bool
    .help   = If True, cryo_fit2 will use a stronger sigma (e.g. 0.021) for secondary structure restraints. \
              If False, it will use the original sigma (e.g. 1)
sigma = 0.021
  .type = float
  .short_caption = The lower this value, the stronger the custom made secondary structure restraints will be. \
                   Oleg recommended 0.021 which is the sigma value for covalent bond. \
                   Doo Nam's benchmark (144 combinations of options) shows that 1.00E-06, 0.1, and 0.2 do not make any difference in base_pair keeping
loose_ss_def = False
    .type   = bool
    .help   = If True, secondary structure definition for nucleic acid is loose. Use this with great caution.  \
              If False, use Oleg's original strict definition. 
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
  
  number_of_MD_in_each_epoch   (default: 4)
                               An epoch here is different from the one in deep learning.
                               Here, the epoch is each iteration of MD from start_temperature to final_temperature.
  
  number_of_steps              (default: 20)
  
  total_number_of_steps        (default: None)
                               If specified, run up to this number of step no matter what.
  
  secondary_structure.enabled  (default: True)
                               Most MD simulations tend to break secondary structure. 
                               Therefore, turning on this option is recommended. 
                               If HELIX/SHEET records are present in supplied .pdb file, 
                               automatic search of the existing secondary structures in the given 
                               input pdb file will not be executed.
  
  secondary_structure.protein.remove_outliers (default: True)
                               False may be useful for very poor low-resolution structures by
                               ignoring some hydrogen "bond" if it exceed certain distance threshold
  
  output_dir                   output folder name prefix (default: output)
  
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
  
  devel                        (default: False)
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


    if (self.params.strong_ss == True):
      write_this = "\nA user turned strong_ss=True\n"
      print (write_this)
      logfile.write(write_this)
      
      eff_file_name = write_custom_geometry(logfile, self.data_manager.get_default_model_name(), self.params.sigma)
      args.append(eff_file_name)
    
    # else:  
    #   if (self.params.sigma != 0.021):
    #     write_this = "\nSpecifying the sigma value when a user turned strong_ss=False is meaningless. \nExit cryo_fit2 now.\n"
    #     print (write_this)
    #     logfile.write(write_this)
    #     exit(1)
    
    print ("user entered resolution", str(self.params.resolution))
    
    print('User input model: %s' % self.data_manager.get_default_model_name(), file=self.logger)
    model_inp = self.data_manager.get_model()
    
    print('User input map: %s' % self.data_manager.get_default_real_map_name(), file=self.logger)
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

    if (self.params.devel == True) :
      self.params.start_temperature = 300
      self.params.final_temperature = 280
      self.params.number_of_MD_in_each_epoch = 2
      self.params.number_of_steps = 1
      self.params.total_number_of_steps = 100

    elif (input_model_file_name_wo_path == "tutorial_cryo_fit2_model.pdb"): 
      self.params.start_temperature = 1000
      self.params.final_temperature = 0
      self.params.number_of_MD_in_each_epoch = 5
      self.params.number_of_steps = 1000
      self.params.total_number_of_steps = 2000

    # regression with total_number_of_steps
    elif (input_model_file_name_wo_path == "tst1_cryo_fit2_model.pdb"): 
      self.params.start_temperature = 300
      self.params.final_temperature = 280
      self.params.number_of_MD_in_each_epoch = 2
      self.params.number_of_steps = 100
      self.params.total_number_of_steps = 1000

    # regression test auto-rerun
    elif (input_model_file_name_wo_path == "tst2_cryo_fit2_model.pdb"): 
      self.params.start_temperature = 300
      self.params.final_temperature = 280
      self.params.number_of_MD_in_each_epoch = 2
      self.params.number_of_steps = 1
    
    
    ########## <begin> Automatic map weight determination
    user_map_weight = ''
    if (self.params.map_weight == None): # a user didn't specify map_weight
      self.params.map_weight = determine_optimal_weight_by_template(self, logfile, map_inp ,'', self.params.weight_boost)
      logfile.write("\nAutomatically optimized map_weight: ")
    else:
      user_map_weight = self.params.map_weight # this user_map_weight will be used later
      logfile.write("\nUser specified map_weight: ")
    
    logfile.write(str(round(self.params.map_weight,1)))
    logfile.write("\n\n")
    ########## <end> Automatic map weight determination
    
    
    '''
    ##################### < begin> explore the optimal combination of parameters
    ######## Based on preliminary benchmarks (~500 combinations with L1 stalk and tRNA), Doonam believes that finding an
    ######## optimum combination of different parameters is a better approach than individually finding each "optimal" parameter
    ## final_temperature is fixed as 0
    ## sigma is fixed as 0.1
    if (self.params.start_temperature == None):
      total_combi_num = 0
      start_temperature_array = []
      for start_temperature in range (300, 901, 300):
        total_combi_num = total_combi_num + 1
        start_temperature_array.append(start_temperature)
      
      number_of_total_cores = know_total_number_of_cores(logfile)
      
      argstuples = [( self, self.params, start_temperature_array[0], logfile), \
                    ( self, self.params, start_temperature_array[1], logfile), \
                    ( self, self.params, start_temperature_array[2], logfile) ]
      for args, res, errstr in easy_mp.multi_core_run( cryo_fit2_by_multi_core, argstuples, number_of_total_cores): # the last argument is nproc
        print ('arguments: %s \n result: %s \n error: %s\n' %(args, res, errstr))
          #print ('arguments: %s \n' %(args))
          #print ('arguments: ', str(args))
    ##################### < end> explore the optimal combination of parameters
    #'''
    
    if (self.params.start_temperature == None):
      self.params.start_temperature = 300
    
    ###############  (begin) core cryo_fit2    
    print ("start_temperature: ", str(self.params.start_temperature))
    print ("final_temperature: ", str(self.params.final_temperature))
    print ("number_of_MD_in_each_epoch: ", str(self.params.number_of_MD_in_each_epoch))
    print ("number_of_steps: ", str(self.params.number_of_steps))
    
    self.params.cool_rate = (self.params.start_temperature-self.params.final_temperature)/(self.params.number_of_MD_in_each_epoch-1)
    print ("cool_rate", str(self.params.cool_rate))
    
    output_dir = get_output_dir_name(self)
    
    ############################# all parameters are determined (either by user or automatic optimization)
    
    cryo_fit2_input_command = "phenix.cryo_fit2 " + self.data_manager.get_default_model_name() + " " \
                            + self.data_manager.get_default_real_map_name() + " " \
                            + "resolution=" + str(self.params.resolution) + " " \
                            + "strong_ss=" + str(self.params.strong_ss) + " " \
                            + "sigma=" + str(self.params.sigma) + " " \
                            + "start_temperature=" + str(self.params.start_temperature) + " " \
                            + "final_temperature=" + str(self.params.final_temperature) + " " \
                            + "number_of_MD_in_each_epoch=" + str(self.params.number_of_MD_in_each_epoch) + " " \
                            + "cool_rate=" + str(round(self.params.cool_rate,1)) + " " \
                            + "number_of_steps=" + str(self.params.number_of_steps) + " " \
                            + "weight_boost=" + str(round(self.params.weight_boost,1)) + " " \
                            + "record_states=" + str(self.params.record_states) + " "
                            #+ "secondary_structure.enabled=" + str(self.params.pdb_interpretation.secondary_structure.enabled) + " " \
                            #+ "secondary_structure.protein.remove_outliers=" + str(self.params.pdb_interpretation.secondary_structure.protein.remove_outliers) + " " \
                            #+ "secondary_structure.nucleic_acid.enabled=" + str(self.params.pdb_interpretation.secondary_structure.nucleic_acid.enabled) + " " \
                            #+ "secondary_structure.nucleic_acid.hbond_distance_cutoff=" + str(self.params.pdb_interpretation.secondary_structure.nucleic_acid.hbond_distance_cutoff) + " " \
                            #+ "secondary_structure.nucleic_acid.angle_between_bond_and_nucleobase_cutoff=" + str(self.params.pdb_interpretation.secondary_structure.nucleic_acid.angle_between_bond_and_nucleobase_cutoff) + " " \
                            #+ "map_weight=" + str(round(self.params.map_weight,1)) + " " \
    if (self.params.total_number_of_steps != None):
      cryo_fit2_input_command = cryo_fit2_input_command + "total_number_of_steps=" + str(self.params.total_number_of_steps) + "\n"
    else:
      cryo_fit2_input_command = cryo_fit2_input_command + "\n"
                              
    print ("\ncryo_fit2_input_command:",cryo_fit2_input_command)
    
    input_command_file = open("cryo_fit2.input_command.txt", "w")
    input_command_file.write(str(cryo_fit2_input_command))
    input_command_file.close()
    
    logfile.write("Input command: ")
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
      weight_boost      = self.params.weight_boost)

    task_obj.validate()
    
    output_dir_final = task_obj.run()
    ############### (end) core cryo_fit2
    
    
    if (self.params.strong_ss == True):
      pymol_ss = input_model_file_name_wo_path + "_ss.pml"
      mv_command_string = "mv " + pymol_ss + " " + eff_file_name + " " + output_dir_final
      libtbx.easy_run.fully_buffered(mv_command_string)
    
    header = "# Geometry restraints used for cryo_fit2\n"
    header += "# %s\n" % date_and_time()
    
    r = model_inp.restraints_as_geo(
        header=header,
        # Stuff for outputting ncs_groups
        #excessive_distance_limit=self.params.ncs.excessive_distance_limit)
        excessive_distance_limit=10)

    # this r is same as the RSR resulted .geo file, therefore I don't need to study about write_geo(m)=True option
    geometry_restraints_file_name = "used_geometry_restraints.geo"
    geo_file = open(geometry_restraints_file_name, "w")
    geo_file.write(r)
    geo_file.close()
    
    # clean up
    mv_command_string = "mv cryo_fit2.input_command.txt " + geometry_restraints_file_name + " " + log_file_name + " " + output_dir_final
    libtbx.easy_run.fully_buffered(mv_command_string)
    
    time_total_end = time.time()
    time_took = show_time(time_total_start, time_total_end)
    
    logfile.write(str("cryo_fit2"))
    logfile.write(str(time_took))
    logfile.write("\n")
    logfile.close()
