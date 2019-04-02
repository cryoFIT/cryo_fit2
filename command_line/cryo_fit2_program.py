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

from mmtbx.refinement.real_space import weight

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
#print ("util_path:",util_path)
from util import *


program_citations = libtbx.phil.parse('''
citation {
  article_id = Pavel_s_2018_realspace_refine mentioned and used this phenix.dynamics
}
''')


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
    args = sys.argv[1:]
    #checked_whether_args_has_eff = check_whether_args_has_eff(args)
    
    print ("user entered resolution", str(self.params.resolution))
    print ("start_temperature", str(self.params.start_temperature))
    print ("final_temperature", str(self.params.final_temperature))
    print ("cool_rate", str(self.params.cool_rate))
    print ("number_of_steps", str(self.params.number_of_steps)) 
    
    print ("self.params.loose_ss_def:",self.params.loose_ss_def)

    time_total_start = time.time()
    
    input_model_file_name = self.data_manager.get_default_model_name()
    print('User input model: %s' % input_model_file_name, file=self.logger)
    model_inp = self.data_manager.get_model()
    
    print('User input map: %s' % self.data_manager.get_default_real_map_name(), file=self.logger)
    map_inp = self.data_manager.get_real_map()
    map_inp.crystal_symmetry()
    
    ####### works
    print ("dir(map_inp):",dir(map_inp))
    #print ("map_inp.show_summary():", map_inp.show_summary())
    #STOP()
    print ("map_inp.unit_cell_crystal_symmetry().unit_cell():",map_inp.unit_cell_crystal_symmetry().unit_cell())
    print ("map_inp.unit_cell_parameters:",map_inp.unit_cell_parameters)
    #print ("str(map_inp.space_group_number):",str(map_inp.space_group_number))
    print ("map_inp.space_group_number:",(map_inp.space_group_number))
    #STOP()
    
    # just shows the address of the object
    #print ("map_inp.unit_cell_crystal_symmetry():",map_inp.unit_cell_crystal_symmetry())
    #print ("map_inp.crystal_symmetry():",map_inp.crystal_symmetry())
    
    # just shows address of the instance
    #print ("map_inp.crystal_symmetry:",map_inp.crystal_symmetry)
    
    
    ########## test
    #map_inp.space_group_number()
    #print ("map_inp.space_group_number():",map_inp.space_group_number())
    #print ("str(map_inp.space_group_number()):",str(map_inp.space_group_number()))
    
    
    ########## not works
    '''
    map_inp.unit_cell_crystal_symmetry()
    print ("map_inp.space_group_number():",map_inp.space_group_number())
    map_inp.unit_cell_parameters()
    print ("map_inp.unit_cell_parameters().unit_cell():",map_inp.unit_cell_parameters().unit_cell())
    '''
    
    '''
    has_nucleic_acid = check_whether_the_pdb_file_has_nucleic_acid(self.data_manager.get_default_model_name())
    if (has_nucleic_acid == True):
      if (self.params.map_weight == None):
        self.params.map_weight = 0.4
      #elif (self.params.map_weight > 0.4): # for 20A saxs derived "cryo-EM" map
      #  self.params.map_weight = 0.4
      elif (self.params.map_weight > 5): # for 30A saxs derived "cryo-EM" map
        self.params.map_weight = 5
    '''
    
    if (self.params.loose_ss_def == True):
      self.params.pdb_interpretation.secondary_structure.nucleic_acid.hbond_distance_cutoff=4
      self.params.pdb_interpretation.secondary_structure.nucleic_acid.angle_between_bond_and_nucleobase_cutoff=30
    
    
    print ("self.params.pdb_interpretation.secondary_structure.enabled:",self.params.pdb_interpretation.secondary_structure.enabled)
    print ("self.params.pdb_interpretation.secondary_structure.protein.remove_outliers:",self.params.pdb_interpretation.secondary_structure.protein.remove_outliers)
    print ("self.params.pdb_interpretation.secondary_structure.nucleic_acid.enabled:",self.params.pdb_interpretation.secondary_structure.nucleic_acid.enabled)
    #print ("self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.enabled:",self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.enabled)
    
    #print ("self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_planarity:",self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_planarity)
    #STOP()
    #print ("self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_hbonds:",self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_hbonds)
    
    log = multi_out()
    out=sys.stdout
    log.register("stdout", out)
    
    splited = input_model_file_name.split("/")
    input_model_file_name_wo_path = splited [len(splited)-1]

    if ((self.params.devel == True) or (input_model_file_name_wo_path == "devel_cryo_fit2_model.pdb") or \
          (input_model_file_name_wo_path == "devel_cryo_fit2_model.cif")):
      # "tst..." lives in modules/cryo_fit2/regression
      self.params.start_temperature = 300
      self.params.final_temperature = 280
      self.params.cool_rate = 10
      self.params.number_of_steps = 1

    if (input_model_file_name_wo_path == "tutorial_cryo_fit2_model.pdb"): 
      self.params.start_temperature = 1000
      self.params.final_temperature = 0
      self.params.cool_rate = 10
      self.params.number_of_steps = 1000
      self.params.pdb_interpretation.secondary_structure.enabled = True
      
    #"_map_wt_" + str(round(self.params.map_weight,1)) + \
    # rename output_dir
    output_dir_prefix = self.params.output_dir
    output_dir = str(output_dir_prefix) + \
                 "_resolution_" + str(self.params.resolution) + \
                 "_start_" + str(self.params.start_temperature) + \
                 "_final_" + str(self.params.final_temperature) + \
                 "_cool_" + str(self.params.cool_rate) + \
                 "_step_" + str(self.params.number_of_steps) #+ \
                 #"_ss_" + str(self.params.pdb_interpretation.secondary_structure.enabled) + \
                 #"_del_outlier_ss_" + str(self.params.pdb_interpretation.secondary_structure.protein.remove_outliers) + \
                 #"_NA_" + str(self.params.pdb_interpretation.secondary_structure.nucleic_acid.enabled) + \
                 #"_hb_dis_" + str(self.params.pdb_interpretation.secondary_structure.nucleic_acid.hbond_distance_cutoff) + \
                 #"_angle_" + str(self.params.pdb_interpretation.secondary_structure.nucleic_acid.angle_between_bond_and_nucleobase_cutoff)
                 #"_bp_planar_" + str(self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_planarity) + \
                 #"_bp_hb_" + str(self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_hbonds)
    
    log_file_name = "cryo_fit2.log"
    
    logfile = open(log_file_name, "w") # since it is 'w', an existing file with the same name will be erased
    #logfile = open(log_file_name, "a") # since it is 'a', new info will be appended to an existing file
    log.register("logfile", logfile)
    
    
    ###############  (begin) when optimizing map_weight once
    if (self.params.map_weight == None): # a user didn't specify map_weight
        self.params.map_weight = determine_optimal_weight_by_template(self, map_inp)
        logfile.write("\nAutomatically optimized map_weight: ")
    else:
      logfile.write("\nUser specified map_weight: ")
    
    logfile.write(str(round(self.params.map_weight,1)))
    logfile.write("\n\n")
                         
    #if (checked_whether_args_has_eff == False):   
    cryo_fit2_input_command = "phenix.cryo_fit2 " + self.data_manager.get_default_model_name() + " " + self.data_manager.get_default_real_map_name() + " " \
                            + "resolution=" + str(self.params.resolution) + " " \
                            + "map_weight=" + str(round(self.params.map_weight,1)) + " " \
                            + "start_temperature=" + str(self.params.start_temperature) + " " \
                            + "final_temperature=" + str(self.params.final_temperature) + " " \
                            + "cool_rate=" + str(self.params.cool_rate) + " " \
                            + "steps=" + str(self.params.number_of_steps) + " " \
                            + "secondary_structure.enabled=" + str(self.params.pdb_interpretation.secondary_structure.enabled) + " " \
                            + "secondary_structure.protein.remove_outliers=" + str(self.params.pdb_interpretation.secondary_structure.protein.remove_outliers) + " " \
                            + "secondary_structure.nucleic_acid.enabled=" + str(self.params.pdb_interpretation.secondary_structure.nucleic_acid.enabled) + " " \
                            + "secondary_structure.nucleic_acid.hbond_distance_cutoff=" + str(self.params.pdb_interpretation.secondary_structure.nucleic_acid.hbond_distance_cutoff) + " " \
                            + "secondary_structure.nucleic_acid.angle_between_bond_and_nucleobase_cutoff=" + str(self.params.pdb_interpretation.secondary_structure.nucleic_acid.angle_between_bond_and_nucleobase_cutoff) + " " \
                            + "\n"
    print ("cryo_fit2_input_command:",cryo_fit2_input_command)
    
    input_command_file = open("cryo_fit2.input_command.txt", "w")
    input_command_file.write(str(cryo_fit2_input_command))
    input_command_file.close()
    
    logfile.write("Input command: ")
    logfile.write(str(cryo_fit2_input_command))
    
    '''
    if (checked_whether_args_has_eff == True):
      logfile.write("User entered custom geometry restraints already.\n")
    else:
      logfile.write("User did not enter custom geometry restraints, make it now.\n")
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
  
      # splited_input_model_file_name = input_model_file_name.split("/")
      # input_model_file_name_wo_path = splited_input_model_file_name[len(splited_input_model_file_name)-1]
      ss_restraints_file_name = input_model_file_name_wo_path + "_ss.pml"
      rewrite_to_custom_geometry(ss_restraints_file_name)
      custom_geom_file_name = ss_restraints_file_name[:-4] + "_custom_geom.eff"
    '''
    task_obj = cryo_fit2_run.cryo_fit2_class(
      model             = model_inp,
      model_name        = self.data_manager.get_default_model_name(),
      map_inp           = map_inp,
      params            = self.params,
      out               = self.logger,
      map_name          = self.data_manager.get_default_real_map_name(),
      logfile           = logfile,
      output_dir        = output_dir)
    
    task_obj.validate()
    
    output_dir_w_CC = task_obj.run()
    #STOP()
    ############### (end) when optimizing map_weight once
    
    '''
    if (checked_whether_args_has_eff == False):
      mv_command_string = "mv " + ss_restraints_file_name + " " + custom_geom_file_name + " " + output_dir_w_CC
      libtbx.easy_run.fully_buffered(mv_command_string)
    '''
    
    '''
    ###############  (begin) when optimizing map_weight many times
    initial_start_temp = self.params.start_temperature
    initial_final_temp = self.params.final_temperature
    
    temp_interval = (initial_start_temp - initial_final_temp)/2
    
    first_iteration = True
    current_start_temp = ''
    user_defined_map_weight = ''
    while (current_start_temp != initial_final_temp):
      
      if (first_iteration == True):
        current_start_temp = initial_start_temp
        first_iteration = False
      else:
        current_start_temp = current_start_temp - temp_interval
      current_final_temp = current_start_temp - temp_interval
      
      if (self.params.map_weight == None): # a user didn't specify map_weight
        self.params.map_weight = determine_optimal_weight_by_template(self, map_inp)
      else:
        user_defined_map_weight = True
      print ("used self.params.map_weight for ",current_start_temp, " and ", current_final_temp, ":", round(self.params.map_weight,1))
    
      
      task_obj = cryo_fit2_run.cryo_fit2_class(
        model             = model_inp,
        model_name        = self.data_manager.get_default_model_name(),
        map_inp           = map_inp,
        params            = self.params,
        out               = self.logger,
        map_name          = self.data_manager.get_default_real_map_name(),
        logfile           = logfile,
        output_dir        = output_dir)
      
      task_obj.validate()
      output_dir_w_CC = task_obj.run()
      if (user_defined_map_weight == False):
        self.params.map_weight == None
    ###############  (end) when optimizing map_weight many times
    '''
    
    header = "# Geometry restraints used for cryo_fit2\n"
    header += "# %s\n" % date_and_time()
      
    r = model_inp.restraints_as_geo(
        header=header,
        # Stuff for outputting ncs_groups
        #excessive_distance_limit=self.params.ncs.excessive_distance_limit)
        excessive_distance_limit=10)
    print ("r:",r)
    
    geometry_restraints_file_name = "used_geometry_restraints.txt"
    geo_file = open(geometry_restraints_file_name, "w")
    geo_file.write(r)
    geo_file.close()
    
    mv_command_string = "mv " + geometry_restraints_file_name + " " + output_dir_w_CC
    libtbx.easy_run.fully_buffered(mv_command_string)
    
    mv_command_string = "mv " + log_file_name + " " + output_dir_w_CC
    libtbx.easy_run.fully_buffered(mv_command_string)
    
    time_total_end = time.time()
    time_took = show_time(time_total_start, time_total_end)
    
    logfile.write(str("cryo_fit2"))
    logfile.write(str(time_took))
    logfile.write("\n")
    logfile.close()
