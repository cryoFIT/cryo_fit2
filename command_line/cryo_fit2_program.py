######to do list
# To run cryo_fit2 at windows as well, replace all unix command like libtbx.easy_run.fully_buffered

from __future__ import division, print_function

from cryo_fit2_util import *

from libtbx import easy_mp
from libtbx import easy_pickle
from libtbx import phil
import iotbx.pdb

import iotbx.phil
from libtbx.phil import change_default_phil_values
import libtbx.phil
import libtbx.phil.command_line

from libtbx.utils import multi_out
from libtbx.utils import Sorry

import mmtbx
import cryo_fit2_run
import os, shutil, sys, time

from mmtbx.refinement.real_space import weight

try:
  from phenix.program_template import ProgramTemplate # at Doonam's laptop, phenix.program_template works
except ImportError:
  from libtbx.program_template import ProgramTemplate
  

def show_time(time_start, time_end):
  time_took = 0 # temporary of course
  if (round((time_end-time_start)/60, 1) < 1):
    time_took = " finished in " + str(round((time_end-time_start), 2)) + " seconds (wallclock)."
  elif (round((time_end-time_start)/60/60, 1) < 1):
    time_took = " finished in " + str(round((time_end-time_start)/60, 2)) + " minutes (wallclock)."
  else:
    time_took = " finished in " + str(round((time_end-time_start)/60/60, 1)) + " hours (wallclock)."
  return time_took
############### end of show_time function


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
number_of_steps = 100
  .type = int
  .short_caption = number_of_steps in phenix.dynamics
map_weight = None
  .type = int
  .short_caption = cryo-EM map weight. A user is recommended NOT to specify this, so that it will be automatically determined.
resolution = None
  .type = int
  .short_caption = cryo-EM map resolution (Angstrom) that needs to be specified by a user
output_dir = output
  .type = path
  .short_caption = Output folder PREFIX
keep_origin = True
    .type = bool
    .help = If True, write out model with origin in original location.  \
            If False, shift origin to (0,0,0). 
    .short_caption = Keep origin of a resulted atomic model
include scope mmtbx.monomer_library.pdb_interpretation.grand_master_phil_str # to use secondary_structure.enabled
''' ############## end of base_master_phil_str


new_default = 'pdb_interpretation.secondary_structure.enabled = True'
modified_master_phil_str = change_default_phil_values(
  base_master_phil_str, new_default, phil_parse=iotbx.phil.parse)

new_default = 'pdb_interpretation.secondary_structure.protein.remove_outliers = True'
modified_master_phil_str = change_default_phil_values(
  base_master_phil_str, new_default, phil_parse=iotbx.phil.parse)

'''
new_default = 'pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_planarity = False'
modified_master_phil_str = change_default_phil_values(
  modified_master_phil_str, new_default, phil_parse=iotbx.phil.parse)
'''

# new_default3 = 'pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_hbonds = True'
# modified_master_phil_str = change_default_phil_values(
#   modified_master_phil_str, new_default3, phil_parse=iotbx.phil.parse)

class Program(ProgramTemplate):

  description = '''
Program for running cryo_fit2.\n

Minimum required inputs:
  Model file
  Map file
  Map resolution information

Example running command:
  phenix.cryo_fit2 model.pdb map.ccp4 resolution=6
  
Options:
  resolution                   (cryo-EM map resolution in Angstrom that needs to be entered by a user)
  map_weight                   (cryo-EM map weight. A user is recommended NOT to specify this, so that it will be automatically determined.)
  start_temperature            (default: 300)
  final_temperature            (default: 0)
  cool_rate                    (default: 10)
  number_of_steps              (default: 100)
  secondary_structure.enabled  (default: True)
                               Most MD simulations tend to break secondary structure. 
                               Therefore, turning on this option is recommended. 
                               If HELIX/SHEET records are present in supplied .pdb file, 
                               automatic search of the existing secondary structures in the given 
                               input pdb file will not be executed.
  secondary_structure.protein.remove_outliers (default: True)
                               False may be useful for very poor low-resolution structures by
                               ignoring some hydrogen "bond" if it exceed certain distance threshold
  secondary_structure.nucleic_acid.base_pair.restrain_planarity  (default: True)
  secondary_structure.nucleic_acid.base_pair.restrain_hbonds  (default: True)
  output_dir                   (output folder name prefix, default: output)
  keep_origin                  (default: True)
'''

  datatypes = ['model', 'real_map', 'phil']
  citations = program_citations
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
      raise Sorry("Map resolution is required. Type \"phenix.cryo_fit2\" to know minimally required options")

  # ---------------------------------------------------------------------------
  def run(self): # this run function actually executed (8/6/2018)
    
    print ("user entered resolution", str(self.params.resolution))
    print ("start_temperature", str(self.params.start_temperature))
    print ("final_temperature", str(self.params.final_temperature))
    print ("cool_rate", str(self.params.cool_rate))
    print ("number_of_steps", str(self.params.number_of_steps)) 
    
    time_total_start = time.time()
    print('User input model: %s' % self.data_manager.get_default_model_name(), file=self.logger)
    model = self.data_manager.get_model()
    
    from iotbx import pdb  #contains hierarchy data structure
    pdb_io = pdb.input(self.data_manager.get_default_model_name())
    hierarchy = pdb_io.construct_hierarchy()
    
    pdb_str_1 = hierarchy.as_pdb_string()
    pdb_str_1 = "CRYST1   44.034   76.843   61.259  90.00  90.00  90.00 P 1\n" + pdb_str_1
    # add this bogus cryst to avoid "Sorry: Crystal symmetry is missing or cannot be extracted." in get_pdb_inputs
    
    if (self.params.map_weight == None): # a user didn't specify map_weight
      self.params.map_weight = determine_optimal_weight(self.params.resolution, pdb_str_1)
      self.params.map_weight = self.params.map_weight/3
      # original determine_optimal_weight seems for RSR, since dynamics changes conformation more, dividing by 2 seems reasonable.
      # indeed, test w/ a helix works perfectly with self.params.map_weight/3
      print ("optimized weight", str(self.params.map_weight))
        
    print('User input map: %s' % self.data_manager.get_default_real_map_name(), file=self.logger)
    map_inp = self.data_manager.get_real_map()
    
    print ("self.params.pdb_interpretation.secondary_structure.enabled:",self.params.pdb_interpretation.secondary_structure.enabled)
    print ("self.params.pdb_interpretation.secondary_structure.protein.remove_outliers:",self.params.pdb_interpretation.secondary_structure.protein.remove_outliers)
    
    #print ("self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_planarity:",self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_planarity)
    #print ("self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_hbonds:",self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_hbonds)
    
    log = multi_out()
    out=sys.stdout
    log.register("stdout", out)
    
    splited = self.data_manager.get_default_model_name().split("/")
    model_name_wo_path = splited [len(splited)-1]
  
    if ((model_name_wo_path == "devel_cryo_fit2_model.pdb") or (model_name_wo_path == "tst_cryo_fit2_model.pdb")):
      # "tst..." lives in modules/cryo_fit2/regression
      self.params.start_temperature = 500
      self.params.final_temperature = 300
      self.params.cool_rate = 1
      self.params.number_of_steps = 10
      self.params.pdb_interpretation.secondary_structure.enabled = True
    
    if (model_name_wo_path == "tutorial_cryo_fit2_model.pdb"): 
      self.params.start_temperature = 1000
      self.params.final_temperature = 0
      self.params.cool_rate = 10
      self.params.number_of_steps = 1000
      self.params.pdb_interpretation.secondary_structure.enabled = True
      
    ss_restraints = self.params.pdb_interpretation.secondary_structure.enabled
    remove_outlier_ss_restraints = self.params.pdb_interpretation.secondary_structure.protein.remove_outliers
    
    # name output_dir
    output_dir_prefix = self.params.output_dir
    output_dir = str(output_dir_prefix) + \
                 "_map_resolution_" + str(self.params.resolution) + \
                 "_map_weight_" + str(round(self.params.map_weight,1)) + \
                 "_start_" + str(self.params.start_temperature) + \
                 "_final_" + str(self.params.final_temperature) + \
                 "_cool_" + str(self.params.cool_rate) + \
                 "_step_" + str(self.params.number_of_steps) + \
                 "_ss_" + str(ss_restraints) + \
                 "_remove_outlier_ss_restraints_" + str(remove_outlier_ss_restraints) 
    
    log_file_name = "cryo_fit2.log"
    
    logfile = open(log_file_name, "w") # since it is 'w', an existing file with the same name will be erased
    log.register("logfile", logfile)
    
    task_obj = cryo_fit2_run.cryo_fit2_class(
      model             = model,
      model_name        = self.data_manager.get_default_model_name(),
      map_inp           = map_inp,
      params            = self.params,
      out               = self.logger,
      map_name          = self.data_manager.get_default_real_map_name(),
      logfile           = logfile,
      output_dir        = output_dir)
    
    task_obj.validate()
    output_dir_w_CC = task_obj.run()
    
    time_total_end = time.time()
    time_took = show_time(time_total_start, time_total_end)

    logfile.write(str("cryo_fit2"))
    logfile.write(str(time_took))
    logfile.write("\n")
    logfile.close()
    
    mv_command_string = "mv " + log_file_name + " " + output_dir_w_CC
    libtbx.easy_run.fully_buffered(mv_command_string)
    