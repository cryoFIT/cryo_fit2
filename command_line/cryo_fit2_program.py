######to do list
# To run cryo_fit2 at windows as well, replace all unix command like libtbx.easy_run.fully_buffered

from __future__ import division, print_function
import iotbx.pdb
import iotbx.phil

from libtbx.phil import change_default_phil_values
import libtbx.phil
import libtbx.phil.command_line
from libtbx import phil
from libtbx.utils import multi_out
from libtbx.utils import Sorry
from libtbx import easy_pickle

import mmtbx
import cryo_fit2_run
import os, shutil, sys, time

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
start_temperature = 500
  .type = int
  .short_caption = Starting temperature of annealing in Kelvin
final_temperature = 300
  .type = int
  .short_caption = Final temperature of annealing in Kelvin
cool_rate = 50
  .type = int
  .short_caption = cooling rate of annealing in Kelvin
number_of_steps = 1000
  .type = int
  .short_caption = number_of_steps in phenix.dynamics
wx = 1
  .type = int
    .short_caption = weight for cryo-EM map
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


class Program(ProgramTemplate):

  description = '''
Program for running cryo_fit2.\n

Minimum required inputs:
  Model file
  Map file

How to run:
  phenix.cryo_fit2 model.pdb map.ccp4
  
Options:
  start_temperature            (default: 500)
  final_temperature            (default: 300)
  cool_rate                    (default: 50)
  number_of_steps              (default: 1,000)
  wx                           (cryo-EM map weight, default: 1)
  keep_origin                  (default: True)
  secondary_structure.enabled  (default: True) Most MD simulations tend to break secondary structure. Therefore, turning on this option is recommended. If HELIX/SHEET records are present in supplied .pdb file, automatic search of the existing secondary structures in the given input pdb file will not be executed.
  secondary_structure.protein.remove_outliers (default: True) False may be useful for very poor low-resolution structures
  output_dir                   (output folder prefix, default: output)
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
  def validate(self):
    print('Validating inputs', file=self.logger)
    self.data_manager.has_models(raise_sorry=True)
    if not (self.data_manager.has_real_maps()):
      raise Sorry("Supply a map file.")

  # ---------------------------------------------------------------------------
  def run(self): # this run function actually executed (8/6/2018)
    
    print ("wx", str(self.params.wx))
    print ("start_temperature", str(self.params.start_temperature))
    print ("final_temperature", str(self.params.final_temperature))
    print ("cool_rate", str(self.params.cool_rate))
    print ("number_of_steps", str(self.params.number_of_steps)) 
    
    time_total_start = time.time()
    print('User input model: %s' % self.data_manager.get_default_model_name(), file=self.logger)
    model = self.data_manager.get_model()
    
    print('User input map: %s' % self.data_manager.get_default_real_map_name(), file=self.logger)
    map_inp = self.data_manager.get_real_map()
    
    print ("self.params.pdb_interpretation.secondary_structure.enabled:",self.params.pdb_interpretation.secondary_structure.enabled)
    print ("self.params.pdb_interpretation.secondary_structure.protein.remove_outliers:",self.params.pdb_interpretation.secondary_structure.protein.remove_outliers)
    
    log = multi_out()
    out=sys.stdout
    log.register("stdout", out)
    
    splited = self.data_manager.get_default_model_name().split("/")
    model_name_wo_path = splited [len(splited)-1]
    
    if (model_name_wo_path == "tst_cryo_fit2_model.pdb"):
        self.params.start_temperature = 500
        self.params.final_temperature = 300
        self.params.cool_rate = 100
        self.params.number_of_steps = 1
    
    ss_restraints = self.params.pdb_interpretation.secondary_structure.enabled
    remove_outlier_ss_restraints = self.params.pdb_interpretation.secondary_structure.protein.remove_outliers
    
    # name output_dir
    output_dir_prefix = self.params.output_dir
    output_dir = str(output_dir_prefix) + \
                 "_start_" + str(self.params.start_temperature) + \
                 "_final_" + str(self.params.final_temperature) + \
                 "_cool_" + str(self.params.cool_rate) + \
                 "_step_" + str(self.params.number_of_steps) + \
                 "_wx_" + str(self.params.wx) + \
                 "_ss_" + str(ss_restraints) + \
                 "_remove_outlier_ss_restraints_" + str(remove_outlier_ss_restraints)
    
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
        
    os.makedirs(output_dir)
    
    log_file_name = os.path.join(output_dir, "cryo_fit2.overall_log.txt")
    
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
    task_obj.run()
    
    time_total_end = time.time()
    time_took = show_time(time_total_start, time_total_end)

    logfile.write(str("cryo_fit2"))
    logfile.write(str(time_took))
    logfile.write("\n")
    logfile.close()
    
    mv_command_string = "mv " + log_file_name + " " + output_dir
    libtbx.easy_run.fully_buffered(mv_command_string)
    