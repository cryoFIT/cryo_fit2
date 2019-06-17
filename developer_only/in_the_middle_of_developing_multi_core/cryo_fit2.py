from __future__ import division, print_function
import iotbx.pdb
import libtbx.phil
import libtbx.phil.command_line
from libtbx import easy_mp
from libtbx import phil
from libtbx.utils import multi_out
from libtbx.utils import Sorry
from libtbx import easy_pickle
import mmtbx.cryo_fit2.cryo_fit2
import os, sys, time

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
# =============================================================================


#'''
def cryo_fit2_by_multi_core(self, wx, start_temp, final_temp, cool_rate, number_of_steps):
  print ("wx", str(wx))
  print ("start_temperature", str(start_temp))
  print ("final_temperature", str(final_temp))
  print ("cool_rate", str(cool_rate))
  print ("number_of_steps", str(number_of_steps))
  
  time_total_start = time.time()
  print('User input model: %s' % self.data_manager.get_default_model_name(), file=self.logger)
  model = self.data_manager.get_model()
  
  print('User input map: %s' % self.data_manager.get_default_real_map_name(), file=self.logger)
  map_inp = self.data_manager.get_real_map()
  
  log = multi_out()
  out=sys.stdout
  log.register("stdout", out)
  
  log_file_name = "cryo_fit2.overall_log.txt"
  logfile = open(log_file_name, "w") # since it is 'w', an existing file with the same name will be erased
  log.register("logfile", logfile)
  
  # redefine params for below
  params.start_temperature = start_temp
  params.final_temperature = final_temp
  params.cool_rate = cool_rate
  params.number_of_steps = number_of_steps
  params.wx = wx
  
  task_obj = mmtbx.cryo_fit2.cryo_fit2.cryo_fit2(
    model             = model,
    model_name        = self.data_manager.get_default_model_name(),
    map_inp           = map_inp,
    params            = self.params,
    out               = self.logger,
    map_name          = self.data_manager.get_default_real_map_name(),
    logfile           = logfile)
  
  task_obj.validate()
  task_obj.run()
  
  time_total_end = time.time()
  time_took = show_time(time_total_start, time_total_end)

  logfile.write(str("cryo_fit2"))
  logfile.write(str(time_took))
  logfile.write("\n")
  
  logfile.write(str("(To see animation, open output/all_states.pdb in pymol and click right bottom play button)"))
  logfile.write("\n")
  
  logfile.close()
############ end of cryo_fit2_by_multi_core()


program_citations = libtbx.phil.parse('''
citation {
  article_id = Pavel_s_2018_realspace_refine_mentioned_this phenix.dynamics
}
''')
# =============================================================================

master_phil_str = '''
include scope libtbx.phil.interface.tracking_params
start_temperature = 1000
  .type = int
  .short_caption = Starting temperature of annealing in Kelvin
final_temperature = 0
  .type = int
  .short_caption = Final temperature of annealing in Kelvin
cool_rate = 500
  .type = int
  .short_caption = cooling rate of annealing in Kelvin
number_of_steps = 100
  .type = int
  .short_caption = number_of_steps in phenix.dynamics
wx = 100
  .type = int
    .short_caption = weight for cryo-EM map
output_dir = output
  .type = path
  .short_caption = Output directory
use_multi_core = False
    .type=bool
    .help = If True, use multiple cores
            If false, use single core
    .short_caption = Keep origin of resulted model
keep_origin = True
    .type=bool
    .help = If True, write out model with origin in original location.  \
            If false, shifted origin to (0,0,0). 
    .short_caption = Keep origin of resulted model
''' # end of master_phil_str

# =============================================================================

class Program(ProgramTemplate):

  description = '''
Program for running cryo_fit2.\n

Minimum required inputs:
  Model file
  Map file

How to run:
  phenix.cryo_fit2 model.pdb map.ccp4
'''

  datatypes = ['model', 'real_map', 'phil']
  citations = program_citations
  master_phil_str = master_phil_str # this is ESSENTIAL to avoid
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

    wx = self.params.wx
    print ("wx", str(wx))
    
    if (use_multi_core == True):
      start_temp = self.params.start_temperature
      print ("start_temperature", str(start_temp))
      
      final_temp = self.params.final_temperature
      final_temp = 0
      print ("final_temperature", str(final_temp))
      
      cool_rate = self.params.cool_rate
      print ("cool_rate", str(cool_rate))
      
      number_of_steps = self.params.number_of_steps
      print ("number_of_steps", str(number_of_steps)) 
        
      argstuples = [( self, wx, start_temp, final_temp, cool_rate, number_of_steps), \
                    ( self, wx*10, start_temp, final_temp, cool_rate, number_of_steps), \
                    ( self, wx*10, start_temp, final_temp, cool_rate, number_of_steps*10), \
                    ( self, wx*100, start_temp, final_temp, cool_rate, number_of_steps*100), \
                    ( self, wx*1000, start_temp, final_temp, cool_rate, number_of_steps*1000) ]
      for args, res, errstr in easy_mp.multi_core_run( cryo_fit2_by_multi_core, argstuples, 4): # the last argument is nproc
        #print ('arguments: %s \nresult: %s \nerror: %s\n' %(args, res, errstr))
        #print ('arguments: %s \n' %(args))
        print ('arguments: ', str(args))
        #print ("wx", str(wx))
    else: # use single core
      print ("start_temperature", str(self.params.start_temperature))
      print ("final_temperature", str(self.params.final_temperature))
      print ("cool_rate", str(self.params.cool_rate))
      print ("number_of_steps", str(self.params.number_of_steps)) 
    
      time_total_start = time.time()
      print('User input model (before cleaning): %s' % self.data_manager.get_default_model_name(), file=self.logger)
      model = self.data_manager.get_model()
      
      print('User input map: %s' % self.data_manager.get_default_real_map_name(), file=self.logger)
      map_inp = self.data_manager.get_real_map()
      
      log = multi_out()
      out=sys.stdout
      log.register("stdout", out)
      
      log_file_name = "cryo_fit2.overall_log.txt"
      logfile = open(log_file_name, "w") # since it is 'w', an existing file with the same name will be erased
      log.register("logfile", logfile)
      
      task_obj = mmtbx.cryo_fit2.cryo_fit2.cryo_fit2(
        model             = model,
        model_name        = self.data_manager.get_default_model_name(),
        map_inp           = map_inp,
        params            = self.params,
        out               = self.logger,
        map_name          = self.data_manager.get_default_real_map_name(),
        logfile           = logfile)
      
      task_obj.validate()
      task_obj.run()
      
      time_total_end = time.time()
      time_took = show_time(time_total_start, time_total_end)
  
      logfile.write(str("cryo_fit2"))
      logfile.write(str(time_took))
      logfile.write("\n")
      
      logfile.close()
      