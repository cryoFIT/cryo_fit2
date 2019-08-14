# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_SIGNALS_DEFAULT=1
#### (source) http://cci.lbl.gov/cctbx_sources/crys3d/command_line/hklview.py

from __future__ import division, print_function
from cctbx.geometry_restraints.base_geometry import Base_geometry
import iotbx.phil, libtbx
from libtbx import group_args
from libtbx.utils import Sorry
from iotbx import map_and_model
import mmtbx.utils, os, sys
from mmtbx.command_line import geometry_minimization # maybe for cctbx_project/cctbx/geometry_restraints/base_geometry.py
from mmtbx.dynamics import simulated_annealing as sa
from mmtbx.superpose import *
import numpy as np
import shutil
import scitbx.math, scitbx.math.superpose, subprocess

os.environ['BOOST_ADAPTBX_FPE_DEFAULT'] = "1"
os.environ['BOOST_ADAPTBX_SIGNALS_DEFAULT'] = "1"

######## <begin> needed to import util.py
path = subprocess.check_output(["which", "phenix.cryo_fit2"])
splited_path = path.split("/")
command_path = ''
for i in range(len(splited_path)-3):
  command_path = command_path + splited_path[i] + "/"
command_path = command_path + "modules/cryo_fit2/"
util_path = command_path + "util/"
sys.path.insert(0, util_path)
from util import *
######## <end> needed to import util.py


class cryo_fit2_class(object):
  def __init__(self, model, model_name, map_inp, params, out, map_name, logfile, output_dir, user_map_weight, weight_multiply):
    self.model             = model
    self.model_name        = model_name
    self.map_inp           = map_inp
    self.params            = params
    self.out               = out
    self.map_name          = map_name
    self.logfile           = logfile
    self.output_dir        = output_dir
    self.desc              = os.path.basename(model_name)
    self.user_map_weight   = user_map_weight
    self.weight_multiply   = weight_multiply

  def __execute(self):
    self.caller(self.write_geo_file,       "Write GEO file")

  def validate(self): # this functions runs
    assert not None in [self.model, self.params, self.out]
    
    # Sanity check for crystal symmetry
    if (self.map_inp is not None):
      self.cs_consensus = mmtbx.utils.check_and_set_crystal_symmetry(
        models = [self.model], map_inps=[self.map_inp])

  def run(self):
    
    hierarchy = self.model.get_hierarchy()
    map_data, grid_unit_cell = None, None
    # sanity check for map and model
    if self.map_inp is not None:
      base = map_and_model.input(
        map_data         = self.map_inp.map_data(),
        model            = self.model,
        crystal_symmetry = self.cs_consensus,
        box              = False)

      hierarchy = base.model().get_hierarchy()
      map_data = base.map_data()
      grid_unit_cell = self.map_inp.grid_unit_cell()
    hierarchy.atoms().reset_i_seq()
  
    # Initialize states accumulator # Pavel's original
    states = mmtbx.utils.states(
     pdb_hierarchy  = self.model.get_hierarchy(),
     xray_structure = self.model.get_xray_structure())

    states.add(sites_cart = self.model.get_xray_structure().sites_cart())
  
    params                         = sa.master_params().extract()    # because of params = sa.master_params().extract() above, core parameters need to be redefined
    params.start_temperature       = self.params.start_temperature
    params.final_temperature       = self.params.final_temperature
    params.cool_rate               = self.params.cool_rate
    params.number_of_steps         = self.params.number_of_steps
    
    total_steps = ''
    if (self.params.total_steps != None):
      total_steps   = self.params.total_steps

    params.update_grads_shift      = 0.
    params.interleave_minimization = False #Pavel will fix the error that occur when params.interleave_minimization=True
    #print ("params:",params) # object like <libtbx.phil.scope_extract object at 0x1146ae210>
    map_inp                        = self.map_inp
    user_map_weight                = self.user_map_weight
    weight_multiply                = self.weight_multiply
    
    cc_before_cryo_fit2 = round(calculate_cc(map_data=map_data, model=self.model, resolution=self.params.resolution), 4)
    write_this = "\nCC before cryo_fit2 (both exploration and final MD): " + str(cc_before_cryo_fit2) + "\n\n"
    print('%s' %(write_this))
    self.logfile.write(str(write_this))

    result = ''
    total_steps_so_far = 0

    if (self.params.record_states == False): # default choice to avoid > 160 GB memory issue with recording all states for L1 stalk
      states = None
    
    if (self.params.reoptimize_map_weight_after_each_cycle_during_final_MD == True):
      cycle_so_far_for_map_weight_reoptimization = 0
    
    splited_model_name = self.model_name[:-4].split("/")
    model_file_name_only = splited_model_name[len(splited_model_name)-1]
    
    #number_of_atoms_in_input_pdb = know_number_of_atoms_in_input_pdb(self.logfile, self.model_name)
    # tRNA      : 1,563
    # L1 stalk  : 3,289
    # Mg channel: 14,940
    
    # number_of_atoms_in_input_pdb seems irrelevant to check_cc_after_these_cycles assignment.
    # but Mg channel with 10k check took 10 days!
    
    
    ########################### <begin> prepare/initialize for iteration
    check_cc_after_these_steps = ''
    if (("tst_cryo_fit2" in model_file_name_only) == True):
      check_cc_after_these_steps = 500 # if too small like 100, it may run forever
    else:
      check_cc_after_these_steps = 100000
  
    reoptimize_map_weight_after_these_cycles = ''
    if (self.params.reoptimize_map_weight_after_each_cycle_during_final_MD == True):
      if (("tst_cryo_fit2" in model_file_name_only) == True):
        reoptimize_map_weight_after_these_cycles = 5
      else:
        reoptimize_map_weight_after_these_cycles = 100 # after 123~171 cycles, full tRNA crashes (when map_weight is multiplied too crazy back then,,,)
        
    if (("tst_cryo_fit2_" in self.model_name) == True): 
      self.params.total_steps_for_exploration = 100
    
    map_weight_before_multiplication = self.params.map_weight
    self.params.map_weight = self.params.map_weight * weight_multiply
    #### This is the only place where weight_multiply is applied (other than reoptimize_map_weight_if_not_specified for final MD)
    
    best_cc_so_far = -999 # tRNA has a negative value of initial cc
    cc_1st_array = []
    cc_2nd_array = []
    result = ''
    total_steps_so_far_for_cc_check = 0 # initialization
    ########################### <end> prepare/initialize for iteration
    
    
    ########################### <begin> iterate until cryo_fit2 derived cc saturates
    for i in range(100000000): # runs well with cryo_fit2.run_tests     #for i in range(1000000000): # fails with cryo_fit2.run_tests with too much memory (bigger than 30 GB)
      
      write_this = "\n" + str(i) + "th iteration with " + str(round(self.params.map_weight,1)) + " self.params.map_weight (after multiplication)\n"
      print (write_this)
      self.logfile.write(str(write_this))
      
      print ("(new iteration) check_cc_after_these_steps:",check_cc_after_these_steps)
      print ("total_steps_so_far_for_cc_check:",total_steps_so_far_for_cc_check)
      
      try:
        if (self.params.progress_on_screen == True): # default choice
          result = sa.run(
            params = params,
            xray_structure     = self.model.get_xray_structure(),
            restraints_manager = self.model.get_restraints_manager(),
            target_map         = map_data,
            real_space         = True,
            wx                 = self.params.map_weight, 
            wc                 = 1, # weight for geometry conformation
            states_collector   = states) 
        else: # (self.params.progress_on_screen = False):
          result = sa.run(
            params = params,
            xray_structure     = self.model.get_xray_structure(),
            restraints_manager = self.model.get_restraints_manager(),
            target_map         = map_data,
            real_space         = True,
            wx                 = self.params.map_weight, 
            wc                 = 1, # weight for geometry conformation
            states_collector   = states,
            log                = self.logfile) # if this is commented, temp= xx dist_moved= xx angles= xx bonds= xx is shown on screen rather than cryo_fit2.log
      except Exception as ex:
        write_this = "exception message:" +  str(ex)
        print (write_this)
        self.logfile.write(str(write_this))
        
        write_this = "Failed during core map weighted phenix.dynamics run."
        print (write_this)
        self.logfile.write(str(write_this))
        return self.output_dir
        
      multiply_this = 1 + ((params.start_temperature-params.final_temperature)/params.cool_rate)
      total_steps_so_far = total_steps_so_far + int(params.number_of_steps*multiply_this)
      
      
      cc_after_small_MD = calculate_cc(map_data=map_data, model=self.model, resolution=self.params.resolution)
      write_this = "CC after this cycle (a small MD iteration): " + str(round(cc_after_small_MD, 7)) + "\n"
      self.logfile.write(str(write_this))
      
      if (self.params.explore == True):
        if (total_steps_so_far < self.params.total_steps_for_exploration):
          write_this = "\ntotal_steps_so_far (" + str(total_steps_so_far) + \
                       ") < total_steps_for_exploration (" + str(self.params.total_steps_for_exploration) + ")\n"
          print('%s' %(write_this))
          self.logfile.write(str(write_this))
          continue
        else:
          write_this = "\ntotal_steps_so_far (" + str(total_steps_so_far) + \
                       ") >= total_steps_for_exploration (" + str(self.params.total_steps_for_exploration) + ")\n"
          print('%s' %(write_this))
          self.logfile.write(str(write_this))
          break
      
      
      ############# all below is for final MD
      total_steps_so_far_for_cc_check = total_steps_so_far_for_cc_check + int(params.number_of_steps*multiply_this)

      
      if (total_steps != ''):
    
        if (total_steps_so_far >= total_steps):
          write_this = "\ntotal_steps_so_far (" + str(total_steps_so_far) + \
                     ") >= A specified total_steps (" + str(total_steps) + ")\n"
          print('%s' %(write_this))
          self.logfile.write(str(write_this))
          break
        if (self.params.reoptimize_map_weight_after_each_cycle_during_final_MD == True):
          cycle_so_far_for_map_weight_reoptimization = cycle_so_far_for_map_weight_reoptimization + 1
      elif (total_steps_so_far_for_cc_check < check_cc_after_these_steps/2):
        if (self.params.reoptimize_map_weight_after_each_cycle_during_final_MD == True):
          cycle_so_far_for_map_weight_reoptimization = cycle_so_far_for_map_weight_reoptimization + 1
        cc_1st_array.append(cc_after_small_MD)
      else:
        if (self.params.reoptimize_map_weight_after_each_cycle_during_final_MD == True):
          cycle_so_far_for_map_weight_reoptimization = cycle_so_far_for_map_weight_reoptimization + 1
        cc_2nd_array.append(cc_after_small_MD)
      
      if (self.params.reoptimize_map_weight_after_each_cycle_during_final_MD == True):
        if (cycle_so_far_for_map_weight_reoptimization >= reoptimize_map_weight_after_these_cycles):
          self.params.map_weight = reoptimize_map_weight_if_not_specified(self, user_map_weight, map_inp)
          self.params.map_weight = self.params.map_weight * weight_multiply
          cycle_so_far_for_map_weight_reoptimization = 0 # reinitialization
          # I confirmed that reoptimizing map_weight_after_each_cycle did change result (cc, SS stat) significantly
      
      if (total_steps_so_far_for_cc_check >= check_cc_after_these_steps):
        
        if (cc_after_small_MD > best_cc_so_far):
          write_this = "current_cc (" + str(cc_after_small_MD) + ") > best_cc_so_far (" + str(best_cc_so_far) + "). Therefore, cryo_fit2 will run longer MD.\n\n"
          print('%s' %(write_this))
          self.logfile.write(str(write_this))

          best_cc_so_far = cc_after_small_MD
          cc_1st_array = [] # reset
          cc_2nd_array = [] # reset
          total_steps_so_far_for_cc_check = 0 # reset
          
          continue 

        else:
          write_this = "current_cc (" + str(cc_after_small_MD) + ") <= best_cc_so_far (" + str(best_cc_so_far) + ")\n"
          print('%s' %(write_this))
          self.logfile.write(str(write_this))

        if ((len(cc_1st_array) == 0) or (len(cc_2nd_array) == 0)):
          print ("(len(cc_1st_array) == 0) or (len(cc_2nd_array) == 0)")
          continue
          
        if (np.mean(cc_2nd_array) > np.mean(cc_1st_array)):
          write_this = "mean of cc_2nd_array (" + str(np.mean(cc_2nd_array)) + ") > mean of cc_1st_array (" + str(np.mean(cc_1st_array)) + ")\n"
          print('%s' %(write_this))
          self.logfile.write(str(write_this))

          cc_1st_array = [] # reset
          cc_2nd_array = [] # reset
          total_steps_so_far_for_cc_check = 0 # reset

        else: #(np.mean(cc_2nd_array) <= np.mean(cc_1st_array)):
          write_this = "mean of cc_2nd_array (" + str(np.mean(cc_2nd_array)) + ") <= mean of cc_1st_array (" + str(np.mean(cc_1st_array)) + ")\n"
          print('%s' %(write_this))
          self.logfile.write(str(write_this))
          
          write_this = "cc values are saturated\ntotal_steps_so_far: " + str(total_steps_so_far) + "\n"
          print('%s' %(write_this))
          self.logfile.write(str(write_this))
          break
######################### <end> iterate until cryo_fit2 derived cc saturates
    
    
    cc_after_cryo_fit2 = calculate_cc(map_data=map_data, model=self.model, resolution=self.params.resolution)
    
    if (self.params.explore == False): # no need to calculate cc after explore
      write_this = "\nCC after final MD of cryo_fit2: " + str(round(cc_after_cryo_fit2, 4)) + "\n\n"
      print('%s' %(write_this))
      self.logfile.write(str(write_this))
    
    output_dir_w_CC = str(self.output_dir) + "_cc_" + str(round(cc_after_cryo_fit2, 3))
    if os.path.exists(output_dir_w_CC):
      shutil.rmtree(output_dir_w_CC)
    os.mkdir(output_dir_w_CC)
    
    if (self.params.record_states == True): 
      all_state_file = os.path.join(output_dir_w_CC, "all_states.pdb")
      states.write(file_name = all_state_file)

    self.model.set_xray_structure(result.xray_structure)

    fitted_file_name = model_file_name_only + "_cryo_fit2_fitted.pdb"
    fitted_file_name_w_path = os.path.join(output_dir_w_CC, fitted_file_name)
    
    ##### this is essential to spit cryo_fitted2 file
    with open(fitted_file_name_w_path, "w") as f:
      f.write(self.model.model_as_pdb())
    f.close()
    
    #print_this ='''
########  How to fix map origin problem in cryo_fit2 #######
    '''
  With 0,0,0 origin map, cryo_fit2 has no problem.
  However, with non-0,0,0 origin cryo-EM map, cryo_fit2 results cryo_fitted pdb model at "wrong" origin
  This is because probably dynamics part uses map at 0,0,0 origin.
  Therefore, cryo_fit2 identifies how much the map origin was moved, then update all xyz coordinates of output pdb file.
  In user's perspective, there is nothing to bother.
  All kinds of mrc files (e.g. "Regular", emdb downloaded, went through phenix.map_box, gaussian filtered by UCSF Chimera and went through relion_image_handler) work fine.
#############################################################
#print (print_this,"\n")
    '''    

    try:
      bp_num_in_fitted_file, H_num_in_fitted_file, E_num_in_fitted_file = \
        count_bp_H_E_in_fitted_file(fitted_file_name_w_path, output_dir_w_CC, self.logfile)
    except Exception as ex:
      write_this = "exception message:" +  str(ex)
      print (write_this)
      self.logfile.write(str(write_this))
      write_this = "(in task_obj loop) An exception occurred in cryo_fit2_run. Maybe cryo_fit2 failed to run (\"nan\") for this condition:" + \
                   " cool_rate (" + str(round(params.cool_rate, 1))          + ")\n" + \
                   " number_of_steps (" + str(params.number_of_steps)        + ")\n" + \
                   " start_temperature (" + str(params.start_temperature)    + ")\n" + \
                   " weight_multiply (" + str(weight_multiply)               + ")\n" + \
                   " final_temperature (" + str(params.final_temperature)    + ")\n" + \
                   " map_weight (" + str(round(self.params.map_weight,2))    + ")\n" + \
                   " total_steps (" + str(total_steps)  + ")"  # total_steps alone is ok without params, self.params
      print (write_this)
      self.logfile.write(str(write_this))
      
      if (os.path.isdir("parameters_exploration/bp_H_E_not_calculated") == False):
        os.mkdir("parameters_exploration/bp_H_E_not_calculated")
      command_string = "mv " + str(output_dir_w_CC) + " parameters_exploration/bp_H_E_not_calculated"
      logfile.write(str(command_string))
      libtbx.easy_run.fully_buffered(command=command_string).raise_if_errors().stdout_lines
            
      return output_dir_w_CC
    
    returned = know_how_much_map_origin_moved(str(self.map_name))
    if (returned != "origin_is_all_zero" and self.params.keep_origin == True):
        write_this = "Restoring original xyz position for a cryo_fit2 fitted atomistic model\n\n"
        print (write_this)
        self.logfile.write(str(write_this))
        return_to_origin_of_pdb_file(fitted_file_name_w_path, returned[0], returned[1], returned[2], returned[3])
        
    if (("tst_cryo_fit2" in fitted_file_name_w_path) == False): 
      calculate_RMSD(self, fitted_file_name_w_path)
      
    output_dir_final = output_dir_w_CC + "_bp_" + str(bp_num_in_fitted_file) \
                      + "_H_" + str(H_num_in_fitted_file) + "_E_" + str(E_num_in_fitted_file)
    if os.path.exists(output_dir_final):
      shutil.rmtree(output_dir_final)

    mv_command_string = "mv " + output_dir_w_CC + " " + output_dir_final
    libtbx.easy_run.fully_buffered(mv_command_string)
    
    
    ############################
    current_dir = os.getcwd()
    os.chdir(output_dir_final)

    command_string = "echo " + str(map_weight_before_multiplication) + " >> used_map_weight_before_multiplication.txt"
    libtbx.easy_run.fully_buffered(command=command_string).raise_if_errors().stdout_lines

    os.chdir(current_dir)    
    ############################


    return output_dir_final
############# end of run function
