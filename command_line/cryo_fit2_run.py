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
  
    params = sa.master_params().extract()    # because of params = sa.master_params().extract() above, core parameters need to be redefined
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
    write_this = "\ncc before cryo_fit2: " + str(cc_before_cryo_fit2) + "\n\n"
    print('%s' %(write_this))
    self.logfile.write(str(write_this))

    result = ''
    total_steps_so_far = 0

    if (self.params.record_states == False): # default choice to avoid > 160 GB memory issue with recording all states for L1 stalk
      states = None
    
    cycle_so_far = 0
    cc_1st_array = []
    cc_2nd_array = []
    
    splited_model_name = self.model_name[:-4].split("/")
    model_file_name_only = splited_model_name[len(splited_model_name)-1]
    
    #number_of_atoms_in_input_pdb = know_number_of_atoms_in_input_pdb(self.logfile, self.model_name)
    # tRNA      : 1,563
    # L1 stalk  : 3,289
    # Mg channel: 14,940
    
    # number_of_atoms_in_input_pdb seems irrelevant to cc_check_after_every_this_cycle assignment.
    # but Mg channel with 10k check took 10 days!
    
    cc_check_after_every_this_cycle = ''
    if (("tst_cryo_fit2" in model_file_name_only) == True):
      cc_check_after_every_this_cycle = 5
    else:
      cc_check_after_every_this_cycle = 500
  
    
  ########################### <begin> iterate until cryo_fit2 derived cc saturates
    best_cc_so_far = -999 # tRNA has a negative value of initial cc
    result = ''
    for i in range(100000000): # runs well with cryo_fit2.run_tests     #for i in range(1000000000): # fails with cryo_fit2.run_tests with too much memory (bigger than 30 GB)
      
      self.params.map_weight = self.params.map_weight * weight_multiply
      # This is the only place where weight_multiply is applied
      # 1x~10x of weight_multiply were not enough for L1 stalk fitting
      # up to 20x of weight_multiply, nucleic acid geometry was ok, 30x broke it
      
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
      except:
        write_this = "Failed during core phenix.dynamics run."
        print (write_this)
        self.logfile.write(str(write_this))
        return self.output_dir
        
      multiply_this = 1 + ((params.start_temperature-params.final_temperature)/params.cool_rate)
      total_steps_so_far = total_steps_so_far + int(params.number_of_steps*multiply_this)
      cc_after_small_MD = calculate_cc(map_data=map_data, model=self.model, resolution=self.params.resolution)

      write_this = "CC after this epoch (a small MD iteration): " + str(round(cc_after_small_MD, 7)) + "\n"
      self.logfile.write(str(write_this))

      if (total_steps != ''):
        if (total_steps_so_far >= total_steps):
          write_this = "\ntotal_steps_so_far (" + str(total_steps_so_far) + \
                       ") >= A specified total_steps (" + str(total_steps) + ")\n"
          print('%s' %(write_this))
          self.logfile.write(str(write_this))
          break
      elif (cycle_so_far < cc_check_after_every_this_cycle/2):
        cycle_so_far = cycle_so_far + 1
        cc_1st_array.append(cc_after_small_MD)
      else:
        cycle_so_far = cycle_so_far + 1
        cc_2nd_array.append(cc_after_small_MD)
      
      if (cycle_so_far >= cc_check_after_every_this_cycle):
        write_this = "cycle_so_far:" + str(cycle_so_far) + "\n"
        print('%s' %(write_this))
        self.logfile.write(str(write_this))

        if (cc_after_small_MD > best_cc_so_far):
          write_this = "current_cc (" + str(cc_after_small_MD) + ") > best_cc_so_far (" + str(best_cc_so_far) + "). Therefore, cryo_fit2 will run longer MD.\n\n"
          print('%s' %(write_this))
          self.logfile.write(str(write_this))

          best_cc_so_far = cc_after_small_MD
          cycle_so_far = 0 # reset
          cc_1st_array = [] # reset
          cc_2nd_array = [] # reset

          if (self.params.reoptimize_map_weight_after_each_epoch == True):
            ### old comment: (let me comment this out since reoptimizing map_weight takes time and may cause conflict during exploration)
            ### new comment: commenting this out didn't help incomplete problem
            #self.params.map_weight = reoptimize_map_weight_if_not_specified(self, user_map_weight, map_inp, str(self.output_dir))
            self.params.map_weight = reoptimize_map_weight_if_not_specified(self, user_map_weight, map_inp)
            # although preliminary (just 1 benchmark), reoptimizing map_weight after each epoch prolongs running time ~5x
            # I confirmed that reoptimizing map_weight_after_each_epoch did change result (cc, SS stat) significantly
          continue 

        write_this = "current_cc (" + str(cc_after_small_MD) + ") <= best_cc_so_far (" + str(best_cc_so_far) + ")\n"
        print('%s' %(write_this))
        self.logfile.write(str(write_this))

        if (np.mean(cc_2nd_array) > np.mean(cc_1st_array)):
          write_this = "mean of cc_2nd_array (" + str(np.mean(cc_2nd_array)) + ") > mean of cc_1st_array (" + str(np.mean(cc_1st_array)) + ")\n"
          print('%s' %(write_this))
          self.logfile.write(str(write_this))
        
          cycle_so_far = 0 # reset
          cc_1st_array = [] # reset
          cc_2nd_array = [] # reset
          
          if (self.params.reoptimize_map_weight_after_each_epoch == True):
            ### old comment: (let me comment this out since reoptimizing map_weight takes time and may cause conflict during exploration)
            ### new comment: commenting this out didn't help incomplete problem
            self.params.map_weight = reoptimize_map_weight_if_not_specified(self, user_map_weight, map_inp)
            # although preliminary (just 1 benchmark), reoptimizing map_weight after each epoch prolongs running time ~5x
            # I confirmed that reoptimizing map_weight_after_each_epoch did change result (cc, SS stat) significantly
        else:
          write_this = "mean of cc_2nd_array (" + str(np.mean(cc_2nd_array)) + ") <= mean of cc_1st_array (" + str(np.mean(cc_1st_array)) + ")\n"
          print('%s' %(write_this))
          self.logfile.write(str(write_this))
          
          if (total_steps_so_far < self.params.total_steps_for_exploration):
            write_this = "\ntotal_steps_so_far (" + str(total_steps_so_far) + \
                       ") < total_steps_for_exploration (" + str(self.params.total_steps_for_exploration) + ")\n"
            print('%s' %(write_this))
            self.logfile.write(str(write_this))
            continue
          else:
            write_this = "cc values are saturated\ntotal_steps_so_far: " + str(total_steps_so_far) + "\n"
            print('%s' %(write_this))
            self.logfile.write(str(write_this))
            break
######################### <end> iterate until cryo_fit2 derived cc saturates
    
    
    cc_after_cryo_fit2 = calculate_cc(map_data=map_data, model=self.model, resolution=self.params.resolution)
    write_this = "\ncc after cryo_fit2: " + str(round(cc_after_cryo_fit2, 4)) + "\n\n"
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
    
    returned = know_how_much_map_origin_moved(str(self.map_name))
    if (returned != "origin_is_all_zero" and self.params.keep_origin == True):
        write_this = "Restoring original xyz position for a cryo_fit2 fitted atomistic model\n\n"
        print (write_this)
        self.logfile.write(str(write_this))
        return_to_origin_of_pdb_file(fitted_file_name_w_path, returned[0], returned[1], returned[2], returned[3])

    try:
      bp_num_in_fitted_file, H_num_in_fitted_file, E_num_in_fitted_file = \
        count_bp_H_E_in_fitted_file(fitted_file_name_w_path, output_dir_w_CC, self.logfile)
    except:
        write_this = "An exception occurred. Maybe cryo_fit2 failed to run (\"nan\") for this condition:" + \
                     " cool_rate (" + str(round(params.cool_rate, 1))   + ")" + \
                     " MD_in_each_epoch (" + str(MD_in_each_epoch)      + ")" + \
                     " number_of_steps (" + str(number_of_steps)        + ")" + \
                     " start_temperature (" + str(start_temperature)    + ")" + \
                     " weight_multiply (" + str(weight_multiply)        + ")" + \
                     " final_temperature (" + str(final_temperature)    + ")" + \
                     " map_weight (" + str(round(params.map_weight,2))  + ")" + \
                     " total_steps (" + str(params.total_steps)  + ")" 
        print (write_this)
        self.logfile.write(str(write_this))
        return output_dir_w_CC
    
    # To save a regression test time
    #if (("tst_cryo_fit2" in self.data_manager.get_default_model_name()) == False): #"AttributeError: 'cryo_fit2_class' object has no attribute 'data_manager'"
    if (("tst_cryo_fit2" in fitted_file_name_w_path) == False): 
      calculate_RMSD(self, fitted_file_name_w_path)
      
    output_dir_final = output_dir_w_CC + "_bp_" + str(bp_num_in_fitted_file) \
                      + "_H_" + str(H_num_in_fitted_file) + "_E_" + str(E_num_in_fitted_file)
    if os.path.exists(output_dir_final):
      shutil.rmtree(output_dir_final)

    mv_command_string = "mv " + output_dir_w_CC + " " + output_dir_final
    libtbx.easy_run.fully_buffered(mv_command_string)

    return output_dir_final
############# end of run function
