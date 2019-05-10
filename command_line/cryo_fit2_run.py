from __future__ import division, print_function
import iotbx.phil, libtbx
from libtbx import group_args
from libtbx.utils import Sorry
from iotbx import map_and_model
import mmtbx.utils, os, sys
from mmtbx.dynamics import simulated_annealing as sa
from mmtbx.superpose import *
import shutil
import scitbx.math, scitbx.math.superpose, subprocess

from mmtbx.command_line import geometry_minimization # maybe for cctbx_project/cctbx/geometry_restraints/base_geometry.py
from cctbx.geometry_restraints.base_geometry import Base_geometry


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
  def __init__(self, model, model_name, map_inp, params, out, map_name, logfile, output_dir, user_map_weight):
    self.model             = model
    self.model_name        = model_name
    self.map_inp           = map_inp
    self.params            = params
    self.out               = out
    self.map_name          = map_name
    self.logfile           = logfile
    self.output_dir        = output_dir
    self.desc = os.path.basename(model_name)
    self.user_map_weight   = user_map_weight
  
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
  
    #'''  
    ########### Pavel original
    # Initialize states accumulator
    states = mmtbx.utils.states(
     pdb_hierarchy  = self.model.get_hierarchy(),
     xray_structure = self.model.get_xray_structure())
    
    states.add(sites_cart = self.model.get_xray_structure().sites_cart())
    #'''
  
  
  
    params = sa.master_params().extract()
    # because of params = sa.master_params().extract() above, core parameters need to be redefined
    params.start_temperature       = self.params.start_temperature
    params.final_temperature       = self.params.final_temperature
    params.cool_rate               = self.params.cool_rate
    params.number_of_steps         = self.params.number_of_steps
    
    total_number_of_steps = ''
    if (self.params.total_number_of_steps != None):
      total_number_of_steps   = self.params.total_number_of_steps
    
    params.update_grads_shift      = 0.
    params.interleave_minimization = False #Pavel will fix the error that occur when params.interleave_minimization=True
    
    map_inp                        = self.map_inp
    user_map_weight                = self.user_map_weight
    
    cc_before_cryo_fit2 = round(calculate_cc(map_data=map_data, model=self.model, resolution=self.params.resolution), 4)
    write_this = "\ncc before cryo_fit2: " + str(cc_before_cryo_fit2) + "\n\n"
    
    print('%s' %(write_this))
    self.logfile.write(str(write_this))
    #print ("params:",params) # object like <libtbx.phil.scope_extract object at 0x1146ae210>
    
    
    ################ <begin> iterate until cryo_fit2 derived cc saturates
    result = ''
    total_number_of_steps_so_far = 0

    if (self.params.record_states == False): # default
      states = None
    
    for i in range(100000000): # runs well with cryo_fit2.run_tests
    #for i in range(1000000000): # fails with cryo_fit2.run_tests with too much memory (bigger than 30 GB)
 
      if (self.params.progress_on_screen == True):
        #'''
        result = sa.run(
          params = params,
          xray_structure     = self.model.get_xray_structure(),
          restraints_manager = self.model.get_restraints_manager(),
          target_map         = map_data,
          real_space         = True,
          wx                 = self.params.map_weight, 
          wc                 = 1, # weight for geometry conformation
          states_collector   = states) # we may need not using this to help 150 GB memory problem?\
                                       #, even when this is commented, all_states.pdb is still produced\
                                       #, even when this is None, all_states.pdb is still produced
        #'''
        
        '''
        result = sa.run(
          params = params,
          xray_structure     = self.model.get_xray_structure(),
          restraints_manager = self.model.get_restraints_manager(),
          target_map         = map_data,
          real_space         = True,
          wx                 = self.params.map_weight, 
          wc                 = 1) # weight for geometry conformation
        '''
        
        '''
        result = sa.run(
          params = params,
          xray_structure     = self.model.get_xray_structure(),
          restraints_manager = self.model.get_restraints_manager(),
          target_map         = map_data,
          real_space         = True,
          wx                 = self.params.map_weight, 
          wc                 = 1, # weight for geometry conformation
          states_collector   = None)
        '''
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
      
      multiply_this = 1 + ((params.start_temperature-params.final_temperature)/params.cool_rate)
      total_number_of_steps_so_far = total_number_of_steps_so_far + params.number_of_steps*multiply_this
      cc_after_cryo_fit2 = calculate_cc(map_data=map_data, model=self.model, resolution=self.params.resolution)

      write_this = "cc after  cryo_fit2: " + str(round(cc_after_cryo_fit2, 4)) + "\n"
      print('%s' %(write_this))
      self.logfile.write(str(write_this))
      
      if (total_number_of_steps != ''):
        if (total_number_of_steps_so_far >= total_number_of_steps):
          write_this = "\ntotal_number_of_steps_so_far " + str(total_number_of_steps_so_far) + \
                       " >= total_number_of_steps " + str(total_number_of_steps) + "\n"
          print('%s' %(write_this))
          self.logfile.write(str(write_this))
          break
      elif (cc_after_cryo_fit2 <= cc_before_cryo_fit2):
        write_this = "\ntotal_number_of_steps_so_far: " + str(total_number_of_steps_so_far) + "\n"
        print('%s' %(write_this))
        self.logfile.write(str(write_this))
        break
      cc_before_cryo_fit2 = cc_after_cryo_fit2 # reassign cc_before_cryo_fit2
    ################ <end> iterate until cryo_fit2 derived cc saturates
    
    '''
    ################ <begin> final cryo_fit2 run for better geometry with a new map_weight
    write_this = "\nFinal cryo_fit2 run for better geometry with a new map_weight\n"
    print('%s' %(write_this))
    self.logfile.write(str(write_this))
    
    if (user_map_weight == ''):
      write_this = "User didn't specify map_weight. Therefore, automatically optimize map_weight for final cryo_fit2 run\n"
      print('%s' %(write_this))
      self.logfile.write(str(write_this))
    
      final = True

      fitted_file_before_final_run = "fitted_file_before_final_run.pdb"
      with open(fitted_file_before_final_run, "w") as f:
        f.write(self.model.model_as_pdb())
      f.close()

      self.params.map_weight = determine_optimal_weight_by_template(self, self.logfile, map_inp, final, fitted_file_before_final_run)
      
      cmd = "rm fitted_file_before_final_run.pdb"
      libtbx.easy_run.fully_buffered(cmd)
      
      write_this = "\nAutomatically optimized "
      print('%s' %(write_this))
      self.logfile.write(write_this)
    else:
      self.params.map_weight = user_map_weight
    
    write_this = "map_weight for final cryo_fit2 run: " + str(round(self.params.map_weight,1)) + "\n"
    print('%s' %(write_this))
    self.logfile.write(str(write_this))
    
    if (self.params.progress_on_screen == True):
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
    
    write_this = "\ncc after cryo_fit2 (final): " + str(round(cc_after_cryo_fit2, 4)) + "\n\n"
    print('%s' %(write_this))
    self.logfile.write(str(write_this))
    
    ################ <end> final cryo_fit2 run for better geometry with a new map_weight
    '''
    
    output_dir_w_CC = str(self.output_dir) + "_cc_" + str(round(cc_after_cryo_fit2, 3))
    if os.path.exists(output_dir_w_CC):
      shutil.rmtree(output_dir_w_CC)
    os.mkdir(output_dir_w_CC)
    
    '''
    all_state_file = os.path.join(output_dir_w_CC, "all_states.pdb")
    states.write(file_name = all_state_file)
    '''
    
    self.model.set_xray_structure(result.xray_structure)
    
    splited_model_name = self.model_name[:-4].split("/")
    model_file_name_only = splited_model_name[len(splited_model_name)-1] 
    fitted_file_name = model_file_name_only + "_cryo_fit2_fitted.pdb"
    fitted_file = os.path.join(output_dir_w_CC, fitted_file_name)
    
    ##### this is essential
    with open(fitted_file, "w") as f:
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
        return_to_origin_of_pdb_file(fitted_file, returned[0], returned[1], returned[2], returned[3])
    
    
    ################################ <begin> RMSD calculation ###########################
    ## (reference) cctbx_project/mmtbx/superpose.py
    fixed = self.model_name
    moving = fitted_file
    
    #print ("\n===== Init =====")
    write_this = "\n===== Init ====="
    print (write_this)
    self.logfile.write(str(write_this))
    # The fixed model can only contain a single model.
    # It will raise an Exception if there is more than one!
    fixed = SuperposePDB(
      fixed,
      selection=self.params.selection_fixed,
      preset=self.params.selection_fixed_preset,
      log=None,
      quiet=False,
      desc=fixed
    )
    
    # The moving pdb can contain many models. These will each be aligned to the
    # fixed model and output as a separate file...
    moving_args = dict(
      selection=self.params.selection_moving,
      preset=self.params.selection_moving_preset,
      desc=moving,
      log=None,
      quiet=False
    )
    for count, moving in enumerate(SuperposePDB.open_models(moving, **moving_args)):
      #print ("\n===== Aligning %s to %s ====="%(fitted_file, self.model_name))
      write_this = "\n===== Aligning %s to %s ====="%(fitted_file, self.model_name)
      print (write_this)
      self.logfile.write(str(write_this))
      if not self.params.selection_moving:
        moving.selectomatic(fixed)
      rmsd, lsq = moving.superpose(fixed)
      write_this = "\nrmsd after cryo_fit2: " + str(round(rmsd,2)) + " angstrom\n\n"
      print (write_this)
      self.logfile.write(str(write_this))
    ################################ <end> RMSD calculation ###########################
    
    return output_dir_w_CC
############# end of run function
