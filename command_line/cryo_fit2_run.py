from __future__ import division, print_function
import iotbx.phil, libtbx
from libtbx import group_args
from libtbx.utils import Sorry
from iotbx import map_and_model
import mmtbx.utils, os
from mmtbx.dynamics import simulated_annealing as sa
import shutil
from cryo_fit2_util import *

class cryo_fit2_class(object):
  def __init__(self, model, model_name, map_inp, params, out, map_name, logfile, output_dir):
    self.model             = model
    self.model_name        = model_name
    self.map_inp           = map_inp
    self.params            = params
    self.out               = out
    self.map_name          = map_name
    self.logfile           = logfile
    self.output_dir        = output_dir
  
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
    
    # Initialize states accumulator
    states = mmtbx.utils.states(
     pdb_hierarchy  = self.model.get_hierarchy(),
     xray_structure = self.model.get_xray_structure())
    states.add(sites_cart = self.model.get_xray_structure().sites_cart())
  
    params = sa.master_params().extract()
    
    params.start_temperature       = self.params.start_temperature
    params.final_temperature       = self.params.final_temperature
    params.cool_rate               = self.params.cool_rate
    params.number_of_steps         = self.params.number_of_steps
    print ("self.params.resolution:", self.params.resolution)
    params.update_grads_shift      = 0.
    params.interleave_minimization = False #Pavel will fix the error that occur when params.interleave_minimization=True
    
    # because of params = sa.master_params().extract() above, map_weight and secondary_structure_enabled are dealt without "params"
    map_weight = self.params.map_weight
    ss_restraints = self.params.pdb_interpretation.secondary_structure.enabled
    remove_outlier_ss_restraints = self.params.pdb_interpretation.secondary_structure.protein.remove_outliers
    
    cryo_fit2_input_command = "phenix.cryo_fit2 " + self.model_name + " " + self.map_name + " " \
                              + "resolution=" + str(self.params.resolution) + " " \
                              + "map_weight=" + str(round(map_weight,1)) + " " \
                              + "start_temperature=" + str(params.start_temperature) + " " \
                              + "final_temperature=" + str(params.final_temperature) + " " \
                              + "cool_rate=" + str(params.cool_rate) + " " \
                              + "number_of_steps=" + str(params.number_of_steps) + " " \
                              + "secondary_structure.enabled=" + str(ss_restraints) + " " \
                              + "secondary_structure.protein.remove_outliers=" + str(remove_outlier_ss_restraints) + "\n"
    print ("cryo_fit2_input_command:",cryo_fit2_input_command)
    
    input_command_file = open("cryo_fit2.input_command.txt", "w")
    input_command_file.write(str(cryo_fit2_input_command))
    input_command_file.close()
    
    self.logfile.write("Input command: ")
    self.logfile.write(str(cryo_fit2_input_command))
    
    cc = round(calculate_cc(map_data=map_data, model=self.model, resolution=self.params.resolution), 3)
    initial_CC = "\nInitial CC: " + str(cc) + "\n"
    
    print('%s' %(initial_CC))
    self.logfile.write(str(initial_CC))
    
    if (self.params.progress_on_screen == True):
      result = sa.run(
        params = params,
        xray_structure     = self.model.get_xray_structure(),
        restraints_manager = self.model.get_restraints_manager(),
        target_map         = map_data,
        real_space         = True,
        wx                 = map_weight, 
        wc                 = 1, # weight for geometry conformation
        states_collector   = states)
    else: # (self.params.progress_on_screen = False):
      result = sa.run(
        params = params,
        xray_structure     = self.model.get_xray_structure(),
        restraints_manager = self.model.get_restraints_manager(),
        target_map         = map_data,
        real_space         = True,
        wx                 = map_weight, 
        wc                 = 1, # weight for geometry conformation
        states_collector   = states,
        log                = self.logfile) # if this is commented, temp= xx dist_moved= xx angles= xx bonds= xx is shown on screen rather than cryo_fit2.log
        
    cc = round(calculate_cc(map_data=map_data, model=self.model, resolution=self.params.resolution), 3)
    final_CC = "Final   CC: " + str(cc) + "\n"
    output_dir_w_CC = str(self.output_dir) + "_CC_" + str(cc)
    if os.path.exists(output_dir_w_CC):
      shutil.rmtree(output_dir_w_CC)
    os.mkdir(output_dir_w_CC)
    
    print('%s' %(final_CC))
    self.logfile.write(str(final_CC))
    
    all_state_file = os.path.join(output_dir_w_CC, "all_states.pdb")
    states.write(file_name = all_state_file)
    
    self.model.set_xray_structure(result.xray_structure)
    
    splited_model_name = self.model_name[:-4].split("/")
    model_file_name_only = splited_model_name[len(splited_model_name)-1] 
    fitted_file_name = model_file_name_only + "_cryo_fit2_fitted.pdb"
    print ("fitted_file_name:", fitted_file_name)
    fitted_file = os.path.join(output_dir_w_CC, fitted_file_name)
    
    with open(fitted_file, "w") as f:
      f.write(self.model.model_as_pdb())
      returned = know_how_much_map_origin_moved(str(self.map_name))
    
    if (returned != "origin_is_all_zero" and self.params.keep_origin == True):
        write_this = "Restoring original position for output model\n"
        print (write_this)
        self.logfile.write(str(write_this))
        return_to_origin_of_pdb_file(fitted_file, returned[0], returned[1], returned[2], returned[3])
    
    return output_dir_w_CC
############# end of run function
