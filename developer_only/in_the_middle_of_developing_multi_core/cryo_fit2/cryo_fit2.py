from __future__ import division, print_function
import iotbx.phil, libtbx
from libtbx import group_args
from libtbx.utils import Sorry
from iotbx import map_and_model
import mmtbx.utils, os
from mmtbx.dynamics import simulated_annealing as sa

def calculate_cc(map_data, model, resolution):
    xrs = model.get_xray_structure()
    fc = xrs.structure_factors(d_min = resolution).f_calc()
    f_map = fc.structure_factors_from_map(
      map            = map_data,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)
    return fc.map_correlation(other = f_map)
####################### end of calculate_cc function

def know_how_much_map_origin_moved(map_file_name):
    print ("Know map origin")
    
    # Compute a target map
    from iotbx import ccp4_map
    ccp4_map = ccp4_map.map_reader(map_file_name)
    print ("\tMap read from", str(map_file_name))
    target_map_data = ccp4_map.map_data()
    
    #print "\tdir(): ", dir(ccp4_map)
    # acc = target_map_data.accessor() # not used, but keep for now
    print ("\ttarget_map_data.origin():", str(target_map_data.origin()))

    emmap_x0 = target_map_data.origin()[0] # tRNA: 0, nucleosome: -98    
    emmap_y0 = target_map_data.origin()[1] # tRNA: 0, nucleosome: -98
    emmap_z0 = target_map_data.origin()[2] # tRNA: 0, nucleosome: -98
    
    if (emmap_x0 == 0 and emmap_y0 == 0 and emmap_z0 == 0):
        return "origin_is_all_zero"
    else:
        print ("\t\tccp4_map.unit_cell_parameters", str(ccp4_map.unit_cell_parameters))
        a,b,c = ccp4_map.unit_cell_parameters[:3]
        widthx = a/target_map_data.all()[0]
        print ("\t\twidthx", str(widthx)) # with nucleosome, I confirmed that widthx doesn't change by origin shift
        
        shifted_in_x = target_map_data.origin()[0]
        shifted_in_y = target_map_data.origin()[1]
        shifted_in_z = target_map_data.origin()[2]
        print ("shifted_in_x", str(shifted_in_x))
        print ("shifted_in_y", str(shifted_in_y))
        print ("shifted_in_z", str(shifted_in_z))
        return widthx, shifted_in_x, shifted_in_y, shifted_in_z     
####################### end of know_how_much_map_origin_moved function

def return_to_origin_of_pdb_file(input_pdb_file_name, widthx, move_x_by, move_y_by, move_z_by):
    move_x_by = move_x_by*widthx
    move_y_by = move_y_by*widthx
    move_z_by = move_z_by*widthx
    print ("move_x_by:",move_x_by)
    print ("move_y_by:",move_y_by)
    print ("move_z_by:",move_z_by)
    f_in = open(input_pdb_file_name)
    
    output_pdb_file_name = input_pdb_file_name[:-4] + "_origin_returned" + ".pdb"
    
    f_out = open(output_pdb_file_name, "w")
    for line in f_in:
      if line[0:4] == "ATOM" or line[0:6] == "HETATM":
        x_coor_former = line[30:38]
        new_x_coor = str(float(x_coor_former) + float(move_x_by))
        new_x_coor = str(round(float(new_x_coor), 3))
        splited = new_x_coor.split(".")
        multi_before_period = 4-len(splited[0])
        multi_after_period = 3-len(splited[1])
        new_line = line[:30] + multi_before_period*" "+splited[0] + "." + splited [1]+multi_after_period*" "
        
        y_coor_former = line[38:46]
        new_y_coor = str(float(y_coor_former) + float(move_y_by))
        new_y_coor = str(round(float(new_y_coor), 3))
        splited = new_y_coor.split(".")
        multi_before_period = 4-len(splited[0])
        multi_after_period = 3-len(splited[1])
        new_line = new_line + multi_before_period*" "+splited[0] + "." + splited [1]+multi_after_period*" "
        
        z_coor_former = line[46:54]
        new_z_coor = str(float(z_coor_former) + float(move_z_by))      
        new_z_coor = str(round(float(new_z_coor), 3))
        splited = new_z_coor.split(".")
        multi_before_period = 4-len(splited[0])
        multi_after_period = 3-len(splited[1])
        new_line = new_line + multi_before_period*" "+splited[0] + "." + splited [1]+multi_after_period*" " \
              + line[54:]
        f_out.write(new_line)
        
      elif line[0:3] == "TER":
        f_out.write(line)
    f_in.close()
    f_out.close()
    command = "mv " + output_pdb_file_name + " " + input_pdb_file_name
    libtbx.easy_run.call(command)
####################### end of return_to_origin_of_pdb_file ()

class cryo_fit2(object):
  def __init__(self, model, model_name, map_inp, params, out, map_name, logfile):
    self.model             = model
    self.model_name        = model_name
    self.map_inp           = map_inp
    self.params            = params
    self.out               = out
    self.map_name          = map_name
    self.logfile           = logfile

  def validate(self): # this functions runs
    assert not None in [self.model, self.params, self.out]
    if (self.model is None):
      raise Sorry("Model is required.")
    if (self.map_inp is None):
      raise Sorry("Map is required.")
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
    
    params.start_temperature = self.params.start_temperature
    params.final_temperature = self.params.final_temperature
    params.cool_rate = self.params.cool_rate
    params.number_of_steps = self.params.number_of_steps
    params.update_grads_shift = 0.
    params.interleave_minimization = False #Pavel will fix the error that occur when params.interleave_minimization=True
    #wx = self.params.wx # expected wx value doesn't work with multi_core_run()
    params.wx = self.params.wx
    
    output_dir = ''
    if (self.params.output_dir is not None) :
        output_dir = self.params.output_dir
    else:
        output_dir = "."

    cryo_fit2_input_command = "phenix.cryo_fit2 " + self.model_name + " " + self.map_name + " " \
                              + "start_temperature=" + str(params.start_temperature) + " " \
                              + "final_temperature=" + str(params.final_temperature) + " " \
                              + "cool_rate=" + str(params.cool_rate) + " " \
                              + "number_of_steps=" + str(params.number_of_steps) + " " \
                              + "wx=" + str(params.wx) + "\n"
    print ("cryo_fit2_input_command:",cryo_fit2_input_command)
    
    
    input_command_file = open("cryo_fit2.input_command.txt", "w")
    input_command_file.write(str(cryo_fit2_input_command))
    input_command_file.close()
    
    self.logfile.write("Input command: ")
    self.logfile.write(str(cryo_fit2_input_command))
    
    cc = round(calculate_cc(map_data=map_data, model=self.model, resolution=3.), 3)
    
    initial_CC = "\ninitial CC: " + str(cc) + "\n"
    STOP()
    print('%s' %(initial_CC))
    self.logfile.write(str(initial_CC))
    
    result = sa.run(
      params = params,
      xray_structure     = self.model.get_xray_structure(),
      restraints_manager = self.model.get_restraints_manager(),
      target_map         = map_data,
      real_space         = True,
      wx                 = params.wx, # weight for cryo-EM map, wx=5 broke helix conformation of tst_00_poor.pdb, wx=100 kept helix well
      wc                 = 1, # weight for stereochemistry/correct conformation
      states_collector   = states)
    
    cc = round(calculate_cc(map_data=map_data, model=self.model, resolution=3.), 3)
    final_CC = "final   CC: " + str(cc) + "\n\n"
    
    print('%s' %(final_CC))
    self.logfile.write(str(final_CC))
    
    if os.path.exists(output_dir):
        print ("")
    else:
        os.makedirs(output_dir)
    
    ''' # temporarily disabled to save time
    all_state_file = os.path.join(output_dir, "all_states.pdb")
    states.write(file_name = all_state_file)
    '''
    
    self.model.set_xray_structure(result.xray_structure)
    
    fitted_file_name = "cryo_fitted_number_of_steps_" + str(params.number_of_steps) + "_wx_" + str(params.wx) + ".pdb"
    fitted_file = os.path.join(output_dir, fitted_file_name)
    
    with open(fitted_file, "w") as f:
      f.write(self.model.model_as_pdb())
      returned = know_how_much_map_origin_moved(str(self.map_name))
    #print ("returned",returned)
    #print ("self.params.keep_origin",self.params.keep_origin)
    
    if (returned != "origin_is_all_zero" and self.params.keep_origin == True):
        write_this = "Restoring original position for output model\n"
        print (write_this)
        self.logfile.write(str(write_this))
        return_to_origin_of_pdb_file(fitted_file, returned[0], returned[1], returned[2], returned[3])
############# end of run function