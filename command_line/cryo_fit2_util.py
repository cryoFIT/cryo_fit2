from __future__ import division, print_function
from cctbx import xray
import iotbx.phil, libtbx
from iotbx import map_and_model
from iotbx.xplor import crystal_symmetry_from_map

from libtbx import group_args
from libtbx.utils import Sorry
from libtbx.utils import null_out

import mmtbx.utils, os
from mmtbx.dynamics import simulated_annealing as sa
from mmtbx.refinement.real_space import weight

import shutil


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

#'''
def check_whether_first_line_starts_w_CRYST1(pdb_file):
    fo = open(pdb_file, "r")
    lines = fo.readlines()
    for line in lines:
        print ("line:",line)
        if line[0:6] == "CRYST1":
            fo.close()
            return True
    fo.close()
    return False
####################### end of check_whether_first_line_starts_w_CRYST1 function
#'''

def count_ATOM_HETATM(pdb_file):
    number_of_ATOM_HETATM = 0
    fo = open(pdb_file, "r")
    lines = fo.readlines()
    for line in lines:
        #print ("line:",line)
        if line[0:4] == "ATOM" or line[0:6] == "HETATM":
            number_of_ATOM_HETATM = number_of_ATOM_HETATM + 1
    fo.close()
    return number_of_ATOM_HETATM
####################### end of count_ATOM_HETATM function


def determine_optimal_weight_by_template(self):
    pi  = get_pdb_inputs_by_pdb_file_name(self)
    f_calc = pi.xrs.structure_factors(d_min = self.params.resolution).f_calc()
    fft_map = f_calc.fft_map(resolution_factor=0.25)
    fft_map.apply_sigma_scaling()
    map_data = fft_map.real_map_unpadded()
    self.params.map_weight = weight.run(
      map_data                    = map_data,
      xray_structure              = pi.xrs,
      pdb_hierarchy               = pi.ph,
      geometry_restraints_manager = pi.grm).weight

    print ("An optimized weight for a map", str(self.params.map_weight))
    return self.params.map_weight
######################### end of determine_optimal_weight_by_template


def determine_optimal_weight_as_macro_cycle_RSR(self, map_inp, model_inp):
    self.structure_monitor = mmtbx.refinement.real_space.rsr_model(
      model             = model_inp,
      target_map_object = map_inp.map_data())
    
    tmp_xrs = self.structure_monitor.model.get_xray_structure().deep_copy_scatterers()
    self.params.map_weight = weight.run(
        map_data                    = map_inp.map_data,
        xray_structure              = tmp_xrs,
        pdb_hierarchy               = model_inp.get_hierarchy(),
        geometry_restraints_manager = model_inp.get_restraints_manager())#,
        #rms_bonds_limit             = self.params.refinement.target_bonds_rmsd,
        #rms_angles_limit            = self.params.refinement.target_angles_rmsd,
        #ncs_groups                  = self.ncs_groups)
    
    print ("An optimized weight for a map", str(self.params.map_weight))
    return self.params.map_weight
####################### end of determine_optimal_weight_as_macro_cycle_RSR()


def get_pdb_inputs_by_pdb_file_name(self):
    try: # works if pdb file has CRYST1
        ppf = mmtbx.utils.process_pdb_file_srv(log=null_out()).process_pdb_files(
        pdb_file_names=[self.data_manager.get_default_model_name()])[0]
    except: # above try results in "Sorry: Crystal symmetry is missing or cannot be extracted."
        try: # works if pdb file has CRYST1
            pdb_inp = iotbx.pdb.input(file_name=self.data_manager.get_default_model_name())
            model = mmtbx.model.manager(model_input = pdb_inp, crystal_symmetry = crystal_symmetry_from_map)
            # tRNA and GAC result in "AttributeError: 'module' object has no attribute '_unit_cell'"
        except:
            print ("\nBoth pdb file and map file lack CRYST1 information.")
            print ("Therefore, map_weight can't be determined automatically.")
            print ("Either add CRYST1 info into .pdb/.cif file, or rerun cryo_fit2 with map_weight.")
            print ("For example, phenix.cryo_fit2 model.pdb map.ccp4 resolution=4 map_weight=5")
            exit(1)
  
    xrs = ppf.xray_structure(show_summary = False)
    restraints_manager = mmtbx.restraints.manager(
      geometry      = ppf.geometry_restraints_manager(show_energies = False),
      normalization = True)
    return group_args(
      ph  = ppf.all_chain_proxies.pdb_hierarchy,
      grm = restraints_manager,
      xrs = xrs)
######################## end of get_pdb_inputs_by_pdb_file_name


def know_how_much_map_origin_moved(map_file_name):
    print ("Know how much map origin moved")
    
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
        return widthx, shifted_in_x, shifted_in_y, shifted_in_z     
############## end of know_how_much_map_origin_moved function


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
################################## end of return_to_origin_of_pdb_file ()


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


