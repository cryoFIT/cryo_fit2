from __future__ import division, print_function
import iotbx.phil, libtbx
from libtbx import group_args
from libtbx.utils import Sorry
from iotbx import map_and_model
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


def determine_optimal_weight(d_min, pdb_str_1): # followed tst_weight.py
    print (d_min, "-"*69)
    pi = get_pdb_inputs(pdb_str=pdb_str_1)
    f_calc = pi.xrs.structure_factors(d_min = d_min).f_calc()
    fft_map = f_calc.fft_map(resolution_factor=0.25)
    fft_map.apply_sigma_scaling()
    map_data = fft_map.real_map_unpadded()
    w = weight.run(
      map_data                    = map_data,
      xray_structure              = pi.xrs,
      pdb_hierarchy               = pi.ph,
      geometry_restraints_manager = pi.grm).weight
    return w
######################## end of determine_optimal_weight
    

def get_pdb_inputs(pdb_str):
  ppf = mmtbx.utils.process_pdb_file_srv().process_pdb_files(
    raw_records=pdb_str.splitlines())[0]
  xrs = ppf.xray_structure(show_summary = False)
  restraints_manager = mmtbx.restraints.manager(
    geometry      = ppf.geometry_restraints_manager(show_energies = False),
    normalization = True)
  return group_args(
    ph  = ppf.all_chain_proxies.pdb_hierarchy,
    grm = restraints_manager,
    xrs = xrs)
######################## end of get_pdb_inputs


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
