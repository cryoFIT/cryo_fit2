from __future__ import division, print_function
from cctbx.uctbx import unit_cell
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


def check_whether_args_has_eff(args):
  for i in range(len(args)):
    if args[i][len(args[i])-4:len(args[i])] == ".eff":
      return True
  return False
######## end of check_whether_args_has_eff(args)


def check_whether_the_pdb_file_has_nucleic_acid(pdb_file):
    fo = open(pdb_file, "r")
    lines = fo.readlines()
    for line in lines:
        #print ("line:",line)
        residue = line[17:20].strip()
        if (residue == "A") or (residue == "U") or (residue == "G") or (residue == "C") or (residue == "T"):
            fo.close()
            return True
    fo.close()
    return False
####################### end of check_whether_the_pdb_file_has_nucleic_acid()



def determine_optimal_weight_by_template(self, logfile, map_inp, current_fitted_file, weight_boost):
  pi = get_pdb_inputs_by_pdb_file_name(self, logfile, map_inp, current_fitted_file)
  f_calc = pi.xrs.structure_factors(d_min = self.params.resolution).f_calc()
  fft_map = f_calc.fft_map(resolution_factor=0.25)
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  self.params.map_weight = weight.run(
    map_data                    = map_data,
    xray_structure              = pi.xrs,
    pdb_hierarchy               = pi.ph,
    geometry_restraints_manager = pi.grm).weight

  #return self.params.map_weight # 1x~10x of weight_boost were not enough for L1 stalk fitting
  return weight_boost*self.params.map_weight # up to 20x of weight_boost, nucleic acid geometry was ok, 30x broke it
######################### end of determine_optimal_weight_by_template


'''
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
######################## end of determine_optimal_weight_as_macro_cycle_RSR()
'''


def get_pdb_inputs_by_pdb_file_name(self, logfile, map_inp, current_fitted_file):
  
  try: # works if pdb file has CRYST1 and has no atoms with unknown nonbonded energy type symbols and resolution is correctly assigned

      ppf = ''
      try:
        ppf = mmtbx.utils.process_pdb_file_srv(log=null_out()).process_pdb_files(
          pdb_file_names=[self.data_manager.get_default_model_name()])[0]
      except:
        ppf = mmtbx.utils.process_pdb_file_srv(log=null_out()).process_pdb_files(
          pdb_file_names=[current_fitted_file])[0]
      
  except:
      # above try either results in "Sorry: Crystal symmetry is missing or cannot be extracted."
      # or
      #   "Sorry: Fatal problems interpreting model file:
      #    Number of atoms with unknown nonbonded energy type symbols: xx
      #    Please edit the model file to resolve the problems and/or supply a
      #    CIF file with matching restraint definitions, along with
      #    apply_cif_modification and apply_cif_link parameter definitions
      #    if necessary."
      #
      try: # try to extract CRYST1 info from map
          write_this = "\nCRYST1 info is not extracted from user input pdb file. Try to extract it from user map instead.\n"
          print (write_this)
          logfile.write(write_this)
          
          file_name_w_user_s_original_pdb_info = prepend_extracted_CRYST1_to_pdb_file(self, logfile, map_inp)
          
          ppf = mmtbx.utils.process_pdb_file_srv(log=null_out()).process_pdb_files(
              pdb_file_names=[self.data_manager.get_default_model_name()])[0]
          
      except:
          write_this = '''
=============================================
map_weight can't be optimized automatically.
=============================================

[Possible reason 1]  User entered wrong resolution. When 4 angstrom resolution was entered for a 9 ansgtrom resolution map, cryo_fit2 can't optimize cryo_em_map_weight.
[Solution for this]  Enter a correct resolution. A user can get the resolution either by EMDB reported value or by running phenix.mtriage


[Possible reason 2]  There could be some residues/atoms with unknown nonbonded energy type symbols in the given atomic model.
[Solution for this]  Fix atoms with unknown nonbonded energy type symbols in the given atomic model.
                     real_space_refine in PHENIX GUI will show users which atoms have unknown nonbonded energy type symbols.
[Note]               cryo_fit2 cleans archaic \"RX\" type nucleic acid name from user input pbd file automatically.
                     If a user manually assigns a map_weight like map_weight=x, then cryo_fit2 will run, but with a message of "Number of atoms with unknown nonbonded energy type symbols:"


[Possible reason 3]  If a user input pdb file lacks CRYST1 header info (https://www.wwpdb.org/documentation/file-format-content/format33/sect8.html)
                     cryo_fit2 automatically assigns it from map to the first line of user input pdb file.
                     However, when a user provided .cif input model file rather than .pdb input model file, this automatic assign of CRYST1 doesn't work.
                     Additionally, when both user input pdb file and map file lack CRYST1 information, cryo_fit2 can't optimize map_weight.
[Solution 1]         Add CRYST1 header information manually into .pdb/.cif file and rerun cryo_fit2
[Solution 2]         Rerun cryo_fit2 with user specified map_weight.
                     For example, phenix.cryo_fit2 model.pdb map.ccp4 resolution=4 map_weight=5
                     However, human entered map_weight may not be optimal, e.g. it may break the geometry or may not be enough to fit into cryo-EM map fully.
[Note]               phenix.cif_as_pdb and phenix.pdb_as_cif interconvert between .pdb format and .cif format 
'''
          print (write_this)
          logfile.write(write_this)
          exit(1)
  
  xrs = ppf.xray_structure(show_summary = False)
  restraints_manager = mmtbx.restraints.manager(
    geometry      = ppf.geometry_restraints_manager(show_energies = False),
    normalization = True)
  
  return group_args(
          ph  = ppf.all_chain_proxies.pdb_hierarchy,
          grm = restraints_manager,
          xrs = xrs)
######################## end of get_pdb_inputs_by_pdb_file_name function



def know_how_much_map_origin_moved(map_file_name):
    print ("#### Know how much map origin moved ####")
    
    # Compute a target map
    from iotbx import ccp4_map
    ccp4_map = ccp4_map.map_reader(map_file_name)
    print ("\tMap read from", str(map_file_name))
    ccp4_map_data = ccp4_map.map_data()
    
    # try
    #print ("\tccp4_map.unit_cell_grid[0]:", ccp4_map.unit_cell_grid[0])
  
    #print ("\tdir(ccp4_map_data): ", dir(ccp4_map_data)) #dir(ccp4_map_data):  ['__abs__', '__add__', '__class__', '__delattr__', '__delitem__', '__dict__', '__div__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__getitem__', '__getitem_fgdit__', '__getstate__', '__gt__', '__hash__', '__iadd__', '__idiv__', '__imul__', '__init__', '__instance_size__', '__isub__', '__itruediv__', '__le__', '__len__', '__lt__', '__module__', '__mul__', '__ne__', '__neg__', '__new__', '__pow__', '__radd__', '__rdiv__', '__reduce__', '__reduce_ex__', '__repr__', '__rmul__', '__rsub__', '__rtruediv__', '__safe_for_unpickling__', '__setattr__', '__setitem__', '__setstate__', '__sizeof__', '__str__', '__sub__', '__subclasshook__', '__truediv__', '__weakref__', 'accessor', 'add_selected', 'all', 'all_approx_equal', 'all_approx_equal_relatively', 'all_eq', 'all_ge', 'all_gt', 'all_le', 'all_lt', 'all_ne', 'angle', 'append', 'as_1d', 'as_double', 'as_float', 'as_numpy_array', 'as_scitbx_matrix', 'as_string', 'assign', 'back', 'capacity', 'clear', 'concatenate', 'copy_selected', 'copy_to_byte_str', 'cos_angle', 'count', 'deep_copy', 'dot', 'eight_point_interpolation', 'eight_point_interpolation_with_gradients', 'element_size', 'extend', 'fill', 'focus', 'focus_size_1d', 'format_max', 'format_mean', 'format_min', 'front', 'id', 'insert', 'iround', 'is_0_based', 'is_padded', 'is_square_matrix', 'is_trivial_1d', 'last', 'mathematica_form', 'matrix_back_substitution', 'matrix_back_substitution_given_transpose', 'matrix_copy_block', 'matrix_copy_column', 'matrix_copy_lower_to_upper_triangle_in_place', 'matrix_copy_lower_triangle', 'matrix_copy_upper_to_lower_triangle_in_place', 'matrix_copy_upper_triangle', 'matrix_determinant_via_lu', 'matrix_diagonal', 'matrix_diagonal_add_in_place', 'matrix_diagonal_product', 'matrix_diagonal_set_in_place', 'matrix_diagonal_sum', 'matrix_forward_substitution', 'matrix_forward_substitution_given_transpose', 'matrix_inversion', 'matrix_inversion_in_place', 'matrix_is_symmetric', 'matrix_lower_bidiagonal', 'matrix_lower_triangle_as_packed_l', 'matrix_lu_back_substitution', 'matrix_lu_decomposition_in_place', 'matrix_multiply', 'matrix_multiply_packed_u', 'matrix_multiply_packed_u_multiply_lhs_transpose', 'matrix_multiply_transpose', 'matrix_norm_1', 'matrix_norm_frobenius', 'matrix_norm_inf', 'matrix_outer_product', 'matrix_packed_l_as_lower_triangle', 'matrix_packed_l_as_symmetric', 'matrix_packed_u_as_symmetric', 'matrix_packed_u_as_upper_triangle', 'matrix_packed_u_diagonal', 'matrix_packed_u_diagonal_add_in_place', 'matrix_packed_u_swap_rows_and_columns_in_place', 'matrix_paste_block_in_place', 'matrix_paste_column_in_place', 'matrix_rot90', 'matrix_swap_columns_in_place', 'matrix_swap_rows_in_place', 'matrix_symmetric_as_packed_l', 'matrix_symmetric_as_packed_u', 'matrix_symmetric_upper_triangle_quadratic_form', 'matrix_symmetric_upper_triangle_swap_rows_and_columns_in_place', 'matrix_trace', 'matrix_transpose', 'matrix_transpose_in_place', 'matrix_transpose_multiply', 'matrix_transpose_multiply_as_packed_u', 'matrix_transpose_multiply_diagonal_multiply_as_packed_u', 'matrix_upper_bidiagonal', 'matrix_upper_triangle_as_packed_u', 'min_max_mean', 'nd', 'norm', 'norm_1', 'norm_inf', 'origin', 'pop_back', 'quadratic_interpolation_with_gradients', 'reserve', 'reshape', 'resize', 'reversed', 'round', 'sample_standard_deviation', 'select', 'set_selected', 'shallow_copy', 'shift_origin', 'size', 'slice_to_byte_str', 'standard_deviation_of_the_sample', 'tricubic_interpolation', 'tricubic_interpolation_with_gradients', 'value_at_closest_grid_point']
    
    #print ("\tdir(ccp4_map): ", dir(ccp4_map)) #['__class__', '__delattr__', '__dict__', '__doc__', '__format__', '__getattribute__', '__hash__', '__init__', '__module__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'crystal_symmetry', 'data', 'dummy', 'grid_unit_cell', 'header_max', 'header_mean', 'header_min', 'header_rms', 'is_similar_map', 'map_data', 'pixel_sizes', 'show_summary', 'space_group_number', 'statistics', 'unit_cell', 'unit_cell_crystal_symmetry', 'unit_cell_grid', 'unit_cell_parameters']
    #print ("\tdir(ccp4_map_data): ", dir(ccp4_map_data)) # ['__abs__', '__add__', '__class__', '__delattr__', '__delitem__', '__dict__', '__div__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__getitem__', '__getitem_fgdit__', '__getstate__', '__gt__', '__hash__', '__iadd__', '__idiv__', '__imul__', '__init__', '__instance_size__', '__isub__', '__itruediv__', '__le__', '__len__', '__lt__', '__module__', '__mul__', '__ne__', '__neg__', '__new__', '__pow__', '__radd__', '__rdiv__', '__reduce__', '__reduce_ex__', '__repr__', '__rmul__', '__rsub__', '__rtruediv__', '__safe_for_unpickling__', '__setattr__', '__setitem__', '__setstate__', '__sizeof__', '__str__', '__sub__', '__subclasshook__', '__truediv__', '__weakref__', 'accessor', 'add_selected', 'all', 'all_approx_equal', 'all_approx_equal_relatively', 'all_eq', 'all_ge', 'all_gt', 'all_le', 'all_lt', 'all_ne', 'angle', 'append', 'as_1d', 'as_double', 'as_float', 'as_numpy_array', 'as_scitbx_matrix', 'as_string', 'assign', 'back', 'capacity', 'clear', 'concatenate', 'copy_selected', 'copy_to_byte_str', 'cos_angle', 'count', 'deep_copy', 'dot', 'eight_point_interpolation', 'eight_point_interpolation_with_gradients', 'element_size', 'extend', 'fill', 'focus', 'focus_size_1d', 'format_max', 'format_mean', 'format_min', 'front', 'id', 'insert', 'iround', 'is_0_based', 'is_padded', 'is_square_matrix', 'is_trivial_1d', 'last', 'mathematica_form', 'matrix_back_substitution', 'matrix_back_substitution_given_transpose', 'matrix_copy_block', 'matrix_copy_column', 'matrix_copy_lower_to_upper_triangle_in_place', 'matrix_copy_lower_triangle', 'matrix_copy_upper_to_lower_triangle_in_place', 'matrix_copy_upper_triangle', 'matrix_determinant_via_lu', 'matrix_diagonal', 'matrix_diagonal_add_in_place', 'matrix_diagonal_product', 'matrix_diagonal_set_in_place', 'matrix_diagonal_sum', 'matrix_forward_substitution', 'matrix_forward_substitution_given_transpose', 'matrix_inversion', 'matrix_inversion_in_place', 'matrix_is_symmetric', 'matrix_lower_bidiagonal', 'matrix_lower_triangle_as_packed_l', 'matrix_lu_back_substitution', 'matrix_lu_decomposition_in_place', 'matrix_multiply', 'matrix_multiply_packed_u', 'matrix_multiply_packed_u_multiply_lhs_transpose', 'matrix_multiply_transpose', 'matrix_norm_1', 'matrix_norm_frobenius', 'matrix_norm_inf', 'matrix_outer_product', 'matrix_packed_l_as_lower_triangle', 'matrix_packed_l_as_symmetric', 'matrix_packed_u_as_symmetric', 'matrix_packed_u_as_upper_triangle', 'matrix_packed_u_diagonal', 'matrix_packed_u_diagonal_add_in_place', 'matrix_packed_u_swap_rows_and_columns_in_place', 'matrix_paste_block_in_place', 'matrix_paste_column_in_place', 'matrix_rot90', 'matrix_swap_columns_in_place', 'matrix_swap_rows_in_place', 'matrix_symmetric_as_packed_l', 'matrix_symmetric_as_packed_u', 'matrix_symmetric_upper_triangle_quadratic_form', 'matrix_symmetric_upper_triangle_swap_rows_and_columns_in_place', 'matrix_trace', 'matrix_transpose', 'matrix_transpose_in_place', 'matrix_transpose_multiply', 'matrix_transpose_multiply_as_packed_u', 'matrix_transpose_multiply_diagonal_multiply_as_packed_u', 'matrix_upper_bidiagonal', 'matrix_upper_triangle_as_packed_u', 'min_max_mean', 'nd', 'norm', 'norm_1', 'norm_inf', 'origin', 'pop_back', 'quadratic_interpolation_with_gradients', 'reserve', 'reshape', 'resize', 'reversed', 'round', 'sample_standard_deviation', 'select', 'set_selected', 'shallow_copy', 'shift_origin', 'size', 'slice_to_byte_str', 'standard_deviation_of_the_sample', 'tricubic_interpolation', 'tricubic_interpolation_with_gradients', 'value_at_closest_grid_point']
    
    acc = ccp4_map_data.accessor() # keep for now
    #print ("\tacc =", acc) # it shows object address "<scitbx_array_family_flex_ext.grid object at 0x110fe3780>"
    #print ("\tdir(acc): ", dir(acc)) # ['__call__', '__class__', '__delattr__', '__dict__', '__doc__', '__eq__', '__format__', '__getattribute__', '__getinitargs__', '__getstate__', '__hash__', '__init__', '__instance_size__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__safe_for_unpickling__', '__setattr__', '__setstate__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'all', 'focus', 'focus_size_1d', 'is_0_based', 'is_padded', 'is_trivial_1d', 'is_valid_index', 'last', 'nd', 'origin', 'set_focus', 'shift_origin', 'show_summary', 'size_1d']
    
    acc_all = ccp4_map_data.accessor().all() # keep for now
    print ("\tacc_all =", acc_all)
    # L1 stalk original emd_6315 (0, 0, 0)
    # L1 stalk boxed map (by DN) (99, 87, 85)
    
    #print ("\tdir(acc_all): ", dir(acc_all)) # ['__add__', '__class__', '__contains__', '__delattr__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__getitem__', '__getnewargs__', '__getslice__', '__gt__', '__hash__', '__init__', '__iter__', '__le__', '__len__', '__lt__', '__mul__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__rmul__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', 'count', 'index']
    
    print ("\tccp4_map_data.origin():", str(ccp4_map_data.origin()))
    
    emmap_x0 = ccp4_map_data.origin()[0] # origin in x axis, tRNA: 0, nucleosome: -98    
    emmap_y0 = ccp4_map_data.origin()[1] # origin in y axis, tRNA: 0, nucleosome: -98
    emmap_z0 = ccp4_map_data.origin()[2] # origin in z axis, tRNA: 0, nucleosome: -98
    print ("\temmap_x0:",emmap_x0,"emmap_y0:",emmap_y0,"emmap_z0:",emmap_z0)
    # L1 stalk original emd_6315 (0, 0, 0)
    # L1 stalk boxed map (by DN) (132, 94, 203)
    #################### Doonam confirmed that this is an answer to fix origin problem
    
    a,b,c = ccp4_map.unit_cell_parameters[:3]
    print ("\ta:",a,"b:",b,"c:",c)
    # L1 stalk both 'original emd_6315' and 'boxed map (by DN)': (377.9999694824219, 377.9999694824219, 377.9999694824219, 90.0, 90.0, 90.0)
    
    print ("\tccp4_map_data.all():",ccp4_map_data.all()) # these are same as map grid for L1 stalk both original and map_boxed map
    # L1 stalk original emd_6315 (360, 360, 360)
    # L1 stalk boxed map (by DN) (99, 87, 85)
    
    
    ##### this works for most maps except phenix.map_box ed maps
    #widthx = a/ccp4_map_data.all()[0]
    # widthx (1.04999) is same as pixel size (1.05) for L1 stalk map_boxed map
    # widthx (3.82) is much larger than pixel size (1.05) for L1 stalk original map
    
    widthx = a/ccp4_map.unit_cell_grid[0]
    print ("\twidthx", str(widthx)) # with nucleosome, I confirmed that widthx doesn't change by origin shift
    
    if (emmap_x0 == 0 and emmap_y0 == 0 and emmap_z0 == 0):
        return "origin_is_all_zero"
    else:
        shifted_in_x = emmap_x0
        shifted_in_y = emmap_y0
        shifted_in_z = emmap_z0
        return widthx, shifted_in_x, shifted_in_y, shifted_in_z     
############## end of know_how_much_map_origin_moved function


def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)
#################### end of line_prepender()


def prepend_extracted_CRYST1_to_pdb_file(self, logfile, map_inp):
    write_this_CRYST1 = "CRYST1"
    unit_cell_parameters_from_map = str(map_inp.unit_cell_crystal_symmetry().unit_cell())
    splited = unit_cell_parameters_from_map.split(",") # ref: https://www.wwpdb.org/documentation/file-format-content/format33/sect8.html
    
    soon_a = splited[0]
    splited_soon_a = soon_a.split("(")
    a = splited_soon_a[1]
    
    multi_before_period = ''
    multi_after_period = ''
    
    splited_a = a.split(".")
    if (len(splited_a) == 1): # just 336
        multi_before_period = 5-len(splited_a[0])
        multi_after_period  = 3
        write_this_CRYST1 = write_this_CRYST1 + multi_before_period*" "+splited_a[0] + multi_after_period*" "    
    else:
        if (len(splited_a[0]) <= 5):
          multi_before_period = 5-len(splited_a[0])
          multi_after_period = 3-len(splited_a[1])
        else:
          multi_before_period = 7-len(splited_a[0])
          multi_after_period = 0-len(splited_a[1])
        write_this_CRYST1 = write_this_CRYST1 + multi_before_period*" "+splited_a[0] + "." + splited_a[1]+multi_after_period*" "    
    
    
    b = splited[1]
    splited_b = b.split(".")
    if (len(splited_b) == 1): # just 336
        multi_before_period = 6-len(splited_b[0])
        multi_after_period  = 3
        write_this_CRYST1 = write_this_CRYST1 + multi_before_period*" "+splited_b[0] + multi_after_period*" " 
    else:
        if (len(splited_b[0]) <= 5):
            multi_before_period = 5-len(splited_b[0])
            multi_after_period = 3-len(splited_b[1])
        else:
            multi_before_period = 7-len(splited_b[0])
            multi_after_period = 0-len(splited_b[1])
        write_this_CRYST1 = write_this_CRYST1 + multi_before_period*" "+splited_b[0] + "." + splited_b[1]+multi_after_period*" "

    
    c = splited[2]
    splited_c = c.split(".")
    if (len(splited_c) == 1): # just 336
        multi_before_period = 6-len(splited_c[0])
        multi_after_period  = 3
        write_this_CRYST1 = write_this_CRYST1 + multi_before_period*" "+splited_c[0] + multi_after_period*" "
    else:
        if (len(splited_c[0]) <= 5):
            multi_before_period = 5-len(splited_c[0])
            multi_after_period = 3-len(splited_c[1])
        else:
            multi_before_period = 7-len(splited_c[0])
            multi_after_period = 0-len(splited_c[1])
        write_this_CRYST1 = write_this_CRYST1 + multi_before_period*" "+splited_c[0] + "." + splited_c[1]+multi_after_period*" "
    

    alpha = splited[3].strip(' ')
    splited_alpha = alpha.split(".")
    if (len(splited_alpha) == 1): # just 90
        multi_before_period = 4-len(splited_alpha[0])
        multi_after_period  = 2
        write_this_CRYST1 = write_this_CRYST1 + multi_before_period*" "+splited_alpha[0] + multi_after_period*" "
    else:
        if (len(splited_alpha[0]) <= 5):
            multi_before_period = 4-len(splited_alpha[0])
            multi_after_period = 2-len(splited_alpha[1])
        else:
            multi_before_period = 5-len(splited_alpha[0])
            multi_after_period = 0-len(splited_alpha[1])
        write_this_CRYST1 = write_this_CRYST1 + multi_before_period*" "+splited_alpha[0] + "." + splited_alpha[1]+multi_after_period*" "
    
    
    beta = splited[4]
    splited_beta = beta.split(".")
    if (len(splited_beta) == 1): # just 90
        multi_before_period = 5-len(splited_beta[0])
        multi_after_period  = 2
        write_this_CRYST1 = write_this_CRYST1 + multi_before_period*" "+splited_beta[0] + multi_after_period*" "
    else:
        if (len(splited_beta[0]) <= 5):
            multi_before_period = 4-len(splited_beta[0])
            multi_after_period = 2-len(splited_beta[1])
        else:
            multi_before_period = 5-len(splited_beta[0])
            multi_after_period = 0-len(splited_beta[1])
        write_this_CRYST1 = write_this_CRYST1 + multi_before_period*" "+splited_beta[0] + "." + splited_beta[1]+multi_after_period*" "
        
    
    soon_gamma = splited[5]
    splited_soon_gamma = soon_gamma.split(")")
    gamma = splited_soon_gamma[0]
    splited_gamma = gamma.split(".")
    if (len(splited_gamma) == 1): # just 90
        multi_before_period = 5-len(splited_gamma[0])
        multi_after_period  = 2
        write_this_CRYST1 = write_this_CRYST1 + multi_before_period*" "+splited_gamma[0] + multi_after_period*" "
    else:
        if (len(splited_gamma[0]) <= 5):
            multi_before_period = 4-len(splited_gamma[0])
            multi_after_period = 2-len(splited_gamma[1])
        else:
            multi_before_period = 5-len(splited_gamma[0])
            multi_after_period = 0-len(splited_gamma[1])
        write_this_CRYST1 = write_this_CRYST1 + multi_before_period*" "+splited_gamma[0] + "." + splited_gamma[1]+multi_after_period*" "
    
    
    print ("Examplar correct CRYST1 format: CRYST1   40.000   80.000   72.000  90.00  90.00  90.00 P 1")
    
    if (map_inp.space_group_number == 19):
      write_this = "space_group_number from user input map = 19. Therefore, assign P 21 21 21 to a user input pdb file\n" # http://img.chem.ucl.ac.uk/sgp/large/019a.htm
      print (write_this)
      logfile.write(write_this)
      write_this_CRYST1 =  write_this_CRYST1 + "  P 21 21 21       # added by cryo_fit2 according to the user cryo-EM map\n" # if I added "# added by cryo_fit2" at the end, it complains "iotbx.pdb.records.FormatError: Corrupt Z value:"
    else:  
      write_this_CRYST1 =  write_this_CRYST1 + "  P 1              # added by cryo_fit2 according to the user cryo-EM map\n" # if I added "# added by cryo_fit2" at the end, it complains "iotbx.pdb.records.FormatError: Corrupt Z value:"
    
    print ("write this CRYST1 information :",write_this_CRYST1)
    
    file_name_w_user_s_original_pdb_info = self.data_manager.get_default_model_name() + ".original"
    command = "cp " + self.data_manager.get_default_model_name() + " " + file_name_w_user_s_original_pdb_info
    libtbx.easy_run.call(command)
    
    line_prepender(self.data_manager.get_default_model_name(), write_this_CRYST1)
    
    write_this = "CRYST1 info from map was prepended to user pdb file, therefore the very original user pdb file is now renamed to " + file_name_w_user_s_original_pdb_info
    print (write_this)
    logfile.write(write_this)
    
    return file_name_w_user_s_original_pdb_info
############################ end of prepend_extracted_CRYST1_to_pdb_file



def remove_R_prefix_in_RNA(input_pdb_file_name): ######### deal very old style of RNA file
    f_in = open(input_pdb_file_name)
    output_pdb_file_name = input_pdb_file_name[:-4] + "_cleaned_for_cryo_fit2.pdb"
    f_out = open(output_pdb_file_name, 'wt')
    cleaned = False
    for line in f_in:
      if (line[0:6] == "HEADER"):
        f_out.write(line)
      else:
        residue = line[17:20]
        trimmed_residue = residue.replace(" ", "")
        if (trimmed_residue == "RA"):
            cleaned = True
            new_line = line[:17] + " A " + line[20:]
            f_out.write(new_line)
        elif (trimmed_residue == "RT"): 
            cleaned = True
            new_line = line[:17] + " T " + line[20:]
            f_out.write(new_line)
        elif (trimmed_residue == "RG"): 
            cleaned = True
            new_line = line[:17] + " G " + line[20:]
            f_out.write(new_line)
        elif (trimmed_residue == "RC"): 
            cleaned = True
            new_line = line[:17] + " C " + line[20:]
            f_out.write(new_line)
        elif (trimmed_residue == "RU"): 
            cleaned = True
            new_line = line[:17] + " U " + line[20:]
            f_out.write(new_line)
        else:
            f_out.write(line)
    f_in.close()
    f_out.close()
    
    if (cleaned == False):
      cmd = "rm " + output_pdb_file_name
      libtbx.easy_run.call(cmd)
    return cleaned, output_pdb_file_name
########################### end of remove_R_prefix_in_RNA function


def reoptimize_map_weight_if_not_specified(self, user_map_weight, map_inp, weight_boost):
  if (user_map_weight == ''):
      write_this = "User didn't specify map_weight. Therefore, automatically optimize map_weight for additional cryo_fit2 MD run\n"
      print('%s' %(write_this))
      self.logfile.write(str(write_this))

      current_fitted_file_name = "current_fitted_file.pdb"
      with open(current_fitted_file_name, "w") as f:
        f.write(self.model.model_as_pdb())
      f.close()
      
      self.params.map_weight = determine_optimal_weight_by_template(self, self.logfile, map_inp, current_fitted_file_name, weight_boost)
      
      cmd = "rm " + current_fitted_file_name
      libtbx.easy_run.fully_buffered(cmd)
      
      write_this = "Automatically optimized map weight: " + str(round(self.params.map_weight,2)) + "\n"
      print('%s' %(write_this))
      self.logfile.write(write_this)
  else:
      self.params.map_weight = user_map_weight

  return self.params.map_weight
############### end of reoptimize_map_weight_if_not_specified function
            

def return_to_origin_of_pdb_file(input_pdb_file_name, widthx, move_x_by, move_y_by, move_z_by):
    print ("widthx:",widthx)
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


def rewrite_pymol_ss_to_custom_geometry_ss(user_input_pymol_ss):
####### reference

##### [pymol ss] dist chain "A" and resi   42  and name  N4  and alt '', chain "A" and resi   28  and name  O6  and alt ''

##### [custom geometry ss]
  '''
geometry_restraints {
  edits {
    bond {
      atom_selection_1 = chain 'A' and resid 28 and name N1
      atom_selection_2 = chain 'A' and resid 42 and name N3
      distance_ideal = 2.8
      sigma = 0.021
    }
    bond {
      atom_selection_1 = chain 'A' and resid 7 and name N1
      atom_selection_2 = chain 'A' and resid 66 and name N3
      distance_ideal = 2.8
      sigma = 0.021
    }
  }
}
  '''

  f_in = open(user_input_pymol_ss)
  out_file = user_input_pymol_ss[:-4] + '_custom_geom.eff'
  f_out = open(out_file, "w")
  f_out.write('''geometry_restraints {
  edits {
''')
  for line in f_in:
    dist_candidate = line[0:4]
    if (dist_candidate == "dist"):
      splited = line.split()
      
      write_this = "    bond {\n"
      f_out.write(write_this)
      
      chain_candidate = splited[2]
      splited_chain_candidate = chain_candidate.split("\"")
      
      resi1 = splited[5].strip()
      atom1 = splited[8]
      write_this = "      atom_selection_1 = chain \'" + splited_chain_candidate[1] + "\' and resid " + resi1 + " and name " + atom1 + "\n"
      f_out.write(write_this)
      
      chain_candidate = splited[13]
      splited_chain_candidate = chain_candidate.split("\"")
      resi2 = splited[16].strip()
      atom2 = splited[19]
      write_this = "      atom_selection_2 = chain \'" + splited_chain_candidate[1] + "\' and resid " + resi2 + " and name " + atom2 + "\n"
      f_out.write(write_this)
      
      
      if ((atom1 == "N6" and atom2 ==  "O4") or (atom1 == "O4" and atom2 ==  "N6")):
        f_out.write("      distance_ideal = 3.0\n")
      elif ((atom1 == "O6" and atom2 ==  "N4") or (atom1 == "N4" and atom2 ==  "O6")):
        f_out.write("      distance_ideal = 2.93\n")
      elif ((atom1 == "N2" and atom2 ==  "O2") or (atom1 == "O2" and atom2 ==  "N2")):
        f_out.write("      distance_ideal = 2.78\n")
      elif ((atom1 == "N1" and atom2 ==  "N3") or (atom1 == "N3" and atom2 ==  "N1")):
        f_out.write("      distance_ideal = 2.85\n")
        '''
        if ((resi1 == "G" and resi2 ==  "C") or (resi1 == "C" and resi2 ==  "G")):
          f_out.write("      distance_ideal = 2.88\n")
        else:
          f_out.write("      distance_ideal = 2.82\n") # between A and T
        '''
      else: # default
        f_out.write("      distance_ideal = 2.91\n")
      # reference modules/cctbx_project/mmtbx/secondary_structure/nucleic_acids.py
      
      f_out.write("      sigma = 0.021\n")
      
      write_this = "    }\n"
      f_out.write(write_this)
  f_out.write('''  }
}
''')
  f_in.close()
  f_out.close()
########## end of rewrite_pymol_ss_to_custom_geometry_ss function


def show_time(time_start, time_end):
  time_took = 0 # temporary of course
  if (round((time_end-time_start)/60, 1) < 1):
    time_took = " finished in " + str(round((time_end-time_start), 1)) + " seconds (wallclock)."
  elif (round((time_end-time_start)/60/60, 1) < 1):
    time_took = " finished in " + str(round((time_end-time_start)/60, 1)) + " minutes (wallclock)."
  else:
    time_took = " finished in " + str(round((time_end-time_start)/60/60, 1)) + " hours (wallclock)."
  return time_took
############### end of show_time function
