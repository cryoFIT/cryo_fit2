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


def add_extracted_CRYST1_to_pdb_file(self,unit_cell_parameters_from_map):
    write_this_CRYST1 = "CRYST1"
    unit_cell_parameters_from_map = str(unit_cell_parameters_from_map)
    print ("unit_cell_parameters_from_map:",unit_cell_parameters_from_map)
    splited = unit_cell_parameters_from_map.split(",")
    # ref: https://www.wwpdb.org/documentation/file-format-content/format33/sect8.html
    #print ("splited:",splited)
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
    #print ("len(splited_alpha):",len(splited_alpha))
    #print ("len(splited_alpha[0]):",len(splited_alpha[0]))
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
        
    write_this_CRYST1 =  write_this_CRYST1 + "  P 1              # added by cryo_fit2 according to user cryo-EM map\n" # if I added "# added by cryo_fit2" at the end, it complains "iotbx.pdb.records.FormatError: Corrupt Z value:"
    print ("correct_CRYST1 format : CRYST1   40.000   80.000   72.000  90.00  90.00  90.00 P 1")
    print ("write_this_CRYST1     :",write_this_CRYST1)
    
    user_s_original_pdb_file = self.data_manager.get_default_model_name() + ".original"
    command = "cp " + self.data_manager.get_default_model_name() + " " + user_s_original_pdb_file
    libtbx.easy_run.call(command)
    
    line_prepender(self.data_manager.get_default_model_name(), write_this_CRYST1)
    return user_s_original_pdb_file
########################### end of add_extracted_CRYST1_to_pdb_file


def check_whether_args_has_eff(args):
  for i in range(len(args)):
    if args[i][len(args[i])-4:len(args[i])] == ".eff":
      return True
  return False
######## end of check_whether_args_has_eff(args)


def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)
#################### end of line_prepender()


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


def determine_optimal_weight_by_template(self, map_inp):
    pi  = get_pdb_inputs_by_pdb_file_name(self, map_inp)
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

def get_pdb_inputs_by_pdb_file_name(self, map_inp):
    #ppf = mmtbx.utils.process_pdb_file_srv(log=null_out()).process_pdb_files(
    #        pdb_file_names=[self.data_manager.get_default_model_name()])[0]
    #'''
    try: # works if pdb file has CRYST1 and has no atoms with unknown nonbonded energy type symbols
        ppf = mmtbx.utils.process_pdb_file_srv(log=null_out()).process_pdb_files(
            pdb_file_names=[self.data_manager.get_default_model_name()])[0]
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
            unit_cell_parameters_from_map = map_inp.unit_cell_crystal_symmetry().unit_cell()
            #print (unit_cell_parameters_from_map)
            user_s_original_pdb_file = add_extracted_CRYST1_to_pdb_file(self,unit_cell_parameters_from_map)
            
            ppf = mmtbx.utils.process_pdb_file_srv(log=null_out()).process_pdb_files(
                pdb_file_names=[self.data_manager.get_default_model_name()])[0]
            
            # return to original pdb file since I can't guarantee space group and z number
            command = "mv " + user_s_original_pdb_file + " " + self.data_manager.get_default_model_name()
            libtbx.easy_run.call(command)
    
        except:
            print ("map_weight can't be optimized automatically.")
            print ("\n(possible reason 1)        There could be some atoms with unknown nonbonded energy type symbols in the given atomic model.")
            print ("(possible reason 1 solution) Fix atoms with unknown nonbonded energy type symbols in the given atomic model.")
            print ("(possible reason 1 solution) real_space_refine will tell you which atoms have unknown nonbonded energy type symbols in the given atomic model.")
            print ("\n(possible reason 2)        Both pdb file and map file lack CRYST1 information.")
            print ("(possible reason 2 solution) Either add CRYST1 info into .pdb/.cif file, or rerun cryo_fit2 with map_weight.")
            print ("(possible reason 2 solution) For example, phenix.cryo_fit2 model.pdb map.ccp4 resolution=4 map_weight=5")
            print ("(possible reason 2 solution) However, human entered map_weight may not be optimal, e.g. it may break the geometry or may not be enough to fit into cryo-EM map fully.")
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
    print ("#### Know how much map origin moved ####")
    
    # Compute a target map
    from iotbx import ccp4_map
    ccp4_map = ccp4_map.map_reader(map_file_name)
    print ("\tMap read from", str(map_file_name))
    target_map_data = ccp4_map.map_data()
    
    #print ("\tdir(ccp4_map): ", dir(ccp4_map)) #['__class__', '__delattr__', '__dict__', '__doc__', '__format__', '__getattribute__', '__hash__', '__init__', '__module__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'crystal_symmetry', 'data', 'dummy', 'grid_unit_cell', 'header_max', 'header_mean', 'header_min', 'header_rms', 'is_similar_map', 'map_data', 'pixel_sizes', 'show_summary', 'space_group_number', 'statistics', 'unit_cell', 'unit_cell_crystal_symmetry', 'unit_cell_grid', 'unit_cell_parameters']
    
    acc = target_map_data.accessor() # keep for now
    #print ("\tacc =", acc) # it shows object address "<scitbx_array_family_flex_ext.grid object at 0x110fe3780>"
    #print ("\tdir(acc): ", dir(acc)) # ['__call__', '__class__', '__delattr__', '__dict__', '__doc__', '__eq__', '__format__', '__getattribute__', '__getinitargs__', '__getstate__', '__hash__', '__init__', '__instance_size__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__safe_for_unpickling__', '__setattr__', '__setstate__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'all', 'focus', 'focus_size_1d', 'is_0_based', 'is_padded', 'is_trivial_1d', 'is_valid_index', 'last', 'nd', 'origin', 'set_focus', 'shift_origin', 'show_summary', 'size_1d']
    
    acc_all = target_map_data.accessor().all() # keep for now
    print ("\tacc_all =", acc_all) # (112, 102, 75) for L1 stalk
    #print ("\tdir(acc_all): ", dir(acc_all)) # ['__add__', '__class__', '__contains__', '__delattr__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__getitem__', '__getnewargs__', '__getslice__', '__gt__', '__hash__', '__init__', '__iter__', '__le__', '__len__', '__lt__', '__mul__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__rmul__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', 'count', 'index']
    
    print ("\ttarget_map_data.origin():", str(target_map_data.origin())) # (120, 76, 202) for L1 stalk
    #STOP()
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


            
              
def rewrite_to_custom_geometry(user_input_pymol_ss):
  f_in = open(user_input_pymol_ss)
  out_file = user_input_pymol_ss[:-4] + '_custom_geom.eff'
  f_out = open(out_file, "w")
  f_out.write('''geometry_restraints {
  edits {
''')
  for line in f_in:
    dist_candidate = line[0:4]
    if (dist_candidate == "dist"):
      #print line
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
########## end of rewrite function


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





#dist chain "A" and resi   42  and name  N4  and alt '', chain "A" and resi   28  and name  O6  and alt ''
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