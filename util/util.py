# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_SIGNALS_DEFAULT=1
#### (source) http://cci.lbl.gov/cctbx_sources/crys3d/command_line/hklview.py

from __future__ import division, print_function
from cctbx.uctbx import unit_cell
from cctbx import xray

import cryo_fit2_run

import decimal, glob, os, platform, subprocess
from os import path

import iotbx.phil, libtbx
from iotbx import map_and_model
from iotbx.xplor import crystal_symmetry_from_map

from libtbx import group_args
from libtbx.str_utils import format_value as fv
from libtbx.utils import Sorry
from libtbx.utils import null_out

from mmtbx.dynamics import simulated_annealing as sa
import mmtbx.maps.map_model_cc
from mmtbx.refinement.real_space import weight
from mmtbx.superpose import *
import mmtbx.utils

import numpy as np
import shutil

os.environ['BOOST_ADAPTBX_FPE_DEFAULT'] = "1"
os.environ['BOOST_ADAPTBX_SIGNALS_DEFAULT'] = "1"

def calculate_overall_cc(map_data, model, resolution):
    xrs = model.get_xray_structure()
    fc = xrs.structure_factors(d_min = resolution).f_calc()
    f_map = fc.structure_factors_from_map(
      map            = map_data,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)
    return fc.map_correlation(other = f_map)
####################### end of calculate_overall_cc function


def calculate_RMSD(self, fitted_file_name_w_path): # (reference) cctbx_project/mmtbx/superpose.py
    fixed = self.model_name
    moving = fitted_file_name_w_path
    
    write_this = "\n===== RMSD calculation ====="
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
      write_this = "\n\n===== Aligning %s to %s ====="%(fitted_file_name_w_path, self.model_name)
      print (write_this)
      self.logfile.write(str(write_this))
      if not self.params.selection_moving:
        moving.selectomatic(fixed)
      rmsd, lsq = moving.superpose(fixed)
      write_this = "\n\nrmsd after cryo_fit2: " + str(round(rmsd,2)) + " angstrom\n"
      print (write_this)
      self.logfile.write(str(write_this))
############ def calculate_RMSD(self):

'''
def check_whether_args_has_eff(args, logfile, location_of_this_code, known_sigma_for_auto_geom):
  for i in range(len(args)):
    if args[i][len(args[i])-4:len(args[i])] == ".eff":
        user_eff_file_name = str(args[i])
        write_this = ''
        if (str(location_of_this_code) == str("prepare_cryo_fit2")):
            write_this = "A user provided " + user_eff_file_name
        else:
            write_this = "cryo_fit2 automatically generated " + user_eff_file_name + " with " + str(known_sigma_for_auto_geom) + " sigma_for_auto_geom"
        write_this = write_this + " that cryo_fit2 will use."
        print (write_this)
        logfile.write(write_this)
        return True, user_eff_file_name
  return False, ""
######## end of check_whether_args_has_eff(args)
'''

def check_whether_the_pdb_file_has_amino_acid(pdb_file):
    file_opened = open(pdb_file, "r")
    lines = file_opened.readlines()
    for line in lines:
        #print ("line:",line)
        residue = line[17:20].strip()
        if (residue == "ALA") or (residue == "CYS") or (residue == "ASP") or (residue == "GLU") \
        or (residue == "PHE") or (residue == "GLY") or (residue == "HIS") or (residue == "ILE") \
        or (residue == "LYS") or (residue == "LEU") or (residue == "MET") or (residue == "ASN") \
        or (residue == "PRO") or (residue == "GLN") or (residue == "ARG") or (residue == "SER") \
        or (residue == "THR") or (residue == "VAL") or (residue == "TRP") or (residue == "TYR"):
            file_opened.close()
            return True
    file_opened.close()
    return False
####################### end of check_whether_the_pdb_file_has_amino_acid()


def check_whether_the_pdb_file_has_nucleic_acid(pdb_file):
    file_opened = open(pdb_file, "r")
    lines = file_opened.readlines()
    for line in lines:
        #print ("line:",line)
        residue = line[17:20].strip()
        if (residue == "A") or (residue == "U") or (residue == "G") or (residue == "C") \
        or (residue == "T"):
            file_opened.close()
            return True
    file_opened.close()
    return False
####################### end of check_whether_the_pdb_file_has_nucleic_acid()



def clean_unusual_residue(input_pdb_file_name):
    report_file = "unusual_residue_removed.txt"
    if (os.path.isfile(report_file) == True):
      os.remove(report_file)
    f_report = open(report_file, 'w')
    
    f_in = open(input_pdb_file_name)
    output_pdb_file_name = input_pdb_file_name[:-4] + "_cleaned_unusual.pdb"
    f_out = open(output_pdb_file_name, 'wt')
    
    cleaned_unusual_residue = False
    
    for line in f_in:
        residue = line[17:20]
        if (line[0:6] == "HETATM"):
            if (residue == "34G"):
                write_this = str(residue) + " removed\n"
                f_report.write(write_this)
                print (write_this)
                cleaned_unusual_residue = True
                continue
            else:
                f_out.write(line)
        else:
          f_out.write(line)
    f_in.close()
    f_out.close()
    f_report.close()
    
    if (cleaned_unusual_residue == False):
      os.remove(output_pdb_file_name)
      os.remove(report_file)
            
    return output_pdb_file_name, cleaned_unusual_residue
############ end of clean_unusual_residue function


def count_bp_H_E_in_fitted_file(fitted_file_name_w_path, output_dir_w_CC, logfile):
    
    if (os.path.isfile(fitted_file_name_w_path) == False):
        write_this = str(fitted_file_name_w_path) + " doesn't exist"
        print (write_this)
        logfile.write(write_this)
        
    starting_dir = os.getcwd()
    os.chdir(output_dir_w_CC)
    
    splited_fitted_file_name_w_path = fitted_file_name_w_path.split("/")
    fitted_file_name_wo_path = splited_fitted_file_name_w_path[len(splited_fitted_file_name_w_path)-1]
    
    command_string = "phenix.secondary_structure_restraints " + fitted_file_name_wo_path
    try:
        libtbx.easy_run.fully_buffered(command_string)
    except:
        write_this = '''
          Maybe the \"fitted\" structure is broken with a following message.
          I confirmed this with both RNA (L1 stalk) and protein (Mg_channel).
          
          /home/doonam/bin/phenix-dev-3311/build/lib/scitbx_array_family_flex_ext.so(boost::python::detail::caller_arity<2u>::impl<scitbx::af::versa<double, scitbx::af::flex_grid<scitbx::af::small<long, 10ul> > > (*)(scitbx::af::versa<double, scitbx::af::flex_grid<scitbx::af::small<long, 10ul> > > const&, double const&), boost::python::default_call_policies, boost::mpl::vector3<scitbx::af::versa<double, scitbx::af::flex_grid<scitbx::af::small<long, 10ul> > >, scitbx::af::versa<double, scitbx::af::flex_grid<scitbx::af::small<long, 10ul> > > const&, double const&> >::operator()(_object*, _object*)+0x123) [0x7f589e621013]
          /home/doonam/bin/phenix-dev-3311/build/lib/scitbx_array_family_flex_ext.so(scitbx::af::boost_python::flex_wrapper<double, boost::python::return_value_policy<boost::python::copy_non_const_reference, boost::python::default_call_policies> >::div_a_s(scitbx::af::versa<double, scitbx::af::flex_grid<scitbx::af::small<long, 10ul> > > const&, double const&)+0x9e) [0x7f589e60cf3e]
          /lib/x86_64-linux-gnu/libc.so.6(+0x354b0) [0x7f58a3f444b0]
          Floating-point error (Python and libc call stacks above)
                This crash may be due to a problem in any imported
                Python module, including modules which are not part
                of the cctbx project. To disable the traps leading
                to this message, define these environment variables
                (e.g. assign the value 1):
                    BOOST_ADAPTBX_FPE_DEFAULT
                    BOOST_ADAPTBX_SIGNALS_DEFAULT
                This will NOT solve the problem, just mask it, but
                may allow you to proceed in case it is not critical.
        '''
        print (write_this)
        logfile.write(write_this)
        return None, None, None
    
    ss_file_name = fitted_file_name_wo_path + "_ss.eff"

    ## count bp
    command_string = "cat " + ss_file_name + " | grep base_pair | wc -l"
    grepped = libtbx.easy_run.fully_buffered(command=command_string).raise_if_errors().stdout_lines
    number_of_bp_in_fitted_pdb = int(grepped[0])
    
    ## count H
    command_string = "cat " + ss_file_name + " | grep \"helix {\" | wc -l"
    grepped = libtbx.easy_run.fully_buffered(command=command_string).raise_if_errors().stdout_lines
    number_of_H_in_fitted_pdb = int(grepped[0])
    
    ## count E
    command_string = "cat " + ss_file_name + " | grep \"sheet {\" | wc -l"
    grepped = libtbx.easy_run.fully_buffered(command=command_string).raise_if_errors().stdout_lines
    number_of_E_in_fitted_pdb = int(grepped[0])
    
    os.chdir(starting_dir)
    
    return number_of_bp_in_fitted_pdb, number_of_H_in_fitted_pdb, number_of_E_in_fitted_pdb
######################## end of count_bp_H_E_in_fitted_file(fitted_file_name_w_path):


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

# this fn runs only once in entire cryo_fit2 running (final MD may run once more, but I'm disabling it)
def determine_optimal_weight_by_template(self, logfile, map_inp, current_fitted_file):
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

  return self.params.map_weight 
######################### end of determine_optimal_weight_by_template


def explore_parameters_by_multi_core(self, params, logfile, user_map_weight, bp_cutoff, H_cutoff, E_cutoff, \
                                     MD_in_each_cycle, number_of_steps, start_temperature, \
                                     weight_multiply):
    #print ("logfile:", str(logfile)) # logfile: <open file 'cryo_fit2.log', mode 'w' at 0x11ac73300>
    
    print ("\nMD parameters that will be explored.")
    print ("MD_in_each_cycle:        ", str(MD_in_each_cycle))
    print ("number_of_steps:         ", str(number_of_steps))
    print ("start_temperature:       ", str(start_temperature))
    print ("weight_multiply:         ", str(weight_multiply), "\n\n")

    print ("params.final_temperature:                       ", str(params.final_temperature))
    print ("params.map_weight:                              ", str(round(params.map_weight,2)))

    if (("tst_cryo_fit2" in self.data_manager.get_default_model_name()) == True):
        params.max_steps_for_exploration = 30 # temporary for development
    
    model_inp = self.data_manager.get_model()
    map_inp   = self.data_manager.get_real_map()

    # Re-assign params for below cryo_fit2 run
    params.MD_in_each_cycle         = MD_in_each_cycle
    params.number_of_steps          = number_of_steps
    params.start_temperature        = start_temperature
    params.weight_multiply          = weight_multiply
    
    params.cool_rate = float((float(params.start_temperature)-float(params.final_temperature))/(int(params.MD_in_each_cycle)-1))
    print ("params.cool_rate:                               ", str(round(params.cool_rate, 1)))
    
    init_output_dir = get_output_dir_name(self)
    
    task_obj = cryo_fit2_run.cryo_fit2_class(
      model             = model_inp,
      model_name        = self.data_manager.get_default_model_name(),
      map_inp           = map_inp,
      params            = self.params,
      out               = self.logger,
      map_name          = self.data_manager.get_default_real_map_name(),
      logfile           = logfile,
      output_dir        = init_output_dir,
      user_map_weight   = user_map_weight,
      weight_multiply   = self.params.weight_multiply)
    
    task_obj.validate()
    
    # Wrapping below task_obj.run() with try/except will not report number of parameter exploration combinations that ran successfully.
    # However, I expected that this may help incomplete running issue.
    # However, it turned out that this try/except loop didn't help incomplete running issue.
    
    # ((((( during explore_parameters_by_multi_core)))))
    # even wrapping try-except around map-weighted phenix.dynamics itself didn't help that the output folder didn't move to "parameters_exploration" folder
    output_dir_final = ''
    try:
        output_dir_final = task_obj.run()
    except:
        
    # (below is not working, so commented out)
    #except Exception as ex:
    #    write_this = "exception message:" +  str(ex)
    #    print (write_this)
    #    logfile.write(str(write_this))

        # this write_this reported appropriately -> really ? -> yes really
        write_this = "(after task_obj loop) An exception occurred in explore_parameters_by_multi_core. Maybe cryo_fit2 failed to run (\"nan\") for this condition:" + \
                     " cool_rate (" + str(round(params.cool_rate, 1))   + ")" + \
                     " MD_in_each_cycle (" + str(params.MD_in_each_cycle)      + ")" + \
                     " number_of_steps (" + str(params.number_of_steps)        + ")" + \
                     " start_temperature (" + str(params.start_temperature)    + ")" + \
                     " weight_multiply (" + str(self.params.weight_multiply)        + ")" + \
                     " final_temperature (" + str(params.final_temperature)    + ")" + \
                     " map_weight (" + str(round(params.map_weight,2))  + ")" + \
                     " max_steps_for_exploration (" + str(params.max_steps_for_exploration)  + ")" 
        print (write_this)
        logfile.write(str(write_this))        
    
    if (output_dir_final != ''):
        if (output_dir_final.find('_bp_') == -1): # here output_dir_final could be ''
            if (os.path.isdir("parameters_exploration/bp_H_E_not_calculated") == False):
                os.mkdir("parameters_exploration/bp_H_E_not_calculated")
            command_string = "mv " + str(output_dir_final) + " parameters_exploration/bp_H_E_not_calculated"
            logfile.write(str(command_string))
            libtbx.easy_run.fully_buffered(command=command_string).raise_if_errors().stdout_lines
            return None, None, None
    else:
        return None, None, None

    #### Regular finish
    splited = output_dir_final.split("_bp_")
    splited2 = splited[1].split("_H_")
    bp = splited2[0]

    splited = output_dir_final.split("_H_")
    splited2 = splited[1].split("_E_")
    H = splited2[0]

    splited = output_dir_final.split("_E_")
    E = splited[len(splited)-1]

    if ( (float(bp) >= float(bp_cutoff)) and (float(H) >= float(H_cutoff)) and (float(E) >= float(E_cutoff))):
        if (os.path.isdir("parameters_exploration/bp_H_E_kept") == False):
            os.mkdir("parameters_exploration/bp_H_E_kept")
        command_string = "mv " + str(output_dir_final) + " parameters_exploration/bp_H_E_kept"
        libtbx.easy_run.fully_buffered(command=command_string).raise_if_errors().stdout_lines
    else:
        if (os.path.isdir("parameters_exploration/bp_H_E_not_kept") == False):
            os.mkdir("parameters_exploration/bp_H_E_not_kept")
        command_string = "mv " + str(output_dir_final) + " parameters_exploration/bp_H_E_not_kept"
        libtbx.easy_run.fully_buffered(command=command_string).raise_if_errors().stdout_lines
    
    return bp, H, E
############ end of explore_parameters_by_multi_core()


def extract_the_best_cc_parameters(logfile):
    starting_dir = os.getcwd()
    if (os.path.isdir("parameters_exploration/bp_H_E_kept") == True):
        os.chdir("parameters_exploration/bp_H_E_kept")
    else:
        write_this = '''MD parameter exploration didn't find parameter combination that kept \
user base pairs, helices and sheets in this user given model file.\n\
Maybe this input pdb has no base_pairs, helices and sheets in the first place?\n\
Otherwise, expand explore combination.\n\
Otherwise, run cryo_fit2 with explore=False\n'''
        logfile.write(write_this)
        print (write_this)
        exit(1)
    
    ## extract the best_cc
    best_cc_so_far = -999
    for check_this_dir in glob.glob("output*"):
        splited = check_this_dir.split("_cc_")
        splited2 = splited[1].split("_bp_")
        cc = splited2[len(splited2)-2]
        if (float(cc) > float(best_cc_so_far)):
            best_cc_so_far = cc
    
    for check_this_dir in glob.glob("output*"):
        splited = check_this_dir.split("_cc_")
        splited2 = splited[1].split("_bp_")
        cc = splited2[len(splited2)-2]
        if (float(cc) == float(best_cc_so_far)):
            splited = check_this_dir.split("_MD_in_each_cycle_")
            splited2 = splited[1].split("_step_")
            optimum_MD_in_each_cycle = splited2[0]
            
            splited = check_this_dir.split("_start_")
            splited2 = splited[1].split("_final_")
            optimum_start_temperature = splited2[0]
            
            splited = check_this_dir.split("_step_")
            splited2 = splited[1].split("_make_ss_for_stronger_ss_")
            optimum_step = splited2[0]
            
            splited = check_this_dir.split("_weight_multiply_")
            splited2 = splited[1].split("_sigma_for_auto_geom_")
            optimum_weight_multiply = splited2[0]
            
            os.chdir(starting_dir)
            
            return optimum_MD_in_each_cycle, optimum_start_temperature, \
                   optimum_step, optimum_weight_multiply
############ end of def extract_the_best_cc_parameters():


def file_size(fname):
    statinfo = os.stat(fname)
    return statinfo.st_size
######## end of file_size(fname)


def float_to_str(f):  # 0.00001 resulted in scientific notation which seems to not function properly, so avoid it
  ctx = decimal.Context() # create a new context for this task

  ctx.prec = 20 # 20 digits should be enough for everyone :D
  """
  Convert the given float to a string,
  without resorting to scientific notation
  """
  d1 = ctx.create_decimal(repr(f))
  return format(d1, 'f')
########## end of def float_to_str(f)


def get_output_dir_name(self):
    # rename output_dir
    output_dir = str(self.params.output_dir) + \
                 "_resolution_" + str(self.params.resolution) + \
                 "_start_" + str(self.params.start_temperature) + \
                 "_final_" + str(self.params.final_temperature) + \
                 "_MD_in_each_cycle_" + str(self.params.MD_in_each_cycle) + \
                 "_step_" + str(self.params.number_of_steps) + \
                 "_make_ss_for_stronger_ss_" + str(self.params.make_ss_for_stronger_ss) + \
                 "_weight_multiply_" + str(round(self.params.weight_multiply,1)) + \
                 "_sigma_for_auto_geom_" + str(self.params.sigma_for_auto_geom) + \
                 "_slack_for_auto_geom_" + str(self.params.slack_for_auto_geom)
                 #"_ss_" + str(self.params.pdb_interpretation.secondary_structure.enabled) + \
                 #"_del_outlier_ss_" + str(self.params.pdb_interpretation.secondary_structure.protein.remove_outliers) + \
                 #"_NA_" + str(self.params.pdb_interpretation.secondary_structure.nucleic_acid.enabled) + \
                 #"_hb_dis_" + str(self.params.pdb_interpretation.secondary_structure.nucleic_acid.hbond_distance_cutoff) + \
                 #"_angle_" + str(self.params.pdb_interpretation.secondary_structure.nucleic_acid.angle_between_bond_and_nucleobase_cutoff)
                 #"_bp_planar_" + str(self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_planarity) + \
                 #"_bp_hb_" + str(self.params.pdb_interpretation.secondary_structure.nucleic_acid.base_pair.restrain_hbonds)
    return output_dir
########## end of def get_output_dir_name(self)


def get_pdb_inputs_by_pdb_file_name(self, logfile, map_inp, current_fitted_file):
  
  try: #  Assigns ppf here well if
       #    input pdb file has CRYST1
       #    and input pdb file has no atoms with unknown nonbonded energy type symbols
       #    and resolution is correctly entered
       #  Therefore, new CRYST1 header will not be added to the input pdb file

      ppf = ''
      try:
        ppf = mmtbx.utils.process_pdb_file_srv(log=null_out()).process_pdb_files(
          pdb_file_names=[self.data_manager.get_default_model_name()])[0]
      except:
        ppf = mmtbx.utils.process_pdb_file_srv(log=null_out()).process_pdb_files(
          pdb_file_names=[current_fitted_file])[0]
      
  except:
      # above try results in
      # either "Sorry: Crystal symmetry is missing or cannot be extracted."
      # or
      #   "Sorry: Fatal problems interpreting model file:
      #    Number of atoms with unknown nonbonded energy type symbols: xx
      #    Please edit the model file to resolve the problems and/or supply a
      #    CIF file with matching restraint definitions, along with
      #    apply_cif_modification and apply_cif_link parameter definitions
      #    if necessary."
      #
      try: # try to extract CRYST1 info from map
          write_this = "\nCRYST1 info is not extracted from user input pdb file. Therefore, cryo_fit2 will try to extract it from user map instead.\n"
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


def get_used_map_weight(file):
    f_in = open(file)
    for line in f_in:
        used_map_weight = line
        f_in.close()
        print ("used_map_weight:",used_map_weight)
        return used_map_weight
      #splited = line.split(" ")
      #cc = splited[4]
################ end of get_users_cc(cc_record)


def hasNumbers(inputString):
    return any(char.isdigit() for char in inputString)
####################### end of hasNumbers fn


def know_bp_H_E_in_a_user_pdb_file(user_pdb_file, logfile):    
    splited_user_pdb_file_w_path = user_pdb_file.split("/")
    user_pdb_file_wo_path = splited_user_pdb_file_w_path[len(splited_user_pdb_file_w_path)-1]
    ss_file_name = user_pdb_file_wo_path + "_ss.eff"
    
    if (os.path.isfile(ss_file_name) == False): # if strong_sigma == False, ss_file may not exist
        command_string = "phenix.secondary_structure_restraints " + user_pdb_file
        libtbx.easy_run.fully_buffered(command_string)
    
    user_pdb_file_path = splited_user_pdb_file_w_path[len(splited_user_pdb_file_w_path)-2]
    
    command_string = "cat " + ss_file_name + " | grep base_pair | wc -l"
    grepped = libtbx.easy_run.fully_buffered(command=command_string).raise_if_errors().stdout_lines
    number_of_bp_in_pdb_file = int(grepped[0])
    
    command_string = "cat " + ss_file_name + " | grep \"helix {\" | wc -l"
    grepped = libtbx.easy_run.fully_buffered(command=command_string).raise_if_errors().stdout_lines
    number_of_H_in_pdb_file = int(grepped[0])
    
    command_string = "cat " + ss_file_name + " | grep \"sheet {\" | wc -l"
    grepped = libtbx.easy_run.fully_buffered(command=command_string).raise_if_errors().stdout_lines
    number_of_E_in_pdb_file = int(grepped[0])
    
    return number_of_bp_in_pdb_file, number_of_H_in_pdb_file, number_of_E_in_pdb_file, ss_file_name
######################## end of def know_bp_H_in_a_user_pdb_file(user_pdb_file)


def know_how_much_map_origin_moved(map_file_name):
    #print ("#### Know how much a map origin moved ####")
    
    # Compute a target map
    from iotbx import ccp4_map
    ccp4_map = ccp4_map.map_reader(map_file_name)
    #print ("\tMap read from", str(map_file_name))
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
    #print ("\tacc_all =", acc_all)
    # L1 stalk original emd_6315 (0, 0, 0)
    # L1 stalk boxed map (by DN) (99, 87, 85)
    
    #print ("\tdir(acc_all): ", dir(acc_all)) # ['__add__', '__class__', '__contains__', '__delattr__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__getitem__', '__getnewargs__', '__getslice__', '__gt__', '__hash__', '__init__', '__iter__', '__le__', '__len__', '__lt__', '__mul__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__rmul__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', 'count', 'index']
    
    #print ("\tccp4_map_data.origin():", str(ccp4_map_data.origin()))
    
    emmap_x0 = ccp4_map_data.origin()[0] # origin in x axis, tRNA: 0, nucleosome: -98    
    emmap_y0 = ccp4_map_data.origin()[1] # origin in y axis, tRNA: 0, nucleosome: -98
    emmap_z0 = ccp4_map_data.origin()[2] # origin in z axis, tRNA: 0, nucleosome: -98
    #print ("\temmap_x0:",emmap_x0,"emmap_y0:",emmap_y0,"emmap_z0:",emmap_z0)
    # L1 stalk original emd_6315 (0, 0, 0)
    # L1 stalk boxed map (by DN) (132, 94, 203)
    #################### Doonam confirmed that this is an answer to fix origin problem
    
    a,b,c = ccp4_map.unit_cell_parameters[:3]
    #print ("\ta:",a,"b:",b,"c:",c)
    # L1 stalk both 'original emd_6315' and 'boxed map (by DN)': (377.9999694824219, 377.9999694824219, 377.9999694824219, 90.0, 90.0, 90.0)
    
    #print ("\tccp4_map_data.all():",ccp4_map_data.all()) # these are same as map grid for L1 stalk both original and map_boxed map
    # L1 stalk original emd_6315 (360, 360, 360)
    # L1 stalk boxed map (by DN) (99, 87, 85)
    
    
    ##### this works for most maps except phenix.map_box ed maps
    #widthx = a/ccp4_map_data.all()[0]
    # widthx (1.04999) is same as pixel size (1.05) for L1 stalk map_boxed map
    # widthx (3.82) is much larger than pixel size (1.05) for L1 stalk original map
    
    widthx = a/ccp4_map.unit_cell_grid[0]
    #print ("\twidthx", str(widthx)) # with nucleosome, I confirmed that widthx doesn't change by origin shift
    
    if (emmap_x0 == 0 and emmap_y0 == 0 and emmap_z0 == 0):
        return "origin_is_all_zero"
    else:
        shifted_in_x = emmap_x0
        shifted_in_y = emmap_y0
        shifted_in_z = emmap_z0
        return widthx, shifted_in_x, shifted_in_y, shifted_in_z     
############## end of know_how_much_map_origin_moved function


def know_number_of_atoms_in_input_pdb(logfile, starting_pdb):
    command_string = "cat " + starting_pdb + " | grep ATOM | wc -l"
    num_ATOMs = libtbx.easy_run.fully_buffered(command=command_string).raise_if_errors().stdout_lines
    number_of_atoms_in_input_pdb = int(num_ATOMs[0])
    write_this = "A user input pdb file, " + starting_pdb + ", has "+ str(number_of_atoms_in_input_pdb) + " atoms\n"
    print (write_this)
    logfile.write(write_this)
    return number_of_atoms_in_input_pdb
################# end of know_number_of_atoms_in_input_pdb()


def leave_one_conformer(logfile, input_pdb_file_name): ######### deal AARG, BARG
    f_in = open(input_pdb_file_name)
    output_pdb_file_name = input_pdb_file_name[:-4] + "_cleaned_for_cryo_fit2.pdb"
    f_out = open(output_pdb_file_name, 'wt')
    cleaned = False
    
    for line in f_in:
        if ((line[0:4] != "ATOM") and (line[0:6] != "HETATM") and (line[0:6] != "ANISOU")):
            f_out.write(line)
        else:
            if (line[16:17] == "A"):
                cleaned = True
                new_line = line[:16] + " " + line[17:20] + line[20:]
                f_out.write(new_line)
            elif (line[16:17] == "B"):
                cleaned = True
                continue
            else:
                f_out.write(line)    
    f_in.close()
    f_out.close()
    
    if (cleaned == False):
        cmd = "rm " + output_pdb_file_name
        libtbx.easy_run.call(cmd)
        return 0
    else:
        file_name_w_user_s_original_pdb_info = input_pdb_file_name + "_before_leaving_one_conformer"
        command = "cp " + str(input_pdb_file_name) + " " + str(file_name_w_user_s_original_pdb_info)
        libtbx.easy_run.call(command)
        command = "mv " + str(output_pdb_file_name) + " " + str(input_pdb_file_name)
        libtbx.easy_run.call(command)
        
        write_this = "cryo_fit2 leaves one conformer. Therefore, the original user pdb file is now renamed to " + file_name_w_user_s_original_pdb_info
        print (write_this)
        logfile.write(write_this)
    
        return 0
########################### end of leave_one_conformer function


def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)
#################### end of line_prepender()


def make_argstuples(self, logfile, user_map_weight, the_pdb_file_has_nucleic_acid, \
                    bp_cutoff, H_cutoff, E_cutoff):
    total_combi_num = 0
    argstuples = []
    ## final_temperature is fixed as 0
    
    if (("tst_cryo_fit2" in self.data_manager.get_default_model_name()) == True):
        # just explore 2 combinations to save regression time
        for MD_in_each_cycle in range (2, 14, 100): 
            for number_of_steps in range (100, 151, 100):
                for start_temperature in np.arange (300.0, 301.0, 300.0):
                    for weight_multiply in range (1, 10, 6): # 2
                        total_combi_num = total_combi_num + 1
                        argstuples.append([self, self.params, logfile, user_map_weight, \
                                           bp_cutoff, H_cutoff, E_cutoff, MD_in_each_cycle, \
                                           number_of_steps, start_temperature, \
                                           weight_multiply])
    else: # regular running                
        for MD_in_each_cycle in range (2, 23, 10): # 3 (e.g. 2, 12, 22) (minimum should be >=2)
            for number_of_steps in range (1, 501, 200): # 5 (e.g. 1, 101, 201, 301, 401)
                #for start_temperature in np.arange (300.0, 901.0, 300.0): # 3 (e.g. 300, 600, 900)
                for start_temperature in np.arange (300.0, 601.0, 300.0): # 900 seems break James' new_small
                    if (("L1_stalk" in self.data_manager.get_default_model_name()) == True):
                        for weight_multiply in range (1, 402, 25): # L1 stalk with sparse map density nearby fails to keep bp with weight_multiply >= 400
                            total_combi_num = total_combi_num + 1 
                            argstuples.append([self, self.params, logfile, user_map_weight, \
                                            bp_cutoff, H_cutoff, E_cutoff, MD_in_each_cycle, \
                                            number_of_steps, start_temperature, \
                                            weight_multiply])
                    else:
                        if (the_pdb_file_has_nucleic_acid == True):
                            for weight_multiply in range (1, 12802, 800): 
                            #for weight_multiply in range (1, 6402, 400): # (with 0.1 sigma) 6001 generated many bp keeping full_tRNA
                            #for weight_multiply in range (1, 3202, 200): # (with 0.5 sigma) 2801 generated many bp keeping full_tRNA
                            #for weight_multiply in range (1, 1602, 100): # (with 0.5 sigma) 1601 generated many bp keeping full_tRNA
                            #for weight_multiply in range (1, 802, 50):   # (with 0.5 sigma) 751 was the best for Mg_Channel, 801 was the best for tRNA 
                                total_combi_num = total_combi_num + 1 
                                argstuples.append([self, self.params, logfile, user_map_weight, \
                                                bp_cutoff, H_cutoff, E_cutoff, MD_in_each_cycle, \
                                                number_of_steps, start_temperature, \
                                                weight_multiply])
                        else: # for protein
                            for weight_multiply in range (1, 802, 50):
                            #for weight_multiply in range (1, 25602, 1600):
                                # (with 0.5 sigma) old_irina kept all ss even with 6402
                                
                                # (with 0.5 sigma) new_small
                                #                            (with 12801) kept all ss
                                #                            (with 1601) kept all ss during exploration, but broke those during final_MD
                                #                            (with 801) increased cc up to 0.40 
                                
                                total_combi_num = total_combi_num + 1 
                                argstuples.append([self, self.params, logfile, user_map_weight, \
                                                bp_cutoff, H_cutoff, E_cutoff, MD_in_each_cycle, \
                                                number_of_steps, start_temperature, \
                                                weight_multiply])
                            
    print ("total_combi_num:",total_combi_num)
    #STOP()
    return total_combi_num, argstuples
##### end of def make_argstuples(logfile):


def prepend_extracted_CRYST1_to_pdb_file(self, logfile, map_inp):
    write_this_CRYST1 = "CRYST1"
    unit_cell_parameters_from_map = str(map_inp.unit_cell_crystal_symmetry().unit_cell())
    splited = unit_cell_parameters_from_map.split(",") # ref: https://www.wwpdb.org/documentation/file-format-content/format33/sect8.html

    # worked for part_tRNA, full_tRNA, L1_stalk, a helix and Mg_Channel
    soon_a = splited[0]
    splited_soon_a = soon_a.split("(")
    a = splited_soon_a[1]
    
    multi_before_period = ''
    multi_after_period = ''
    
    splited_a = a.split(".")
    if (len(splited_a) == 1): # just 336 without numbers after .
        multi_before_period = 5-len(splited_a[0])
        multi_after_period  = 3
        write_this_CRYST1 = write_this_CRYST1 + multi_before_period*" "+splited_a[0] + multi_after_period*" "    
    else:
        multi_before_period = 4-len(splited_a[0])
        multi_after_period = 4-len(splited_a[1])
        write_this_CRYST1 = write_this_CRYST1 + multi_before_period*" "+splited_a[0] + "." + splited_a[1]+multi_after_period*" "    
    
    
    b = splited[1]
    splited_b = b.split(".")
    if (len(splited_b) == 1): # just 336
        multi_before_period = 6-len(splited_b[0])
        multi_after_period  = 3
        write_this_CRYST1 = write_this_CRYST1 + multi_before_period*" "+splited_b[0] + multi_after_period*" " 
    else:
        multi_before_period = 4-len(splited_b[0])
        multi_after_period = 4-len(splited_b[1])
        write_this_CRYST1 = write_this_CRYST1 + multi_before_period*" "+splited_b[0] + "." + splited_b[1]+multi_after_period*" "

    
    c = splited[2]
    splited_c = c.split(".")
    if (len(splited_c) == 1): # just 336
        multi_before_period = 6-len(splited_c[0])
        multi_after_period  = 3
        write_this_CRYST1 = write_this_CRYST1 + multi_before_period*" "+splited_c[0] + multi_after_period*" "
    else:
        multi_before_period = 4-len(splited_c[0])
        multi_after_period = 4-len(splited_c[1])
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
    
    print ("Examplar correct CRYST1 format                : CRYST1   40.000   80.000   72.000  90.00  90.00  90.00 P 1")
    
    if (map_inp.space_group_number == 19):
      write_this = "space_group_number from user input map = 19. Therefore, assign P 21 21 21 to a user input pdb file\n" # http://img.chem.ucl.ac.uk/sgp/large/019a.htm
      print (write_this)
      logfile.write(write_this)
      write_this_CRYST1 =  write_this_CRYST1 + "  P 21 21 21       # added by cryo_fit2 according to the user cryo-EM map\n" # if I added "# added by cryo_fit2" at the end, it complains "iotbx.pdb.records.FormatError: Corrupt Z value:"
    else:  
      write_this_CRYST1 =  write_this_CRYST1 + "  P 1              # added by cryo_fit2 according to the user cryo-EM map\n" # if I added "# added by cryo_fit2" at the end, it complains "iotbx.pdb.records.FormatError: Corrupt Z value:"
    
    print ("Cryo_fit2 will prepend this CRYST1 information:",write_this_CRYST1, " to a user pdb file")
    
    file_name_w_user_s_original_pdb_info = self.data_manager.get_default_model_name() + "_before_adding_CRYST1"
    command = "cp " + self.data_manager.get_default_model_name() + " " + file_name_w_user_s_original_pdb_info
    libtbx.easy_run.call(command)
    
    line_prepender(self.data_manager.get_default_model_name(), write_this_CRYST1)
    
    write_this = "CRYST1 info from map was prepended to user pdb file. Therefore, the original user pdb file is now renamed to " + file_name_w_user_s_original_pdb_info
    print (write_this)
    logfile.write(write_this)
    
    return file_name_w_user_s_original_pdb_info
############################ end of prepend_extracted_CRYST1_to_pdb_file


def remove_prefix_in_AA_name(input_pdb_file_name): ######### deal AARG, BARG
    f_in = open(input_pdb_file_name)
    output_pdb_file_name = input_pdb_file_name[:-4] + "_cleaned_for_cryo_fit2.pdb"
    f_out = open(output_pdb_file_name, 'wt')
    cleaned = False
    
    for line in f_in:
        if ((line[0:4] != "ATOM") and (line[0:6] != "HETATM") and (line[0:6] != "ANISOU")):
            f_out.write(line)
        else:
            if (line[16:17] != ' '):
                cleaned = True
                new_line = line[:16] + " " + line[17:20] + line[20:]
                f_out.write(new_line)
            else:
                f_out.write(line)    
    f_in.close()
    f_out.close()
    
    if (cleaned == False):
        cmd = "rm " + output_pdb_file_name
        libtbx.easy_run.call(cmd)
        return 0
    else:
        file_name_w_user_s_original_pdb_info = input_pdb_file_name + "_before_removing_prefix"
        command = "cp " + str(input_pdb_file_name) + " " + str(file_name_w_user_s_original_pdb_info)
        libtbx.easy_run.call(command)
        command = "mv " + str(output_pdb_file_name) + " " + str(input_pdb_file_name)
        libtbx.easy_run.call(command)
        
        write_this = "cryo_fit2 leaves one confomer. Therefore, the original user pdb file is now renamed to " + file_name_w_user_s_original_pdb_info
        print (write_this)
        logfile.write(write_this)
    
        return 0
########################### end of remove_prefix_in_AA_name function


def remove_R_prefix_in_RNA(input_pdb_file_name): ######### deal very old style of RNA file
    f_in = open(input_pdb_file_name)
    output_pdb_file_name = input_pdb_file_name[:-4] + "_cleaned_for_cryo_fit2.pdb"
    f_out = open(output_pdb_file_name, 'wt')
    cleaned = False
    for line in f_in:
        if ((line[0:4] != "ATOM") and (line[0:6] != "HETATM")):
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


# this fn runs only during final_MD
def reoptimize_map_weight_if_not_specified(self, user_map_weight, map_inp):
  if (user_map_weight == ''):
      write_this = "\nA user didn't specify a map_weight. Therefore, cryo_fit2 will optimize map_weight for additional MD run\n"
      print('%s' %(write_this))
      self.logfile.write(str(write_this))

      current_fitted_file_name = "current_fitted_file.pdb"
      with open(current_fitted_file_name, "w") as f:
        f.write(self.model.model_as_pdb())
      f.close()
      
      self.params.map_weight = determine_optimal_weight_by_template(self, self.logfile, map_inp, current_fitted_file_name)
      
      cmd = "rm " + current_fitted_file_name
      libtbx.easy_run.fully_buffered(cmd)
      
      write_this = "\nAn automatically optimized map weight (during final MD): " + str(round(self.params.map_weight,2)) + "\n"
      print('%s' %(write_this))
      self.logfile.write(write_this)
  else:
      self.params.map_weight = user_map_weight

  return self.params.map_weight
############### end of reoptimize_map_weight_if_not_specified function


def report_map_model_cc(self, map_inp, model, crystal_symmetry, logfile):
    # reference: modules/phenix/phenix/programs/map_model_cc.py
    
    self.params.map_model_cc.resolution = self.params.resolution
    # without this, self.params.map_model_cc.resolution = None
    
    base = map_and_model.input(
      map_data         = map_inp.map_data(),
      model            = model,
      crystal_symmetry = crystal_symmetry,
      box              = True,
      #ignore_symmetry_conflicts= 
      #                   self.params.map_model_cc.ignore_symmetry_conflicts)
      ignore_symmetry_conflicts = True)
    
    task_obj = mmtbx.maps.map_model_cc.map_model_cc(
      map_data         = base.map_data(),
      pdb_hierarchy    = base.model().get_hierarchy(),
      crystal_symmetry = base.model().crystal_symmetry(),
      params           = self.params.map_model_cc)
    
    task_obj.validate()
    task_obj.run()
    self.results = task_obj.get_results()
    self.results.fsc = mmtbx.maps.mtriage.get_fsc(
      map_data = base.map_data(),
      model    = base.model(),
      params   = self.params.map_model_cc)
    r = self.results
    #
    write_this = "CC_mask  : " + str(round(float(r.cc_mask), 4)) + "\n"
    print (write_this)
    logfile.write(str(write_this))
    
    write_this = "CC_Volume  : " + str(round(float(r.cc_volume), 4)) + "\n"
    print (write_this)
    logfile.write(str(write_this))
    
    write_this = "CC_peaks  : " + str(round(float(r.cc_peaks), 4)) + "\n"
    print (write_this)
    logfile.write(str(write_this))
    
    write_this = "CC_box  : " + str(round(float(r.cc_box), 4)) + "\n"
    print (write_this)
    logfile.write(str(write_this))

########## end of def report_map_model_cc():



def return_list_of_eff_from_args(args):
    list_of_eff = []
    for i in range(len(args)):
      if args[i][len(args[i])-4:len(args[i])] == ".eff":
          list_of_eff.append(str(args[i]))
    return list_of_eff
######## end of return_list_of_eff_from_args(args)


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


def rewrite_pymol_ss_to_custom_geometry_ss(user_input_pymol_ss, sigma_for_auto_geom, slack_for_auto_geom):
####### reference

#################### DISTANCE
##### [pymol ss]
#   dist chain "A" and resi   42  and name  N4  and alt '', chain "A" and resi   28  and name  O6  and alt ''

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
  out_file = user_input_pymol_ss[:-4] + '_cryo_fit2_auto.eff'
  f_out = open(out_file, "w")
  f_out.write('''geometry_restraints {
  edits {
''')
  for line in f_in:
    dist_angle_candidate = line[0:5]
    splited = line.split()
    if (len(splited) > 12):
        if (splited[12] == "'B'"): # use A conformer only for angle
            continue
    if (len(splited) > 22):
        if (splited[22] == "'B'"): # use A conformer only for dist
            continue
    if (len(splited) > 23):
        if (splited[23] == "'B',"): # use A conformer only for angle
            continue
    if (len(splited) > 34):
        if (splited[34] == "'B'"): # use A conformer only for angle
            continue
    if (dist_angle_candidate == "dist "):
      write_this = "    bond {\n"
      f_out.write(write_this)
      
      #cctbx_project/mmtbx/monomer_library/pdb_interpretation.py
      write_this = "      action = *add\n"
      
      # for protein and implicit and explicit 0.05, reflected in used_geom
      
      # for protein and explicit 0.03 and 6, not reflected in used_geom
      # because "Bond-like restraints:" overrides my custom geom
      # these "Bond-like restraints:" are reflected in used_geom even when make_ss_for_stronger_ss=False
      
      # for RNA (0.01~6), reflected in used_geom
      
      #write_this = "      action = add\n" # ??
      #write_this = "      action = change\n"  #exception message:geometry_restraints.edits.bond.action = change not implemented
      #write_this = "      action = delete\n"  #exception message:geometry_restraints.edits.bond.action = delete not implemented
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
      
      ########## for nucleic acids, atoms have numbers like N6, O4
      ########## for proteins, atoms do not have numbers like N, O
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
      else:
        # default H-bond length for nucleic acid base pairs and helix and sheet
          f_out.write("      distance_ideal = 2.911\n") # Doonam observed 2.8, 3.2 in a helix
      ########## [reference] modules/cctbx_project/mmtbx/secondary_structure/nucleic_acids.py
      ########## [reference] https://www.phenix-online.org/documentation/reference/secondary_structure.html#proteins
      
      write_this = "      sigma = " + float_to_str(sigma_for_auto_geom) + "\n" 
      f_out.write(write_this) 
      
      write_this = "      slack = " + float_to_str(slack_for_auto_geom) + "\n" 
      f_out.write(write_this) 

      write_this = "    }\n"
      f_out.write(write_this)
      ############# end of if (dist_angle_candidate == "dist"):
        
######################## ANGLE
 ##### [pymol ss]
# angle a0, chain "A" and resi   42  and name  C4  and alt '', chain "A" and resi   42  and name  N4  and alt '', chain "A" and resi   28  and name  O6  and alt ''

##### [custom geometry ss] http://www.phenix-online.org/pipermail/phenixbb/2014-September/021173.html

#    angle {
#      atom_selection_1 = chain 'A' and resid 42 and name C4
#      atom_selection_2 = chain 'A' and resid 42 and name N4
#      atom_selection_3 = chain 'A' and resid 28 and name O6
#      angle_ideal = 117.3
#      sigma = 0.021
#    }

    elif (dist_angle_candidate == "angle"):
        
        write_this = "    angle {\n"
        f_out.write(write_this)
        
        write_this = "      action = *add\n"
        # for protein (0.03, 0.05, 6), reflected in used_geom
        # for RNA (0.01~6), reflected in used_geom
        f_out.write(write_this)
            
        chain_candidate = splited[3]
        splited_chain_candidate = chain_candidate.split("\"")
        resi1 = splited[6].strip()
        atom1 = splited[9]
        write_this = "      atom_selection_1 = chain \'" + splited_chain_candidate[1] + "\' and resid " + resi1 + " and name " + atom1 + "\n"
        f_out.write(write_this)
      
        chain_candidate = splited[14]
        splited_chain_candidate = chain_candidate.split("\"")
        resi2 = splited[17].strip()
        atom2 = splited[20]
        write_this = "      atom_selection_2 = chain \'" + splited_chain_candidate[1] + "\' and resid " + resi2 + " and name " + atom2 + "\n"
        f_out.write(write_this)
        
        chain_candidate = splited[25]
        splited_chain_candidate = chain_candidate.split("\"")
        resi3 = splited[28].strip()
        atom3 = splited[31]
        write_this = "      atom_selection_3 = chain \'" + splited_chain_candidate[1] + "\' and resid " + resi3 + " and name " + atom3 + "\n"
        f_out.write(write_this)
        
        atom1 = splited[9]
        if (hasNumbers(atom1) == False): #this_line_is_protein
            if (atom1 == "C" and atom2 ==  "O" and atom3 ==  "N"): 
                f_out.write("      angle_ideal = 146.0\n") # observed from DN helix: 145, 148.8
            elif (atom1 == "CA" and atom2 ==  "N" and atom3 ==  "O"): 
                f_out.write("      angle_ideal = 111.0\n") # observed from DN helix: 110.6, 111.3
            else: # should be (atom1 == "C" and atom2 ==  "N" and atom3 ==  "O"): 
                f_out.write("      angle_ideal = 120.0\n") # observed from DN helix: 117.8, 122.9
        else: # this is for nucleic acid
            if (atom1 == "C4" and atom2 ==  "N4" and atom3 ==  "O6"): 
                f_out.write("      angle_ideal = 117.3\n") # derived from Oleg slide and modules/cctbx_project/mmtbx/secondary_structure/nucleic_acids.py
            elif (atom1 == "C6" and atom2 ==  "O6" and atom3 ==  "N4"):
                f_out.write("      angle_ideal = 122.8\n") # derived from Oleg slide and modules/cctbx_project/mmtbx/secondary_structure/nucleic_acids.py
            elif (atom1 == "C2" and atom2 ==  "N3" and atom3 ==  "N1"): 
                f_out.write("      angle_ideal = 119.1\n") # either 119.1 or 116.3 # derived from Oleg slide and modules/cctbx_project/mmtbx/secondary_structure/nucleic_acids.py
            elif (atom1 == "C2" and atom2 ==  "N1" and atom3 ==  "N3"):
                f_out.write("      angle_ideal = 116.3\n") # either 119.1 or 116.3 # derived from Oleg slide and modules/cctbx_project/mmtbx/secondary_structure/nucleic_acids.py
            elif (atom1 == "C2" and atom2 ==  "O2" and atom3 ==  "N2"):
                f_out.write("      angle_ideal = 120.7\n") # derived from Oleg slide and modules/cctbx_project/mmtbx/secondary_structure/nucleic_acids.py
            elif (atom1 == "C2" and atom2 ==  "N2" and atom3 ==  "O2"):
                f_out.write("      angle_ideal = 122.2\n") # derived from Oleg slide and modules/cctbx_project/mmtbx/secondary_structure/nucleic_acids.py
            else: 
                f_out.write("      angle_ideal = 120.0\n") # just my guess
        write_this = "      sigma = " + float_to_str(sigma_for_auto_geom) + "\n" 
        f_out.write(write_this)
        
        # geometry_restraints.edits.angle.slack is not recognized
        # "Sorry: Some PHIL parameters are not recognized by phenix.cryo_fit2."
        #  write_this = "      slack = " + float_to_str(slack_for_auto_geom) + "\n" 
        #  f_out.write(write_this)
        
        write_this = "    }\n"
        f_out.write(write_this)
        ############# end of if (dist_angle_candidate == "angle"):
     
  f_out.write('''  }
}
''')
  f_in.close()
  f_out.close()
  
  return out_file
########## end of rewrite_pymol_ss_to_custom_geometry_ss function


def secondary_structure_restraints_DN(args, params=None, out=sys.stdout, log=sys.stderr):

  # parse command-line arguments
  if (params is None):
    pcl = iotbx.phil.process_command_line_with_files(
      args=args,
      master_phil_string=master_phil_str,
      pdb_file_def="file_name")
    work_params = pcl.work.extract()
  # or use parameters defined by GUI
  else:
    work_params = params
  pdb_files = work_params.file_name

  work_params.secondary_structure.enabled=True
  assert work_params.format in ["phenix", "phenix_refine", "phenix_bonds",
      "pymol", "refmac", "kinemage", "pdb"]
  if work_params.quiet :
    out = cStringIO.StringIO()

  pdb_combined = iotbx.pdb.combine_unique_pdb_files(file_names=pdb_files)
  pdb_structure = iotbx.pdb.input(source_info=None,
    lines=flex.std_string(pdb_combined.raw_records))
  cs = pdb_structure.crystal_symmetry()

  corrupted_cs = False
  if cs is not None:
    if [cs.unit_cell(), cs.space_group()].count(None) > 0:
      corrupted_cs = True
      cs = None
    elif cs.unit_cell().volume() < 10:
      corrupted_cs = True
      cs = None

  if cs is None:
    if corrupted_cs:
      print >> out, "Symmetry information is corrupted, "
    else:
      print >> out, "Symmetry information was not found, "
    print >> out, "putting molecule in P1 box."
    from cctbx import uctbx
    atoms = pdb_structure.atoms()
    box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
      sites_cart=atoms.extract_xyz(),
      buffer_layer=3)
    atoms.set_xyz(new_xyz=box.sites_cart)
    cs = box.crystal_symmetry()

  defpars = mmtbx.model.manager.get_default_pdb_interpretation_params()
  defpars.pdb_interpretation.automatic_linking.link_carbohydrates=False
  defpars.pdb_interpretation.c_beta_restraints=False
  defpars.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
  model = mmtbx.model.manager(
      model_input=pdb_structure,
      crystal_symmetry=cs,
      pdb_interpretation_params=defpars,
      stop_for_unknowns=False)
  pdb_hierarchy = model.get_hierarchy()
  geometry = model.get_restraints_manager().geometry
  if len(pdb_hierarchy.models()) != 1 :
    raise Sorry("Multiple models not supported.")
  ss_from_file = None
  if (hasattr(pdb_structure, "extract_secondary_structure") and
      not work_params.ignore_annotation_in_file):
    ss_from_file = pdb_structure.extract_secondary_structure()
  m = manager(pdb_hierarchy=pdb_hierarchy,
    geometry_restraints_manager=geometry,
    sec_str_from_pdb_file=ss_from_file,
    params=work_params.secondary_structure,
    verbose=work_params.verbose)

  # bp_p = nucleic_acids.get_basepair_plane_proxies(
  #     pdb_hierarchy,
  #     m.params.secondary_structure.nucleic_acid.base_pair,
  #     geometry)
  # st_p = nucleic_acids.get_stacking_proxies(
  #     pdb_hierarchy,
  #     m.params.secondary_structure.nucleic_acid.stacking_pair,
  #     geometry)
  # hb_b, hb_a = nucleic_acids.get_basepair_hbond_proxies(pdb_hierarchy,
  #     m.params.secondary_structure.nucleic_acid.base_pair)
  result_out = cStringIO.StringIO()
  # prefix_scope="refinement.pdb_interpretation"
  # prefix_scope=""
  prefix_scope=""
  if work_params.format == "phenix_refine":
    prefix_scope = "refinement.pdb_interpretation"
  elif work_params.format == "phenix":
    prefix_scope = "pdb_interpretation"
  ss_phil = None
  working_phil = m.as_phil_str(master_phil=sec_str_master_phil)
  phil_diff = sec_str_master_phil.fetch_diff(source=working_phil)

  if work_params.format in ["phenix", "phenix_refine"]:
    comment = "\n".join([
      "# These parameters are suitable for use in e.g. phenix.real_space_refine",
      "# or geometry_minimization. To use them in phenix.refine add ",
      "# 'refinement.' if front of pdb_interpretation."])
    if work_params.format == "phenix_refine":
      comment = "\n".join([
      "# These parameters are suitable for use in phenix.refine only.",
      "# To use them in other Phenix tools remove ",
      "# 'refinement.' if front of pdb_interpretation."])
    print >> result_out, comment
    if (prefix_scope != ""):
      print >> result_out, "%s {" % prefix_scope
    if work_params.show_all_params :
      working_phil.show(prefix="  ", out=result_out)
    else :
      phil_diff.show(prefix="  ", out=result_out)
    if (prefix_scope != ""):
      print >> result_out, "}"
  elif work_params.format == "pdb":
    print >> result_out, m.actual_sec_str.as_pdb_str()
  elif work_params.format == "phenix_bonds" :
    raise Sorry("Not yet implemented.")
  elif work_params.format in ["pymol", "refmac", "kinemage"] :
    m.show_summary(log=out)
    (hb_proxies, hb_angle_proxies, planarity_proxies,
        parallelity_proxies) = m.create_all_new_restraints(
        pdb_hierarchy=pdb_hierarchy,
        grm=geometry,
        log=out)
    if hb_proxies.size() > 0:
      if work_params.format == "pymol" :
        file_load_add = "load %s" % work_params.file_name[0]
        # surprisingly, pymol handles filenames with whitespaces without quotes...
        print >> result_out, file_load_add
        bonds_in_format = hb_proxies.as_pymol_dashes(
            pdb_hierarchy=pdb_hierarchy)
      elif work_params.format == "kinemage" :
        bonds_in_format = hb_proxies.as_kinemage(
            pdb_hierarchy=pdb_hierarchy)
      else :
        bonds_in_format = hb_proxies.as_refmac_restraints(
            pdb_hierarchy=pdb_hierarchy)
      print >> result_out, bonds_in_format
    if hb_angle_proxies.size() > 0:
      if work_params.format == "pymol":
        angles_in_format = hb_angle_proxies.as_pymol_dashes(
            pdb_hierarchy=pdb_hierarchy)
        print >> result_out, angles_in_format
  result = result_out.getvalue()
  out_prefix = os.path.basename(work_params.file_name[0])
  if work_params.output_prefix is not None:
    out_prefix = work_params.output_prefix
  filename = "%s_ss.eff" % out_prefix
  if work_params.format == "pymol":
    filename = "%s_ss.pml" % out_prefix
  outf = open(filename, "w")
  outf.write(result)
  outf.close()
  print >> out, result

  return os.path.abspath(filename)
############ end of secondary_structure_restraints_DN


def show_time(app, time_start, time_end):
  time_took = 0 # temporary of course
  if (round((time_end-time_start)/60, 1) < 1):
    time_took = str(app) + " finished in " + str(round((time_end-time_start), 1)) + " seconds (wallclock)."
  elif (round((time_end-time_start)/60/60, 1) < 1):
    time_took = str(app) + " finished in " + str(round((time_end-time_start)/60, 1)) + " minutes (wallclock)."
  else:
    time_took = str(app) + " finished in " + str(round((time_end-time_start)/60/60, 1)) + " hours (wallclock)."
  return time_took
############### end of show_time function


def write_custom_geometry(logfile, input_model_file_name, sigma_for_auto_geom, slack_for_auto_geom):

  ######## produce pymol format secondary structure restraints #########
  # I heard that running phenix commandline directly is not ideal from Nigel.
  # Therefore, I had used code directly rather than executing phenix executables at commandline such as calculating rmsd
  # However, I think that running phenix.secondary_structure_restraints is the best option here.
  # The reason is that I need to copy most of the codes in cctbx_project/mmtbx/command_line/secondary_structure_restraints.py
  #to use codes directly instead of running executables at commandline
  write_this = "Cryo_fit2 is generating pymol based secondary structure restraints for the user input model file to enforce a stronger sigma_for_auto_geom "
  if (sigma_for_auto_geom != None):
    write_this = write_this + "(e.g. " + str(sigma_for_auto_geom) + ") "
    
  if (slack_for_auto_geom != None):
    write_this = write_this + "and slack_for_auto_geom (e.g. " + str(slack_for_auto_geom) + ")"
    
  write_this = write_this + ".\n"
  print(write_this)
  logfile.write(write_this)
    
  if (sigma_for_auto_geom == None):
    write_this = "sigma_for_auto_geom = None till now, debug now\n"
    print(write_this)
    logfile.write(write_this)
    exit(1)
    
  make_pymol_ss_restraints = "phenix.secondary_structure_restraints " + input_model_file_name + " format=pymol"
  libtbx.easy_run.fully_buffered(make_pymol_ss_restraints)
  
  splited_input_model_file_name = input_model_file_name.split("/")
  input_model_file_name_wo_path = splited_input_model_file_name[len(splited_input_model_file_name)-1]
  ss_restraints_file_name = input_model_file_name_wo_path + "_ss.pml"
  
  if (path.isfile(ss_restraints_file_name) == False):
    write_this = '''phenix.secondary_structure_restraints can't run with a user input file.
    
To identify the cause of this error, run phenix.secondary_structure_restraints with a user input file.

If the error message is like
    "Sorry: number of groups of duplicate atom labels:  76
    total number of affected atoms:          152
    group "ATOM    .*.  CA  GLU F  55 .*.     C  "
          "ATOM    .*.  CA  GLU F  55 .*.     C  ""
then, provide input pdb file after solving duplicity issue.

    Most problems will be solved by running
        python <user phenix>/modules/cryo_fit2/util/solve_duplicate_atoms_by_adding_a_prefix/solve_duplicate_atoms_by_adding_a_prefix.py

    For multi-conformations, run
        phenix.pdbtools <user>.pdb remove_alt_confs=True
        (if MODEL #, ENDMDL are present, remove those lines before running phenix.pdbtools)


If the error message is like
    "Sorry: Multiple models not supported."
then provide input pdb file after leaving one model only.
    '''
    print(write_this)
    logfile.write(write_this)
    exit(1)
    
  ##### rewrite_pymol_ss_to_custom_geometry_ss
  eff_file_name = rewrite_pymol_ss_to_custom_geometry_ss(ss_restraints_file_name, sigma_for_auto_geom, slack_for_auto_geom)
  
  return eff_file_name
########### end of write_custom_geometry(input_model_file_name, sigma_for_auto_geom)


def write_geo(self, model_inp, file_name):
    geo_str = model_inp.restraints_as_geo(
        header="# Geometry restraints, cryo_fit2\n")
    
    with open(file_name, 'w') as f:
      f.write(geo_str)
    
    return True
###### end of write_geo
