# LIBTBX_SET_DISPATCHER_NAME cryo_fit2.run_tests
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

import glob, iotbx.pdb.hierarchy, os, subprocess, sys, time
from libtbx import phil
import libtbx.phil.command_line
from libtbx.utils import Sorry
from subprocess import check_output
import libtbx.load_env
import shutil
import time

########## <begin> import util py files
cryo_fit2_repository_dir = libtbx.env.dist_path("cryo_fit2") # Locate phenix.cryo_fit.run_tests executable
util_path = cryo_fit2_repository_dir + "/util/"
print (util_path)
sys.path.insert(0, util_path)
from util import *
########## <end> import util py files


def test_fn ():
    time_total_start = time.time()
    
    assert len(os.listdir(os.getcwd()))==0, 'run in an empty directory' # added by Nigel so that this test runs in a clear path

#    print "This phenix.cryo_fit2.run_tests executable comes from ", cryo_fit2_repository_dir
    splited = cryo_fit2_repository_dir.split("/")
    regression_path = ''
    for i in range(len(splited)-1):
      regression_path = regression_path + splited[i] + "/"

    ############# test 1, simplest biomolecule with total_steps ###############
    regression_path = os.path.join(cryo_fit2_repository_dir,
                                     'regression')
    print "regression_path:", regression_path
    os.chdir(regression_path)

    time_start = time.time()
    command_string = "python tst1_cryo_fit2_peptide_total_steps.py" % locals()
    print "command_string:", command_string
    rc1 = libtbx.easy_run.call(command=command_string)
    assert rc1==0 # make sure there is no error with this test
    time_end = time.time()
    print "Minutes took for this tst test:", ( round(((time_end-time_start)/60),1) )
    
    # remove no longer needed folder and input_command file
    rm_command_string = "rm -r cryo_fit2.input_command.txt output_*"
    libtbx.easy_run.fully_buffered(rm_command_string)



    ######### test 2, simplest biomolecule with regular repeating of small MD running ###############
    regression_path = os.path.join(cryo_fit2_repository_dir,
                                     'regression')
    print "regression_path:", regression_path
    os.chdir(regression_path)

    time_start = time.time()
    command_string = "python tst2_cryo_fit2_peptide_auto-rerun.py" % locals()
    print "command_string:", command_string
    rc2 = libtbx.easy_run.call(command=command_string)
    assert rc2==0 # make sure there is no error with this test
    time_end = time.time()
    print "Minutes took for this tst test:", ( round(((time_end-time_start)/60),1) )
    
    # remove no longer needed folder and input_command file
    rm_command_string = "rm -r cryo_fit2.input_command.txt output_*"
    libtbx.easy_run.fully_buffered(rm_command_string)


    ############# test 3, simplest biomolecule test stong_ss ###############
    regression_path = os.path.join(cryo_fit2_repository_dir,
                                     'regression')    
    print "regression_path:", regression_path
    os.chdir(regression_path)

    time_start = time.time()
    command_string = "python tst3_cryo_fit2_peptide_make_ss_for_stronger_ss.py" % locals()
    print "command_string:", command_string
    rc3 = libtbx.easy_run.call(command=command_string)
    assert rc3==0 # make sure there is no error with this test
    time_end = time.time()
    print "Minutes took for this tst test:", ( round(((time_end-time_start)/60),1) )
    
    # remove no longer needed folder and input_command file
    rm_command_string = "rm -r cryo_fit2.input_command.txt output_*"
    libtbx.easy_run.fully_buffered(rm_command_string)
    
    
    ############# test 4
    regression_path = os.path.join(cryo_fit2_repository_dir,
                                     'regression')
    print "regression_path:", regression_path
    os.chdir(regression_path)
    
    time_start = time.time()
    command_string = "python tst4_cryo_fit2_RNA_parameters_exploration.py" % locals()
    print "command_string:", command_string
    rc4 = libtbx.easy_run.call(command=command_string)
    assert rc4==0 # make sure there is no error with this test
    time_end = time.time()
    print "Minutes took for this tst test:", ( round(((time_end-time_start)/60),1) )
    
    # remove no longer needed folder and input_command file
    rm_command_string = "rm -r cryo_fit2.input_command.txt output_* parameters_exploration"
    libtbx.easy_run.fully_buffered(rm_command_string)
    
    
    ############# test 5
    regression_path = os.path.join(cryo_fit2_repository_dir,
                                     'regression')
    print "regression_path:", regression_path
    os.chdir(regression_path)

    time_start = time.time()
    command_string = "python tst5_cryo_fit2_peptide_parameters_exploration.py" % locals()
    print "command_string:", command_string
    rc5 = libtbx.easy_run.call(command=command_string)
    assert rc5==0 # make sure there is no error with this test
    
    time_end = time.time()
    print "Minutes took for this tst test:", ( round(((time_end-time_start)/60),1) )
    
    # remove a no longer needed folder and an input command file
    rm_command_string = "rm -r cryo_fit2.input_command.txt output_* parameters_exploration"
    libtbx.easy_run.fully_buffered(rm_command_string)
    
    
    
    ############# test 6
    regression_path = os.path.join(cryo_fit2_repository_dir,
                                     'regression')
    print "regression_path:", regression_path
    os.chdir(regression_path)

    time_start = time.time()
    command_string = "python tst6_cryo_fit2_protein_RNA_full.py" % locals()
    print "command_string:", command_string
    rc5 = libtbx.easy_run.call(command=command_string)
    assert rc5==0 # make sure there is no error with this test
    
    time_end = time.time()
    print "Minutes took for this tst test:", ( round(((time_end-time_start)/60),1) )
    
    # remove a no longer needed folder and an input command file
    rm_command_string = "rm -r cryo_fit2.input_command.txt output_* parameters_exploration"
    libtbx.easy_run.fully_buffered(rm_command_string)
    
    
    time_total_end = time.time()
    time_took = show_time("All regression tests", time_total_start, time_total_end)
    print (time_took)
    
    if ((rc1 == 0) and (rc2 == 0) and (rc3 == 0) and (rc4 == 0) and (rc5 == 0)):
      return 0 # success
    else:
      return 1 # fail
###### end of test_fn

  
if (__name__ == "__main__") :
    test_fn()
    
