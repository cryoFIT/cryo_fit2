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

cryo_fit2_repository_dir = libtbx.env.dist_path("cryo_fit2") # Locate phenix.cryo_fit.run_tests executable

if (__name__ == "__main__") :

    assert len(os.listdir(os.getcwd()))==0, 'run in an empty directory' # added by Nigel so that this test runs in a clear path

#    print "This phenix.cryo_fit2.run_tests executable comes from ", cryo_fit2_repository_dir

    splited = cryo_fit2_repository_dir.split("/")
    regression_path = ''
    for i in range(len(splited)-1):
      regression_path = regression_path + splited[i] + "/"


    ############# test 1, simplest biomolecule with total_number_of_steps ###############
    regression_path = os.path.join(cryo_fit2_repository_dir,
                                     'regression')
    print "regression_path:", regression_path
    os.chdir(regression_path)

    command_string = "python tst1_cryo_fit2_test_total_number_of_steps.py" % locals()
    print "command_string:", command_string
    rc = libtbx.easy_run.call(command=command_string)
    assert rc==0 # make sure there is no error with this test
    
    # remove no longer needed folder and input_command file
    rm_command_string = "rm -r cryo_fit2.input_command.txt output_*"
    libtbx.easy_run.fully_buffered(rm_command_string)



    ######### test 2, simplest biomolecule with regular repeating of small MD running ###############
    regression_path = os.path.join(cryo_fit2_repository_dir,
                                     'regression')
    print "regression_path:", regression_path
    os.chdir(regression_path)

    command_string = "python tst2_cryo_fit2_test_auto-rerun.py" % locals()
    print "command_string:", command_string
    rc = libtbx.easy_run.call(command=command_string)
    assert rc==0 # make sure there is no error with this test
    
    # remove no longer needed folder and input_command file
    rm_command_string = "rm -r cryo_fit2.input_command.txt output_*"
    libtbx.easy_run.fully_buffered(rm_command_string)




