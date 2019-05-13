<this is under active development now, please consider to use cryo_fit1 instead>

@ .../phenix-..../modules

git clone git@github.com:cryoFIT/cryo_fit2.git

@ any folder

export PATH=:/Users/..../phenix-dev-xxxx/build/bin:$PATH

libtbx.configure cryo_fit2

enable phenix.cryo_fit2 to run


=========================

Development history

future    : Automatically support Rob's multi-threading for search ideal parameters

future    : Allow distance or pair restraints

future    : Automatically calculate all 4 phenix cc

future    : Support mtz input as well

future    : Automatically estimate compactness, then assign colder start_temparature for a wider (less compact) starting structure

future    : Automatically identify less fitted local region and fit that region only as real_space_refine2 does

future   : Automatically use sophisticated strong geometry restraints (Oleg's) for nucleic acids and protein to solve floppy protein local structures (Mg_channel, adenylate kinase). I may need to use cryo_fit2_command.py to automatically use .eff file

future    : Automatically support Nigel's phenix.eLBOW (as Doonam sees cctbx_project/mmtbx/command_line/dynamics.py's "phenix.dynamics model.pdb ligands.cif", and https://www.phenix-online.org/documentation/reference/dynamics.html 's " 
ny necessary restraints (CIF) files", maybe eLBOW is already supported?). Find DDB from Venki's ribosome then confirm that cryo_fit2 can work with .cif

current   : update tst2 regression to test short auto-rerun of md

05/13/2019: Automatically reoptimize map weight for every 5k MD iteration

05/10/2019: Fine-tuneed default parameters (decrease number_of_steps, increase cool_rate) to better maintain starting geometry. Seems to work for L1_stalk

05/04/2019: Added a new option: total_number_of_steps

05/02/2019: Automatically re-run cryo_fit2 until cc becomes a plateau

04/18/2019: Automatically update archaic nucleic acid names in user's pdb file

04/10/2019: Automatically deal origin issue of both "regular" (from emdb) maps and other maps (from phenix.map_box, Chimera's gaussian filter and relion image handler). Users don't need to pay any attention with respect to origin issue

03/05/2019: Automatically calculate RMSD after cryo_fit2

01/23/2019: Automatically optimize map weight following phenix.real_space_refine style

11/01/2018: Refactored into Billy's template style (as phenix.EMRINGER)

08/04/2018: Pavel and Doonam made initial form of cryo_fit2 using phenix.dynamics
