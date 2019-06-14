note for speed: Although PHENIX group prefers to make cryo_fit2 works on laptop only as well, MD simulation itself takes a lot of computational time as in https://www.olcf.ornl.gov/2019/05/20/summit-charts-a-course-to-uncover-the-origins-of-genetic-diseases. 

Additionally, http://milou.science.uu.nl/services/HADDOCK2.2/haddock.php says that its docking may take 1~2 days. 

https://cluspro.bu.edu even needs 10 days with restraints.


Cryo_fit2 will be avilable on laptop only, but as of now, it may require 1~2 days on laptop if we model big biomolecules.


future    : Automatically support Rob's multi-threading for search ideal parameters

future    : Allow distance or pair restraints

future    : Automatically calculate all 4 phenix cc

future    : Support mtz input as well

future    : Automatically estimate compactness, then assign colder start_temparature for a wider (less compact) starting structure

future    : Automatically identify less fitted local region and fit that region only as real_space_refine2 does

future    : Automatically support Nigel's phenix.eLBOW
	    As Doonam sees cctbx_project/mmtbx/command_line/dynamics.py's "phenix.dynamics model.pdb ligands.cif", and https://www.phenix-online.org/documentation/reference/dynamics.html 's " any necessary restraints (CIF) files", he thought maybe eLBOW is already supported for cryo_fit2.
	    However, even when I provided eLBOW generated .cif ligand, cryo_fit2 shows "Number of atoms with unknown nonbonded energy type symbols:"

	    Currently, with unknown ligand (e.g. DDB), automatic map_weight optimization doesn't work.
	    With map_weight=x, cryo_fit2 runs but with "Number of atoms with unknown nonbonded energy type symbols:"



current    : Automatically use sophisticated strong geometry restraints (Oleg's) for nucleic acids and protein to solve floppy protein local structures (Mg_channel, adenylate kinase).
             tRNA is the smallest molecule (fastest one to run) but more challenging than L1 stalk.
	     /home/doonam/research/run/phenix/cryo_fit2/Mg_channel/20_steps_20_boost_strong_ss is promising for protein
	     
plan       : if combination works, use Rob's multi threading to automatically explore


here I confirmed that distance, angle, bonds moves from phenix.dynamics are independent from each temperature step (rather than cumulative)
...
  temp= 1000.0 dist_moved=  0.20 angles=  2.63 bonds= 0.021
  temp=  900.0 dist_moved=  0.20 angles=  2.46 bonds= 0.020
  temp=  800.0 dist_moved=  0.22 angles=  2.38 bonds= 0.019
  temp=  700.0 dist_moved=  0.23 angles=  2.31 bonds= 0.019
  temp=  600.0 dist_moved=  0.24 angles=  2.23 bonds= 0.018
  temp=  500.0 dist_moved=  0.24 angles=  2.14 bonds= 0.017
  temp=  400.0 dist_moved=  0.25 angles=  2.09 bonds= 0.016
  temp=  300.0 dist_moved=  0.25 angles=  2.03 bonds= 0.016
  temp=  200.0 dist_moved=  0.24 angles=  1.93 bonds= 0.015
  temp=  100.0 dist_moved=  0.22 angles=  1.85 bonds= 0.014
  temp=    0.0 dist_moved=  0.22 angles=  1.85 bonds= 0.014
...


05/17/2019: Updated tst2 regression to test short auto-rerun of MD

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
