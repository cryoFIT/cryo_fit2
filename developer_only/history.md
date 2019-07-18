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

07/13/2019: Now eff file works correctly when provided as an argument (If the file was generated and then appended during cryo_fit2, sigma not works)

06/27/2019: Cryo_fit2 uses available number of cores automatically using from libtbx.introspection import number_of_processors
	     
06/26/2019: Cryo_fit2 uses Rob's multi threading to automatically explore different MD parameters.

06/26/2019: As I checked with several extreme cases, different parameters combination exploration does have decent reproducibility.

05/17/2019: Updated tst2 regression to test short auto-rerun of MD

05/13/2019: Automatically reoptimize map weight for every 5k MD iteration

05/10/2019: Fine-tuneed default parameters (decrease number_of_steps, increase cool_rate) to better maintain starting geometry. Seems to work for L1_stalk

05/04/2019: Added a new option: total_steps

05/02/2019: Automatically re-run cryo_fit2 until cc becomes a plateau

04/18/2019: Automatically update archaic nucleic acid names in user's pdb file

04/10/2019: Automatically deal origin issue of both "regular" (from emdb) maps and other maps (from phenix.map_box, Chimera's gaussian filter and relion image handler). Users don't need to pay any attention with respect to origin issue

03/05/2019: Automatically calculate RMSD after cryo_fit2

01/23/2019: Automatically optimize map weight following phenix.real_space_refine style

11/01/2018: Refactored into Billy's template style (as phenix.EMRINGER)

08/04/2018: Pavel and Doonam made initial form of cryo_fit2 using phenix.dynamics
