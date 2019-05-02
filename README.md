<this is under active development now, please consider to use cryo_fit1 instead>

@ .../phenix-..../modules

git clone git@github.com:cryoFIT/cryo_fit2.git

@ any folder

export PATH=:/Users/..../phenix-dev-xxxx/build/bin:$PATH

libtbx.configure cryo_fit2

enable phenix.cryo_fit2 to run


=========================

Development history

future    : Automatically support Nigel's phenix.eLBOW

future    : Automatically support Rob's multi-threading for search ideal parameters

future    : Automatically calculate all 4 phenix cc

future    : Support mtz input as well

future    : Automatically estimate compactness, then assign colder start_temparature for a wider (less compact) starting structure

current   : Automatically use sophiscated strong geometry restraints (Oleg's) for nucleic acids

current   : Automatically solve floppy protein local structures (Mg_channel, adenylate kinase)

05/02/2019: Automatically re-run cryo_fit2 until cc becomes a plateau

04/18/2019: Automatically update archaic nucleic acid names in user's pdb file

04/10/2019: Automatically deal origin issue of both "regular" (from emdb) maps and other maps (from phenix.map_box, Chimera's gaussian filter and relion image handler). Users don't need to pay any attention with respect to origin issue

03/05/2019: Automatically calculate RMSD before and after cryo_fit2

01/23/2019: Automatically optimize map weight following phenix.real_space_refine style

11/01/2018: Refactored into Billy's template style (as phenix.EMRINGER)

08/04/2018: Pavel and Doonam made initial form of cryo_fit2 using phenix.dynamics
