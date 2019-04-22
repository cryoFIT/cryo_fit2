<this is under active development now, please consider to use cryo_fit1 instead>

@ .../phenix-..../modules

git clone git@github.com:cryoFIT/cryo_fit2.git

@ any folder

export PATH=:/Users/..../phenix-dev-xxxx/build/bin:$PATH

libtbx.configure cryo_fit2

enable phenix.cryo_fit2 to run


=========================
Development history
future    : Support Nigel's phenix.eLBOW
future    : Support Rob's multi-threading for search ideal parameters
future    : Calculate all 4 phenix cc
future    : Support mtz input as well

current   : Automatically use sophiscated strong geometry restraints (Oleg's) for nucleic acids

04/18/2019: Update archaic nucleic acid names automatically
04/10/2019: Both "regular" (from emdb) maps and other maps (from phenix.map_box, Chimera's gaussian filter and relion image handler) have no problem of origin issue from user's perspective
03/05/2019: Added RMSD calculation
01/23/2019: Uses phenix.real_space_refine style automatic map weight determination
11/01/2018: Refactored into Billy's template style (as phenix.EMRINGER)
08/04/2018: Pavel and Doonam made initial form of cryo_fit2 using phenix.dynamics