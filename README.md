<this is under active development now, please consider to use cryo_fit1 instead>

@ .../phenix-..../modules

git clone git@github.com:cryoFIT/cryo_fit2.git

@ any folder

export PATH=:/Users/..../phenix-dev-xxxx/build/bin:$PATH

libtbx.configure cryo_fit2

enable phenix.cryo_fit2 to run


=========================
Development history
future    : Support multi-threading for search ideal parameters, calculate for 4 phenix cc
current   : map origin problem fixing (when map oriin is 0,0,0 -> no problem, when map origin is not 0,0,0 -> problem)
current   : Adding additional strong geometry restraints for nucleic acids
03/05/2019: Added RMSD calculation
01/23/2019: Uses phenix.real_space_refine style automatic map weight determination
11/01/2018: Refactored into Billy's template style (as phenix.EMRINGER)
08/04/2018: Pavel and Doonam made initial form of cryo_fit2 using phenix.dynamics