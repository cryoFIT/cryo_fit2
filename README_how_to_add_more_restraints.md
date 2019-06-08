[Rationale]
By default, cryo_fit2 already enforces Oleg's secondary structure restraints for protein modeling (e.g. secondary_structure.enabled=True secondary_structure.protein.remove_outliers=True).

By default, cryo_fit2 already enforces Oleg's secondary structure restraints for nucleic acid modeling (e.g. secondary_structure.nucleic_acid.enabled=True).

However, these restraints are applied loosely.

Therefore, if preferred, cryo_fit2 users can enforce stronger weight to these secondary structure restraints (e.g. sigma = 0.021)

Then, cryo_fit2 can better preserve starting secondary structure of protein and base-pairing of nucleic acids.
Then, cryo_fit2 can better fit to cryo-EM map faster.

This unsually strong sigma (e.g. 0.021) may distort local geometry little awkwardly.
Therefore, at coolor temperature, cryo_fit2 refines structure with a original sigma (e.g. 1) again before spitting a fitted structure to a user.



[How to add additional restraints]
cryo_fit2 <user>.map <user>.pdb resolution=<x> strong_ss=True



[How to see phenix generated secondary restraints]
Open <user>.pdb and <user>pdb_ss.pml at pymol
For example, at pymol commandline, @initial_3jcf_fitted_to_final_6553_no_Mg.pdb_ss.pml
