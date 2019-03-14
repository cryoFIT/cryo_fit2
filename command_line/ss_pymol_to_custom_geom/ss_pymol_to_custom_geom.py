import sys

def rewrite(args):
  user_input_pymol_ss = args[0]
  f_in = open(user_input_pymol_ss)
  out_file = user_input_pymol_ss[:-4] + '_custom_geom.eff'
  f_out = open(out_file, "w")
  f_out.write('''geometry_restraints {
  edits {
''')
  for line in f_in:
    dist_candidate = line[0:4]
    if (dist_candidate == "dist"):
      print line
      splited = line.split()
      
      write_this = "    bond {\n"
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
      else: # default
        f_out.write("      distance_ideal = 2.91\n")
      # reference modules/cctbx_project/mmtbx/secondary_structure/nucleic_acids.py
      
      f_out.write("      sigma = 0.021\n")
      
      write_this = "    }\n"
      f_out.write(write_this)
  f_out.write('''  }
}
''')
  f_in.close()
  f_out.close()

########## end of rewrite function

#dist chain "A" and resi   42  and name  N4  and alt '', chain "A" and resi   28  and name  O6  and alt ''
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

if (__name__ == "__main__"):
  args = sys.argv[1:]
  if len(args) == 0:
    print "Provide *.pdb_ss.pml file as a starting ss pymol script file"
    print "Example: ss_pymol_to_custom_geom.py tRNA_initial_box.pdb_ss.pml"
    exit(1)
  rewrite(args)
  print "OK"