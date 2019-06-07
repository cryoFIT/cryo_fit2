import sys

def remove_sigma_in_custom_geom(args):
  custom_geom_eff = args[0]
  f_in = open(custom_geom_eff)
  out_file = custom_geom_eff[:-4] + '_no_sigma.eff'
  f_out = open(out_file, "w")
  
  for line in f_in:
    splited = line.split()
    if (splited[0] != "sigma"):
      f_out.write(line)
  f_in.close()
  f_out.close()

########## end of remove_sigma_in_custom_geom function


if (__name__ == "__main__"):
  args = sys.argv[1:]
  if len(args) == 0:
    print "Provide xxx_custom_geom.eff"
    print "Example: remove_sigma_in_custom_geom.py tRNA_initial_box_H.pdb_ss_custom_geom.eff"
    exit(1)
  remove_sigma_in_custom_geom(args)
  print "OK"