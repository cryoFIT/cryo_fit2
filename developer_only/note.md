note for speed: Although PHENIX group prefers to make cryo_fit2 works on laptop only as well, MD simulation itself takes a lot of computational time as in https://www.olcf.ornl.gov/2019/05/20/summit-charts-a-course-to-uncover-the-origins-of-genetic-diseases. 

Additionally, http://milou.science.uu.nl/services/HADDOCK2.2/haddock.php says that its docking may take 1~2 days. 

https://cluspro.bu.edu even needs 10 days with restraints.

Cryo_fit2 will be avilable on laptop only, but as of now, it may require 1~2 days on laptop if we model big biomolecules.


07/26/2019: Doo Nam confirmed that reoptimize_map_weight_after_each_cycle_during_final_MD_during_final_MD is effective to prevent nan error during core cryo-EM map based core dynamics run for full-tRNA.

06/25/2019:
Here, I confirmed that distance, angle, bonds moves from phenix.dynamics are independent from each temperature step (rather than cumulative)
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
