
run ./hydmol.py
run color_by_mutation.py

load HA_1819_NORTH_H1N1_QNT.pdb,qnt
load HA_1819_NORTH_H1N1_WHO.pdb,who

color gray

color_by_mutation who,qnt


set_view (\
     0.225687623,   -0.260932356,    0.938605309,\
    -0.972999334,   -0.012562572,    0.230465233,\
    -0.048343621,   -0.965275466,   -0.256722957,\
     0.000000000,    0.000000000, -339.883789062,\
   -11.484329224,  -25.482124329,    7.594523907,\
   268.851959229,  410.915618896,  -20.000000000 )
### cut above here and paste into script ###



show surface

set ray_opaque_background, 0

png surf_who_qnt18.png, dpi=100, ray=1

### cut below here and paste into script ###
set_view (\
     0.857454479,   -0.306516141,    0.413301498,\
    -0.446074635,   -0.042415179,    0.893988967,\
    -0.256491333,   -0.950919747,   -0.173099101,\
     0.000016525,   -0.000023562, -202.731155396,\
   -22.866436005,  -20.811021805,    7.496198654,\
   131.701934814,  273.765777588,  -20.000000000 )
### cut above here and paste into script ###

set ray_opaque_background, 0

png surf_who_qnt18.png, dpi=100, ray=1




quit