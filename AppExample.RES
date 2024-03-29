BURBOUN flow analyses of Block with fractURes via BOUNdary and finite elements

BURBOUN application example, manuscript subm. to Comp. Geosci. January 2024

# blocks to be analyzed:     1
# elements per block side:  4
# elements per (inner) fracture:  4

BLOCK #:    1
COORDINATES: Xbl:  0.000    Ybl:  0.000     Xtr:  10.00    Ytr:  10.00    
CONDUCTIVITY: 0.1000E-07,# of fractures:  21,# of intersections:  10

EXTREMITIES COORDINATES, FRACTURE TRANSMISSIVITY  
   FRAC.#     X1          Y1          X2          Y2           T
        1  0.000       6.910       2.350       4.900      0.3034E-06
        2  0.000       3.720       2.350       4.900      0.3824E-05
        3  0.000       1.380       5.000       2.640      0.1099E-04
        4  2.350       4.900       2.650       4.640      0.1314E-06
        5  2.350       4.900       5.740       6.600      0.3214E-05
        6  2.650       4.640       5.000       2.640      0.2653E-06
        7  2.650       4.640       5.260       5.230      0.5103E-05
        8  5.000       2.640       6.590       1.280      0.3016E-06
        9  5.000       2.640       6.040       2.900      0.1275E-04
       10  6.590       1.280       8.090       0.000      0.3067E-06
       11  5.260       5.230       6.040       2.900      0.5143E-06
       12  5.260       5.230       6.620       5.540      0.4582E-05
       13  6.040       2.900       6.590       1.280      0.6376E-06
       14  6.040       2.900       8.290       3.470      0.1028E-04
       15  2.980       10.00       5.740       6.600      0.1798E-06
       16  5.740       6.600       6.620       5.540      0.1934E-06
       17  6.620       5.540       6.880       5.210      0.5684E-07
       18  6.620       5.540       10.00       6.300      0.4524E-05
       19  6.880       5.210       8.290       3.470      0.1973E-06
       20  6.880       5.210       10.00       5.110      0.2076E-05
       21  8.290       3.470       10.00       1.380      0.2158E-06
FRACTURE INCIDENCES PER INTERSECTION
INTE. #    F1   Ex.    F2   Ex.    F3   Ex.    F4   Ex.
      1     1     2     2     2     4     1     5     1
      2     4     2     6     1     7     1     0     0
      3     6     2     3     2     8     1     9     1
      4     8     2    10     1    13     2     0     0
      5     7     2    11     1    12     1     0     0
      6    11     2     9     2    13     1    14     1
      7    15     2     5     2    16     1     0     0
      8    16     2    12     2    17     1    18     1
      9    17     2    20     1    19     1     0     0
     10    14     2    19     2    21     1     0     0

-----
LOADING #1
BLOCK-EDGE DISCHARGES                            BUDGET
 -.2737E-07 0.1401E-06 0.3379E-07 -.1479E-06 -.1462E-08

NODAL VALUES
FLUXES AT THE BLOCK EDGES
BLOCK 1
    1 0.9785E-09    2 0.7334E-09    3 0.7597E-09    4 -.4771E-09
BLOCK 2
    5 0.3455E-08    6 -.7688E-09    7 -.5110E-09    8 0.6074E-11
BLOCK 3
    9 -.1169E-08   10 -.1331E-08   11 -.1506E-08   12 -.1187E-08
BLOCK 4
   13 0.1106E-10   14 0.5554E-09   15 -.1641E-08   16 0.4037E-10
FLUXES NORMAL TO THE FRACTURES
FRACTURE #  1
   17 -.2909E-08   18 0.8157E-09   19 0.6214E-10   20 0.2106E-10
FRACTURE #  2
   21 0.4069E-08   22 -.5476E-09   23 0.3499E-09   24 0.8485E-09
FRACTURE #  3
   25 -.2264E-10   26 0.3026E-09   27 0.3993E-09   28 0.6199E-09
FRACTURE #  4
   29 -.6070E-10   30 -.1065E-08   31 -.1758E-08   32 -.3873E-08
FRACTURE #  5
   33 0.2029E-08   34 0.1740E-08   35 0.1638E-08   36 0.1799E-08
FRACTURE #  6
   37 -.4879E-09   38 -.6062E-10   39 -.3084E-11   40 0.3337E-10
FRACTURE #  7
   41 -.2268E-08   42 -.1877E-08   43 -.1734E-08   44 -.1608E-08
FRACTURE #  8
   45 -.2087E-09   46 -.3683E-09   47 -.4943E-09   48 -.7621E-09
FRACTURE #  9
   49 0.9979E-09   50 0.9468E-09   51 0.9652E-09   52 0.9752E-09
FRACTURE # 10
   53 -.1167E-08   54 -.5065E-09   55 -.4381E-08   56 0.1084E-07
FRACTURE # 11
   57 0.1429E-09   58 0.8790E-10   59 0.4025E-10   60 -.8010E-11
FRACTURE # 12
   61 -.1564E-08   62 -.1662E-08   63 -.1697E-08   64 -.1761E-08
FRACTURE # 13
   65 -.5534E-10   66 -.1098E-09   67 -.1818E-09   68 -.4732E-09
FRACTURE # 14
   69 0.9946E-09   70 0.1014E-08   71 0.1051E-08   72 0.1271E-08
FRACTURE # 15
   73 0.2758E-09   74 0.3737E-09   75 0.5628E-09   76 0.1171E-08
FRACTURE # 16
   77 0.2448E-08   78 0.1257E-08   79 0.8686E-09   80 0.3923E-09
FRACTURE # 17
   81 -.3701E-09   82 -.6229E-10   83 0.4849E-10   84 0.1371E-09
FRACTURE # 18
   85 -.9777E-09   86 -.3977E-09   87 -.3351E-09   88 0.9801E-09
FRACTURE # 19
   89 0.4213E-09   90 0.5120E-09   91 0.6620E-09   92 0.1207E-08
FRACTURE # 20
   93 -.4103E-09   94 -.5091E-09   95 -.4235E-09   96 0.3018E-09
FRACTURE # 21
   97 0.1601E-08   98 0.9577E-09   99 0.3155E-08  100 -.1436E-07
HYDRAULIC HEADS AT THE FRACTURES
FRACTURE #  1
  101 0.3431      102 0.4080      103 0.4734      104 0.5390    
FRACTURE #  2
  105 0.6208      106 0.6068      107 0.5928      108 0.5789    
FRACTURE #  3
  109 0.8565      110 0.8455      111 0.8345      112 0.8236    
FRACTURE #  4
  113 0.5662      114 0.5548      115 0.5433      116 0.5316    
FRACTURE #  5
  117 0.5661      118 0.5549      119 0.5442      120 0.5339    
FRACTURE #  6
  121 0.5627      122 0.6358      123 0.7087      124 0.7817    
FRACTURE #  7
  125 0.5230      126 0.5174      127 0.5116      128 0.5057    
FRACTURE #  8
  129 0.8236      130 0.8343      131 0.8446      132 0.8545    
FRACTURE #  9
  133 0.8174      134 0.8158      135 0.8143      136 0.8128    
FRACTURE # 10
  137 0.8774      138 0.9130      139 0.9478      140 0.9809    
FRACTURE # 11
  141 0.5413      142 0.6186      143 0.6959      144 0.7733    
FRACTURE # 12
  145 0.4985      146 0.4900      147 0.4815      148 0.4729    
FRACTURE # 13
  149 0.8179      150 0.8298      151 0.8416      152 0.8534    
FRACTURE # 14
  153 0.8113      154 0.8099      155 0.8085      156 0.8072    
FRACTURE # 15
  157 0.6358E-01  158 0.1922      159 0.3234      160 0.4587    
FRACTURE # 16
  161 0.5206      162 0.5049      163 0.4900      164 0.4757    
FRACTURE # 17
  165 0.4747      166 0.4869      167 0.4991      168 0.5113    
FRACTURE # 18
  169 0.4563      170 0.4317      171 0.4070      172 0.3823    
FRACTURE # 19
  173 0.5527      174 0.6240      175 0.6961      176 0.7693    
FRACTURE # 20
  177 0.5139      178 0.5070      179 0.4998      180 0.4926    
FRACTURE # 21
  181 0.8124      182 0.8266      183 0.8435      184 0.8619    

-----
LOADING #2
BLOCK-EDGE DISCHARGES                            BUDGET
 -.8126E-07 -.5155E-06 -.4278E-08 0.5990E-06 -.2012E-08

NODAL VALUES
FLUXES AT THE BLOCK EDGES
BLOCK 1
    1 0.3782E-09    2 0.1064E-08    3 0.1959E-08    4 -.4276E-09
BLOCK 2
    5 0.7015E-08    6 0.3017E-08    7 -.3793E-09    8 0.1226E-08
BLOCK 3
    9 0.1140E-10   10 0.4897E-09   11 0.1134E-08   12 0.2305E-09
BLOCK 4
   13 -.9595E-09   14 -.6563E-09   15 -.6054E-09   16 -.6046E-09
FLUXES NORMAL TO THE FRACTURES
FRACTURE #  1
   17 0.4621E-10   18 -.2103E-10   19 0.8132E-10   20 -.5653E-10
FRACTURE #  2
   21 -.7667E-10   22 -.1559E-09   23 -.6341E-09   24 -.3461E-08
FRACTURE #  3
   25 -.4858E-09   26 -.1229E-08   27 -.2248E-08   28 -.2814E-08
FRACTURE #  4
   29 -.8200E-09   30 0.5376E-08   31 0.9294E-08   32 0.2138E-07
FRACTURE #  5
   33 -.8744E-08   34 -.5952E-08   35 -.4677E-08   36 -.4330E-08
FRACTURE #  6
   37 0.3056E-08   38 0.5497E-09   39 0.1003E-09   40 0.8257E-11
FRACTURE #  7
   41 0.1115E-07   42 0.8267E-08   43 0.7059E-08   44 0.6209E-08
FRACTURE #  8
   45 -.1573E-09   46 -.4489E-09   47 -.7673E-09   48 -.1545E-08
FRACTURE #  9
   49 -.2508E-08   50 -.2652E-08   51 -.2684E-08   52 -.2598E-08
FRACTURE # 10
   53 -.2707E-08   54 -.1114E-08   55 -.3944E-08   56 0.8588E-08
FRACTURE # 11
   57 -.4383E-09   58 -.6107E-09   59 -.4673E-09   60 -.1035E-09
FRACTURE # 12
   61 0.5707E-08   62 0.5374E-08   63 0.4780E-08   64 0.3180E-08
FRACTURE # 13
   65 -.8510E-12   66 -.3877E-09   67 -.7398E-09   68 -.1668E-08
FRACTURE # 14
   69 -.3124E-08   70 -.4037E-08   71 -.4874E-08   72 -.6520E-08
FRACTURE # 15
   73 -.1360E-08   74 -.4179E-09   75 -.1032E-08   76 -.2280E-08
FRACTURE # 16
   77 -.5106E-08   78 -.2766E-08   79 -.1983E-08   80 -.1173E-08
FRACTURE # 17
   81 0.1267E-08   82 0.3051E-08   83 0.4732E-08   84 0.1070E-07
FRACTURE # 18
   85 -.5922E-09   86 -.3209E-09   87 -.5224E-09   88 0.2748E-08
FRACTURE # 19
   89 0.1184E-08   90 -.8675E-09   91 -.1901E-08   92 -.4946E-08
FRACTURE # 20
   93 0.6267E-08   94 0.4631E-08   95 0.2881E-08   96 0.2176E-09
FRACTURE # 21
   97 -.6584E-08   98 -.3489E-08   99 0.9525E-09  100 -.1699E-07
HYDRAULIC HEADS AT THE FRACTURES
FRACTURE #  1
  101 0.1871E-01  102 0.5619E-01  103 0.9367E-01  104 0.1312    
FRACTURE #  2
  105 0.1880E-01  106 0.5638E-01  107 0.9395E-01  108 0.1314    
FRACTURE #  3
  109 0.2102E-01  110 0.6299E-01  111 0.1048      112 0.1462    
FRACTURE #  4
  113 0.1999      114 0.2998      115 0.4000      116 0.5010    
FRACTURE #  5
  117 0.1642      118 0.1909      119 0.2158      120 0.2394    
FRACTURE #  6
  121 0.5013      122 0.4043      123 0.3091      124 0.2142    
FRACTURE #  7
  125 0.5634      126 0.5868      127 0.6109      128 0.6357    
FRACTURE #  8
  129 0.1881      130 0.2307      131 0.2728      132 0.3141    
FRACTURE #  9
  133 0.1699      134 0.1760      135 0.1821      136 0.1882    
FRACTURE # 10
  137 0.3948      138 0.5142      139 0.6323      140 0.7488    
FRACTURE # 11
  141 0.5914      142 0.4776      143 0.3633      144 0.2486    
FRACTURE # 12
  145 0.6592      146 0.6813      147 0.7035      148 0.7259    
FRACTURE # 13
  149 0.2093      150 0.2452      151 0.2811      152 0.3167    
FRACTURE # 14
  153 0.1952      154 0.2030      155 0.2107      156 0.2182    
FRACTURE # 15
  157 0.2980      158 0.2921      159 0.2821      160 0.2647    
FRACTURE # 16
  161 0.3134      162 0.4363      163 0.5575      164 0.6774    
FRACTURE # 17
  165 0.7535      166 0.7867      167 0.8205      168 0.8553    
FRACTURE # 18
  169 0.7700      170 0.8357      171 0.9014      172 0.9670    
FRACTURE # 19
  173 0.7928      174 0.6328      175 0.4715      176 0.3068    
FRACTURE # 20
  177 0.8881      178 0.9187      179 0.9507      180 0.9835    
FRACTURE # 21
  181 0.3276      182 0.5296      183 0.7245      184 0.9156    

-----
LOADING #3
BLOCK-EDGE DISCHARGES                            BUDGET
 0.3108E-07 -.1459E-06 -.3385E-07 0.1488E-06 0.1908E-09

NODAL VALUES
FLUXES AT THE BLOCK EDGES
BLOCK 1
    1 -.9980E-09    2 -.7488E-09    3 -.9622E-09    4 -.1059E-08
BLOCK 2
    5 0.7359E-09    6 0.6989E-09    7 -.5473E-09    8 0.8117E-10
BLOCK 3
    9 0.1162E-08   10 0.1291E-08   11 0.1975E-08   12 0.1264E-08
BLOCK 4
   13 0.1201E-09   14 0.1103E-08   15 -.7781E-09   16 0.1842E-09
FLUXES NORMAL TO THE FRACTURES
FRACTURE #  1
   17 -.3199E-08   18 0.7444E-09   19 -.5515E-10   20 0.1867E-10
FRACTURE #  2
   21 0.2733E-08   22 -.4135E-09   23 -.2549E-09   24 -.8526E-09
FRACTURE #  3
   25 0.1442E-09   26 -.2145E-09   27 -.3899E-09   28 -.6156E-09
FRACTURE #  4
   29 0.2253E-10   30 0.1064E-08   31 0.1761E-08   32 0.3880E-08
FRACTURE #  5
   33 -.2007E-08   34 -.1710E-08   35 -.1621E-08   36 -.1797E-08
FRACTURE #  6
   37 0.4959E-09   38 0.6681E-10   39 0.5461E-11   40 -.3335E-10
FRACTURE #  7
   41 0.2269E-08   42 0.1878E-08   43 0.1734E-08   44 0.1608E-08
FRACTURE #  8
   45 0.2111E-09   46 0.3728E-09   47 0.4988E-09   48 0.7764E-09
FRACTURE #  9
   49 -.9982E-09   50 -.9460E-09   51 -.9645E-09   52 -.9740E-09
FRACTURE # 10
   53 0.1036E-08   54 0.5372E-09   55 0.3886E-09   56 0.5718E-10
FRACTURE # 11
   57 -.1427E-09   58 -.8782E-10   59 -.4023E-10   60 0.8026E-11
FRACTURE # 12
   61 0.1564E-08   62 0.1661E-08   63 0.1695E-08   64 0.1760E-08
FRACTURE # 13
   65 0.5388E-10   66 0.1034E-09   67 0.1670E-09   68 0.4574E-09
FRACTURE # 14
   69 -.9962E-09   70 -.1019E-08   71 -.1056E-08   72 -.1285E-08
FRACTURE # 15
   73 -.1792E-08   74 0.1634E-11   75 -.5556E-09   76 -.1166E-08
FRACTURE # 16
   77 -.2455E-08   78 -.1259E-08   79 -.8702E-09   80 -.3931E-09
FRACTURE # 17
   81 0.3724E-09   82 0.6286E-10   83 -.4786E-10   84 -.1359E-09
FRACTURE # 18
   85 0.9735E-09   86 0.4109E-09   87 -.1132E-09   88 0.1718E-08
FRACTURE # 19
   89 -.4201E-09   90 -.5101E-09   91 -.6577E-09   92 -.1210E-08
FRACTURE # 20
   93 0.4091E-09   94 0.5111E-09   95 0.4297E-09   96 -.1868E-09
FRACTURE # 21
   97 -.1597E-08   98 -.1117E-08   99 -.1295E-09  100 -.1440E-08
HYDRAULIC HEADS AT THE FRACTURES
FRACTURE #  1
  101 0.6596      102 0.5932      103 0.5270      104 0.4610    
FRACTURE #  2
  105 0.3789      106 0.3930      107 0.4070      108 0.4210    
FRACTURE #  3
  109 0.1435      110 0.1546      111 0.1657      112 0.1766    
FRACTURE #  4
  113 0.4337      114 0.4451      115 0.4566      116 0.4683    
FRACTURE #  5
  117 0.4338      118 0.4450      119 0.4557      120 0.4660    
FRACTURE #  6
  121 0.4373      122 0.3642      123 0.2914      124 0.2185    
FRACTURE #  7
  125 0.4770      126 0.4826      127 0.4883      128 0.4943    
FRACTURE #  8
  129 0.1767      130 0.1660      131 0.1557      132 0.1458    
FRACTURE #  9
  133 0.1829      134 0.1844      135 0.1860      136 0.1875    
FRACTURE # 10
  137 0.1230      138 0.8743E-01  139 0.5227E-01  140 0.1740E-01
FRACTURE # 11
  141 0.4587      142 0.3815      143 0.3042      144 0.2269    
FRACTURE # 12
  145 0.5015      146 0.5100      147 0.5185      148 0.5271    
FRACTURE # 13
  149 0.1824      150 0.1705      151 0.1587      152 0.1470    
FRACTURE # 14
  153 0.1890      154 0.1905      155 0.1919      156 0.1932    
FRACTURE # 15
  157 0.9388      158 0.8089      159 0.6770      160 0.5414    
FRACTURE # 16
  161 0.4793      162 0.4950      163 0.5099      164 0.5243    
FRACTURE # 17
  165 0.5252      166 0.5131      167 0.5009      168 0.4887    
FRACTURE # 18
  169 0.5436      170 0.5682      171 0.5929      172 0.6176    
FRACTURE # 19
  173 0.4473      174 0.3762      175 0.3042      176 0.2311    
FRACTURE # 20
  177 0.4861      178 0.4931      179 0.5002      180 0.5074    
FRACTURE # 21
  181 0.1889      182 0.1766      183 0.1621      184 0.1467    

-----
LOADING #4
BLOCK-EDGE DISCHARGES                            BUDGET
 0.8497E-07 0.5096E-06 0.4214E-08 -.5981E-06 0.7385E-09

NODAL VALUES
FLUXES AT THE BLOCK EDGES
BLOCK 1
    1 -.3977E-09    2 -.1080E-08    3 -.2161E-08    4 -.1108E-08
BLOCK 2
    5 -.2824E-08    6 -.3087E-08    7 -.6791E-09    8 -.1139E-08
BLOCK 3
    9 -.1756E-10   10 -.5293E-09   11 -.6647E-09   12 -.1536E-09
BLOCK 4
   13 0.1091E-08   14 0.2314E-08   15 -.1814E-08   16 0.8292E-09
FLUXES NORMAL TO THE FRACTURES
FRACTURE #  1
   17 -.6154E-08   18 0.1581E-08   19 -.7434E-10   20 0.9628E-10
FRACTURE #  2
   21 0.6879E-08   22 -.8053E-09   23 0.7291E-09   24 0.3457E-08
FRACTURE #  3
   25 0.6073E-09   26 0.1317E-08   27 0.2257E-08   28 0.2818E-08
FRACTURE #  4
   29 0.7819E-09   30 -.5378E-08   31 -.9290E-08   32 -.2137E-07
FRACTURE #  5
   33 0.8766E-08   34 0.5981E-08   35 0.4695E-08   36 0.4331E-08
FRACTURE #  6
   37 -.3048E-08   38 -.5435E-09   39 -.9790E-10   40 -.8246E-11
FRACTURE #  7
   41 -.1115E-07   42 -.8266E-08   43 -.7059E-08   44 -.6209E-08
FRACTURE #  8
   45 0.1597E-09   46 0.4535E-09   47 0.7719E-09   48 0.1559E-08
FRACTURE #  9
   49 0.2507E-08   50 0.2653E-08   51 0.2685E-08   52 0.2599E-08
FRACTURE # 10
   53 0.2577E-08   54 0.1144E-08   55 -.4818E-10   56 0.2310E-08
FRACTURE # 11
   57 0.4385E-09   58 0.6108E-09   59 0.4673E-09   60 0.1036E-09
FRACTURE # 12
   61 -.5708E-08   62 -.5375E-08   63 -.4781E-08   64 -.3181E-08
FRACTURE # 13
   65 -.6600E-12   66 0.3813E-09   67 0.7249E-09   68 0.1653E-08
FRACTURE # 14
   69 0.3122E-08   70 0.4031E-08   71 0.4869E-08   72 0.6507E-08
FRACTURE # 15
   73 -.1562E-09   74 0.7933E-09   75 0.1040E-08   76 0.2285E-08
FRACTURE # 16
   77 0.5099E-08   78 0.2764E-08   79 0.1982E-08   80 0.1172E-08
FRACTURE # 17
   81 -.1265E-08   82 -.3050E-08   83 -.4731E-08   84 -.1070E-07
FRACTURE # 18
   85 0.5879E-09   86 0.3341E-09   87 0.7402E-10   88 -.4963E-10
FRACTURE # 19
   89 -.1183E-08   90 0.8695E-09   91 0.1905E-08   92 0.4943E-08
FRACTURE # 20
   93 -.6268E-08   94 -.4629E-08   95 -.2875E-08   96 -.1026E-09
FRACTURE # 21
   97 0.6588E-08   98 0.3329E-08   99 0.2073E-08  100 0.1198E-08
HYDRAULIC HEADS AT THE FRACTURES
FRACTURE #  1
  101 0.9841      102 0.9450      103 0.9067      104 0.8688    
FRACTURE #  2
  105 0.9810      106 0.9434      107 0.9059      108 0.8684    
FRACTURE #  3
  109 0.9790      110 0.9371      111 0.8954      112 0.8540    
FRACTURE #  4
  113 0.7999      114 0.7001      115 0.5999      116 0.4989    
FRACTURE #  5
  117 0.8356      118 0.8090      119 0.7840      120 0.7604    
FRACTURE #  6
  121 0.4987      122 0.5957      123 0.6910      124 0.7860    
FRACTURE #  7
  125 0.4366      126 0.4132      127 0.3890      128 0.3643    
FRACTURE #  8
  129 0.8121      130 0.7696      131 0.7275      132 0.6862    
FRACTURE #  9
  133 0.8304      134 0.8243      135 0.8181      136 0.8121    
FRACTURE # 10
  137 0.6056      138 0.4862      139 0.3677      140 0.2495    
FRACTURE # 11
  141 0.4086      142 0.5225      143 0.6369      144 0.7516    
FRACTURE # 12
  145 0.3407      146 0.3186      147 0.2964      148 0.2741    
FRACTURE # 13
  149 0.7910      150 0.7551      151 0.7193      152 0.6837    
FRACTURE # 14
  153 0.8051      154 0.7973      155 0.7897      156 0.7822    
FRACTURE # 15
  157 0.7043      158 0.7089      159 0.7183      160 0.7354    
FRACTURE # 16
  161 0.6865      162 0.5636      163 0.4425      164 0.3225    
FRACTURE # 17
  165 0.2464      166 0.2132      167 0.1795      168 0.1447    
FRACTURE # 18
  169 0.2300      170 0.1642      171 0.9852E-01  172 0.3284E-01
FRACTURE # 19
  173 0.2072      174 0.3674      175 0.5287      176 0.6935    
FRACTURE # 20
  177 0.1120      178 0.8130E-01  179 0.4929E-01  180 0.1647E-01
FRACTURE # 21
  181 0.6737      182 0.4736      183 0.2810      184 0.9297E-01
