# SCCS $Id: survexp.usr.s,v 5.6 1999-01-29 15:52:22 therneau Exp $
#
# Create the US total hazards table, by race
#   The raw numbers below are q* 10^5.  Note that there are 24 leap years/100
#
temp <- c(1538,1040,556,4812,487,265,190,153,138,124,114,106,102,100,101,103,
    114,127,143,158,172,186,199,212,223,232,238,241,243,245,251,259,268,
    279,291,306,323,342,363,387,414,443,476,513,554,600,650,706,766,833,
    904,981,1064,1155,1253,1360,1476,1602,1737,1881,2034,2195,2366,2548,
    2743,2952,3177,3420,3685,3975,4293,4643,5028,5454,5924,6443,7014,7637,
    8313,9040,9818,10647,11530,12471,13472,14537,15668,16859,18104,19395,
    20727,22091,23482,24894,26322,27760,29202,30642,32076,33496,34898,
    36275,37623,38935,40205,41429,42599,43712,44760,45738,46640,47462,
    48476,1187,757,447,3789,432,220,161,128,110,96,85,77,72,70,70,72,77,
    86,96,107,117,126,136,145,154,162,170,176,182,188,195,203,211,220,230,
    240,252,264,278,292,309,326,346,368,393,420,451,485,523,564,608,655,
    706,762,822,888,961,1040,1128,1224,1330,1446,1574,1714,1867,2035,2217,
    2419,2643,2893,3174,3489,3841,4233,4669,5150,5680,6259,6889,7569,8300,
    9083,9921,10819,11780,12809,13906,15070,16294,17573,18902,20276,21690,
    23141,24624,26136,27671,29226,30796,32379,33968,35561,37152,38739,
    40316,41880,43427,44951,46450,47919,49353,50750,52104,1541,1479,1365,
    9863,2036,959,516,389,330,284,249,225,212,207,212,224,248,285,329,372,
    411,445,479,511,537,555,560,554,543,535,537,548,566,588,613,637,660,
    683,708,734,763,796,832,872,916,963,1015,1072,1135,1203,1278,1361,
    1450,1544,1643,1745,1849,1960,2079,2212,2361,2529,2715,2914,3123,3340,
    3562,3794,4043,4314,4614,4949,5322,5736,6196,6704,7260,7854,8470,9093,
    9709,10309,10904,11511,12147,12828,13572,14392,15301,16312,17438,
    18694,20091,21643,23363,25264,27360,29664,32189,34948,37954,41221,
    44761,48588,52715,57155,61921,67028,72486,72486,72486,72486,72486,
    1270,1148,1016,8717,2038,986,546,393,310,251,212,190,183,189,203,225,
    256,299,349,400,448,495,545,593,636,670,693,709,720,726,731,733,730,
    725,719,712,706,701,697,696,699,706,717,732,751,773,799,828,864,907,
    959,1020,1091,1171,1260,1358,1464,1577,1696,1819,1947,2078,2213,2353,
    2499,2652,2814,2987,3171,3370,3585,3820,4080,4373,4708,5092,5527,5998,
    6485,6966,7421,7837,8232,8633,9066,9556,10131,10810,11615,12567,13686,
    14993,16508,18253,20249,22515,25073,27943,31146,34703,38635,42962,
    47705,52884,58521,64637,71251,78385,86059,86059,86059,86059,86059,
    1921,1688,1150,8228,937,432,269,216,186,163,147,137,134,138,149,167,
    194,231,274,320,369,422,483,544,602,650,685,711,733,754,780,810,840,
    872,906,943,983,1025,1071,1121,1174,1230,1293,1362,1440,1528,1629,
    1740,1859,1986,2118,2255,2394,2536,2679,2823,2966,3108,3248,3386,3520,
    3650,3779,3910,4045,4189,4343,4508,4685,4875,5077,5294,5532,5799,6104,
    6455,6857,7309,7803,8336,8902,9495,10107,10730,11353,11969,12573,
    13173,13783,14415,15083,15799,16574,17417,18340,19352,20463,21685,
    23027,24500,26113,27877,29802,31900,34178,36649,39322,42208,45317,
    48658,52244,56082,60185,1486,1230,974,6584,796,372,248,209,175,146,
    123,108,101,104,116,140,182,241,307,371,424,465,501,532,559,583,603,
    616,627,640,657,680,705,733,764,799,837,880,924,971,1020,1070,1123,
    1181,1246,1320,1405,1499,1602,1711,1824,1942,2062,2187,2315,2447,2583,
    2721,2858,2992,3121,3242,3358,3472,3586,3703,3825,3954,4090,4233,4384,
    4544,4718,4912,5129,5376,5655,5963,6294,6641,6998,7362,7737,8127,8540,
    8981,9457,9971,10529,11135,11793,12509,13287,14132,15048,16040,17112,
    18270,19517,20858,22299,23843,25496,27261,29143,31148,33280,35543,
    37941,40481,43165,46000,48988,1101,875,274,3069,211,139,106,90,81,75,
    68,63,60,60,62,67,76,90,104,120,132,144,153,162,170,174,175,174,171,
    168,169,173,176,182,190,201,215,229,249,269,294,322,355,392,432,477,
    526,579,638,700,771,847,925,1012,1107,1211,1328,1453,1588,1731,1882,
    2041,2207,2381,2569,2772,2985,3207,3445,3706,4000,4320,4662,5028,5434,
    5886,6385,6920,7501,8120,8788,9489,10211,10994,11851,12801,13877,
    15052,16305,17591,18897,20205,21518,22903,24284,25763,27303,28959,
    30553,32227,33911,35617,37357,39112,40866,42600,44312,46014,47709,
    49403,51100,52810,54520,832,613,219,2355,189,112,87,69,61,54,47,44,40,
    40,39,40,43,48,54,59,63,67,69,74,76,80,82,85,88,92,97,101,109,114,123,
    131,141,150,161,173,188,204,222,241,264,287,313,343,372,407,441,480,
    519,561,608,663,721,784,854,930,1019,1117,1224,1340,1468,1608,1751,
    1900,2063,2254,2490,2764,3068,3410,3786,4203,4650,5129,5651,6221,6850,
    7537,8269,9062,9910,10830,11815,12855,13969,15142,16407,17746,19168,
    20653,22232,23835,25571,27393,29261,31159,33050,34954,36895,38839,
    40752,42600,44367,46076,47750,49417,51100,52810,54529,1449,1154,570,
    5089,466,260,176,144,124,106,96,87,83,84,89,99,115,137,164,192,220,
    249,282,314,344,369,387,399,409,420,435,452,471,492,515,543,574,608,
    646,687,732,778,826,879,940,1011,1093,1185,1285,1394,1510,1635,1768,
    1909,2058,2217,2388,2571,2762,2953,3139,3320,3498,3676,3855,4036,4216,
    4394,4576,4765,4967,5178,5393,5620,5867,6140,6441,6764,7108,7474,7862,
    8256,8657,9086,9566,10117,10692,11278,11944,12761,13800,15093,16595,
    18255,20022,21845,23758,25796,27907,30040,32147,34259,36411,38552,
    40632,42600,44425,46141,47794,49431,51100,52810,54529,1123,863,465,
    4087,388,215,166,127,107,90,76,65,58,55,57,64,79,100,125,150,173,192,
    210,227,244,261,276,290,303,318,334,352,370,390,413,439,470,504,542,
    582,626,671,719,770,826,890,962,1042,1127,1216,1308,1402,1498,1599,
    1707,1824,1953,2093,2239,2387,2532,2673,2813,2954,3098,3245,3396,3548,
    3704,3865,4033,4203,4373,4553,4748,4969,5216,5486,5773,6074,6385,6692,
    7001,7327,7689,8102,8490,8840,9270,9898,10840,12145,13734,15535,17477,
    19489,21619,23913,26300,28710,31070,33427,35831,38208,40489,42600,
    44494,46218,47845,49448,51100,52810,54529,1449,1154,570,5089,466,260,
    176,144,124,106,96,87,83,84,89,99,115,137,164,192,220,249,282,314,344,
    369,387,399,409,420,435,452,471,492,515,543,574,608,646,687,732,778,
    826,879,940,1011,1093,1185,1285,1394,1510,1635,1768,1909,2058,2217,
    2388,2571,2762,2953,3139,3320,3498,3676,3855,4036,4216,4394,4576,4765,
    4967,5178,5393,5620,5867,6140,6441,6764,7108,7474,7862,8256,8657,9086,
    9566,10117,10692,11278,11944,12761,13800,15093,16595,18255,20022,
    21845,23758,25796,27907,30040,32147,34259,36411,38552,40632,42600,
    44425,46141,47794,49431,51100,52810,54529,1123,863,465,4087,388,215,
    166,127,107,90,76,65,58,55,57,64,79,100,125,150,173,192,210,227,244,
    261,276,290,303,318,334,352,370,390,413,439,470,504,542,582,626,671,
    719,770,826,890,962,1042,1127,1216,1308,1402,1498,1599,1707,1824,1953,
    2093,2239,2387,2532,2673,2813,2954,3098,3245,3396,3548,3704,3865,4033,
    4203,4373,4553,4748,4969,5216,5486,5773,6074,6385,6692,7001,7327,7689,
    8102,8490,8840,9270,9898,10840,12145,13734,15535,17477,19489,21619,
    23913,26300,28710,31070,33427,35831,38208,40489,42600,44494,46218,
    47845,49448,51100,52810,54529,1079,709,191,2592,153,101,81,69,62,57,
    53,49,45,42,42,47,59,75,93,111,126,139,149,159,169,174,172,165,156,
    149,145,145,149,156,163,171,181,193,207,225,246,270,299,332,368,409,
    454,504,558,617,686,766,856,955,1058,1162,1264,1368,1475,1593,1730,
    1891,2074,2271,2476,2690,2912,3143,3389,3652,3930,4225,4538,4871,5230,
    5623,6060,6542,7066,7636,8271,8986,9788,10732,11799,12895,13920,14861,
    16039,17303,18665,20194,21877,23601,25289,26973,28612,30128,31416,
    32915,34450,36018,37616,39242,40891,42562,44250,45951,47662,49378,
    51095,52810,54519,828,499,146,1964,135,81,63,55,47,41,37,33,30,28,28,
    29,32,36,41,47,51,54,55,56,58,60,62,63,65,68,71,74,79,85,91,97,105,
    113,122,133,145,158,174,190,209,229,252,276,303,331,362,396,432,473,
    517,560,601,642,687,740,805,886,981,1088,1203,1325,1454,1592,1742,
    1909,2100,2319,2567,2836,3129,3462,3845,4278,4742,5245,5827,6509,7294,
    8213,9231,10264,11235,12151,13625,15237,16936,18731,20611,22560,24536,
    26481,28322,29988,31416,32915,34450,36018,37616,39242,40891,42562,
    44250,45951,47662,49378,51095,52810,54519,1594,988,443,4699,337,197,
    134,102,87,76,68,63,60,60,64,72,85,102,120,140,162,186,210,236,262,
    283,298,307,316,327,339,353,370,389,409,431,455,483,513,546,585,633,
    688,749,814,875,931,984,1038,1101,1183,1292,1422,1565,1710,1854,1994,
    2131,2273,2427,2589,2762,2947,3137,3335,3554,3801,4072,4365,4665,4953,
    5213,5448,5690,5944,6177,6375,6548,6673,6803,7037,7460,8065,8836,9668,
    10452,11038,11410,12280,13313,14588,16219,18166,20304,22519,24791,
    27050,29270,31416,32915,34450,36018,37616,39242,40891,42562,44250,
    45951,47662,49378,51095,52810,54519,1272,742,385,3828,289,160,115,92,
    77,66,56,49,43,40,39,41,46,53,63,73,84,94,104,116,128,140,150,160,171,
    182,197,214,234,256,279,303,326,350,374,402,434,471,514,561,611,656,
    696,733,769,814,875,957,1058,1167,1279,1392,1504,1617,1731,1852,1983,
    2130,2287,2459,2632,2784,2901,2995,3072,3162,3296,3500,3762,4066,4372,
    4646,4853,5007,5127,5274,5510,5897,6420,7060,7731,8349,8813,9127,
    10205,11467,12972,14790,16879,19137,21496,23959,26478,28995,31416,
    32915,34450,36018,37616,39242,40891,42562,44250,45951,47662,49378,
    51095,52810,54519,1594,988,443,4699,337,197,134,102,87,76,68,63,60,60,
    64,72,85,102,120,140,162,186,210,236,262,283,298,307,316,327,339,353,
    370,389,409,431,455,483,513,546,585,633,688,749,814,875,931,984,1038,
    1101,1183,1292,1422,1565,1710,1854,1994,2131,2273,2427,2589,2762,2947,
    3137,3335,3554,3801,4072,4365,4665,4953,5213,5448,5690,5944,6177,6375,
    6548,6673,6803,7037,7460,8065,8836,9668,10452,11038,11410,12280,13313,
    14588,16219,18166,20304,22519,24791,27050,29270,31416,32915,34450,
    36018,37616,39242,40891,42562,44250,45951,47662,49378,51095,52810,
    54519,1272,742,385,3828,289,160,115,92,77,66,56,49,43,40,39,41,46,53,
    63,73,84,94,104,116,128,140,150,160,171,182,197,214,234,256,279,303,
    326,350,374,402,434,471,514,561,611,656,696,733,769,814,875,957,1058,
    1167,1279,1392,1504,1617,1731,1852,1983,2130,2287,2459,2632,2784,2901,
    2995,3072,3162,3296,3500,3762,4066,4372,4646,4853,5007,5127,5274,5510,
    5897,6420,7060,7731,8349,8813,9127,10205,11467,12972,14790,16879,
    19137,21496,23959,26478,28995,31416,32915,34450,36018,37616,39242,
    40891,42562,44250,45951,47662,49378,51095,52810,54519,892,527,139,
    2006,116,83,72,59,54,51,48,44,39,34,32,39,55,80,107,134,156,172,181,
    190,201,205,203,195,184,173,165,162,165,170,176,183,192,203,217,235,
    256,281,310,340,372,409,452,501,555,612,673,739,812,892,980,1081,1194,
    1318,1452,1594,1745,1906,2077,2258,2451,2657,2879,3120,3386,3674,3977,
    4284,4597,4916,5262,5655,6118,6647,7231,7843,8472,9103,9749,10466,
    11273,12127,13012,13942,15033,16321,17666,18947,20145,21344,22684,
    24152,25767,27426,29014,30431,31784,33085,34324,35479,36553,37550,
    38471,39320,40101,40818,41475,42075,42624,686,373,115,1532,101,67,54,
    47,40,36,32,29,27,25,24,26,31,38,46,55,61,64,64,64,65,65,66,67,68,70,
    72,75,79,84,90,97,104,113,122,133,146,161,177,193,211,230,254,280,308,
    336,366,397,430,466,507,550,596,646,699,758,819,884,953,1027,1110,
    1203,1309,1429,1563,1713,1883,2075,2288,2513,2759,3048,3396,3803,4255,
    4740,5264,5829,6440,7128,7893,8702,9539,10427,11465,12685,13944,15144,
    16303,17570,19046,20617,22206,23758,25298,26762,28133,29413,30615,
    31742,32794,33772,34679,35517,36289,36999,37651,38248,38793,1522,742,
    272,3408,217,155,109,94,82,74,67,60,53,48,48,58,80,113,151,190,230,
    270,309,357,410,452,473,475,468,464,464,474,494,515,535,558,587,621,
    657,697,742,791,843,898,955,1016,1080,1149,1222,1298,1381,1472,1573,
    1683,1802,1927,2054,2182,2314,2453,2602,2763,2939,3127,3324,3532,3744,
    3959,4171,4389,4636,4935,5292,5714,6169,6617,7001,7318,7636,8001,8350,
    8666,8942,9160,9353,9587,9937,10418,11257,12156,13089,13980,14832,
    15687,16620,17656,18827,20064,21270,21795,22278,22723,23132,23506,
    23848,24160,24445,24705,24941,25155,25350,25526,25686,1170,592,236,
    2765,189,130,90,68,62,52,45,39,35,33,33,36,44,54,67,80,93,103,111,121,
    132,141,150,157,164,173,183,195,210,225,242,262,286,313,343,373,405,
    438,472,507,542,581,624,672,725,779,836,892,950,1013,1081,1153,1229,
    1308,1392,1483,1581,1689,1809,1937,2073,2226,2392,2566,2738,2909,3093,
    3308,3561,3863,4194,4519,4790,5004,5208,5446,5704,5997,6321,6656,6991,
    7342,7716,8122,8747,9465,10282,11201,12222,13355,14548,15672,16605,
    17401,18220,18719,19180,19605,19996,20355,20684,20985,21259,21510,
    21738,21945,22134,22305,22460,1431,694,257,3408,217,155,109,94,82,74,
    67,60,53,48,48,58,80,113,151,190,230,270,309,357,410,452,473,475,468,
    464,464,474,494,515,535,558,587,621,657,697,742,791,843,898,955,1016,
    1080,1149,1222,1298,1381,1472,1573,1683,1802,1927,2054,2182,2314,2453,
    2602,2763,2939,3127,3324,3532,3744,3959,4171,4389,4636,4935,5292,5714,
    6169,6617,7001,7318,7636,8001,8350,8666,8942,9160,9353,9587,9937,
    10418,11257,12156,13089,13980,14832,15687,16620,17656,18827,20064,
    21270,21795,22278,22723,23132,23506,23848,24160,24445,24705,24941,
    25155,25350,25526,25686,1097,555,223,2765,189,130,90,68,62,52,45,39,
    35,33,33,36,44,54,67,80,93,103,111,121,132,141,150,157,164,173,183,
    195,210,225,242,262,286,313,343,373,405,438,472,507,542,581,624,672,
    725,779,836,892,950,1013,1081,1153,1229,1308,1392,1483,1581,1689,1809,
    1937,2073,2226,2392,2566,2738,2909,3093,3308,3561,3863,4194,4519,4790,
    5004,5208,5446,5704,5997,6321,6656,6991,7342,7716,8122,8747,9465,
    10282,11201,12222,13355,14548,15672,16605,17401,18220,18719,19180,
    19605,19996,20355,20684,20985,21259,21510,21738,21945,22134,22305,
    22460,438,256,139,1231,92,66,53,43,39,37,34,30,24,19,19,28,46,71,96,
    118,137,151,163,175,186,193,193,189,183,177,172,168,167,166,165,166,
    169,175,184,196,209,224,240,261,287,316,348,382,420,463,514,573,639,
    706,775,850,934,1027,1125,1227,1338,1464,1605,1762,1933,2119,2316,
    2523,2738,2968,3218,3495,3805,4148,4516,4901,5295,5703,6146,6642,7180,
    7762,8394,9099,9886,10733,11613,12523,13507,14592,15691,16774,17875,
    19058,20389,21864,23453,25061,26617,28001,29311,30545,31703,32784,
    33791,34724,35588,36384,37117,37790,38407,38971,39486,361,190,111,965,
    77,51,37,30,28,26,23,21,18,17,16,19,25,32,40,47,52,54,55,56,57,57,58,
    58,58,58,59,60,63,65,68,72,77,83,90,99,109,119,130,143,158,174,192,
    211,231,254,280,310,343,376,410,447,488,532,579,628,681,742,811,889,
    975,1067,1162,1259,1359,1470,1595,1740,1907,2092,2294,2517,2760,3027,
    3315,3637,4015,4467,4995,5589,6239,6949,7713,8539,9463,10491,11534,
    12559,13617,14831,16231,17709,19198,20690,22228,23729,25173,26551,
    27859,29094,30255,31342,32355,33297,34168,34973,35715,36397,37022,869,
    421,230,2061,139,101,82,66,58,51,45,39,34,30,31,39,55,76,98,119,140,
    162,186,212,239,262,279,291,302,314,325,335,346,356,367,379,395,413,
    436,461,491,523,557,595,638,687,743,807,877,952,1036,1128,1225,1323,
    1424,1531,1648,1774,1905,2039,2174,2312,2459,2619,2794,2981,3172,3361,
    3545,3733,3936,4171,4445,4754,5084,5421,5742,6046,6356,6699,7083,7538,
    8088,8772,9578,10433,11190,11781,12406,13154,13945,14805,15729,16621,
    17527,18599,19866,21229,22554,23274,23944,24563,25135,25662,26146,
    26590,26996,27367,27706,28014,28295,28550,28782,741,331,207,1739,120,
    80,63,46,41,34,29,26,24,23,24,26,31,36,43,49,55,60,66,72,78,84,90,96,
    102,109,115,121,127,133,140,148,159,172,187,205,224,245,266,290,316,
    345,378,413,451,492,537,585,636,688,740,796,858,927,1001,1079,1161,
    1247,1340,1441,1552,1668,1785,1898,2007,2120,2253,2422,2631,2875,3138,
    3408,3659,3888,4114,4363,4648,5001,5442,5992,6626,7279,7834,8251,8685,
    9238,9881,10652,11547,12514,13529,14624,15791,17016,18279,19170,20022,
    20825,21577,22279,22930,23534,24091,24605,25077,25510,25907,26269,
    26600,768,372,206,2297,148,110,86,70,63,55,49,43,37,33,34,41,57,78,99,
    120,142,165,191,221,251,279,300,315,330,346,362,377,392,408,424,441,
    460,483,509,539,572,609,648,691,739,794,857,929,1007,1090,1181,1280,
    1384,1488,1594,1709,1835,1972,2116,2262,2408,2556,2711,2877,3058,3252,
    3452,3651,3846,4044,4260,4511,4804,5141,5501,5866,6202,6508,6814,7154,
    7537,7999,8566,9268,10087,10953,11714,12302,12872,13559,14282,15071,
    15928,16761,17617,18648,19888,21236,22554,23274,23944,24563,25135,
    25662,26146,26590,26996,27367,27706,28014,28295,28550,28782,660,298,
    186,1927,127,87,66,48,44,37,31,27,25,24,24,27,31,37,43,49,56,62,68,74,
    81,88,95,102,109,118,126,133,140,148,157,168,180,194,211,231,252,275,
    298,324,352,385,421,462,505,552,602,655,710,765,821,882,950,1026,1107,
    1192,1280,1372,1470,1577,1695,1817,1936,2050,2158,2272,2408,2587,2810,
    3072,3354,3639,3899,4132,4360,4615,4909,5282,5754,6350,7041,7751,8335,
    8744,9106,9591,10168,10886,11738,12656,13619,14672,15816,17027,18279,
    19170,20022,20825,21577,22279,22930,23534,24091,24605,25077,25510,
    25907,26269,26600,302,134,100,862,66,49,37,32,28,26,24,22,19,16,17,24,
    39,59,81,102,118,127,132,136,141,145,148,150,151,153,156,162,169,177,
    185,193,201,210,219,230,240,250,260,271,283,298,317,341,370,404,441,
    479,518,564,620,683,753,831,913,1004,1109,1231,1366,1503,1641,1788,
    1947,2118,2297,2483,2689,2926,3200,3509,3848,4215,4598,4993,5414,5875,
    6372,6920,7533,8246,9049,9891,10715,11519,12436,13522,14695,15927,
    17219,18617,20159,21773,23376,24893,26329,27914,29399,30869,32413,
    34033,35735,37522,39398,41368,43436,45608,47888,50282,52797,249,101,
    77,667,59,37,29,23,21,19,17,16,15,14,14,17,21,27,34,40,45,47,48,49,50,
    51,51,51,51,51,53,55,58,62,66,70,74,78,82,88,94,102,111,121,131,143,
    157,173,193,215,240,265,291,321,356,394,434,476,521,571,628,693,764,
    837,912,993,1081,1177,1278,1383,1500,1634,1791,1969,2168,2386,2618,
    2860,3111,3387,3707,4090,4542,5053,5606,6215,6878,7607,8445,9402,
    10431,11512,12673,14015,15536,17101,18626,20148,21737,23434,25091,
    26715,28318,30017,31818,33727,35750,37895,40169,42579,45134,47842,
    50712,681,227,176,1712,121,77,63,48,42,38,34,30,24,20,20,31,55,87,121,
    153,180,201,217,234,252,267,276,283,287,292,301,315,333,351,369,389,
    411,435,462,489,516,538,558,579,603,631,665,705,752,805,864,926,991,
    1061,1138,1225,1326,1441,1566,1698,1838,1984,2134,2284,2440,2613,2809,
    3025,3250,3475,3706,3950,4218,4524,4867,5230,5586,5921,6251,6603,6988,
    7444,7989,8644,9370,10106,10749,11273,11827,12507,13318,14333,15560,
    16977,18502,19999,21198,22061,22903,24048,25250,26513,27838,29230,
    30692,32226,33837,35529,37306,39171,41130,43186,45345,575,174,148,
    1428,103,62,47,35,34,29,25,22,20,19,20,22,26,31,37,43,49,55,59,65,70,
    76,82,88,94,100,107,116,126,137,147,157,167,178,189,200,213,226,239,
    254,271,289,310,335,363,396,433,473,515,561,612,665,721,781,843,910,
    987,1077,1176,1278,1381,1491,1610,1740,1877,2018,2162,2312,2473,2650,
    2851,3075,3315,3564,3817,4082,4377,4726,5146,5647,6209,6803,7370,7896,
    8452,9105,9831,10667,11633,12752,13979,15196,16268,17245,18338,19682,
    21089,22557,23911,25346,26866,28478,30187,31998,33918,35953,38110,
    40397,42821,811,263,201,1977,134,85,69,54,48,43,39,34,27,22,23,35,62,
    97,137,173,205,231,252,274,298,317,330,338,344,351,363,382,407,433,
    458,485,513,544,577,612,645,675,702,730,762,799,841,890,946,1010,1081,
    1155,1232,1314,1404,1504,1619,1748,1885,2029,2181,2340,2506,2672,2842,
    3033,3250,3489,3739,3989,4244,4511,4801,5131,5500,5885,6255,6599,6931,
    7285,7675,8145,8713,9403,10169,10937,11578,12062,12515,13086,13796,
    14736,15912,17285,18751,20162,21224,21915,22659,23792,24982,26231,
    27542,28920,30365,31884,33478,35152,36909,38755,40693,42727,44864,675,
    200,169,1644,113,69,51,40,37,32,27,24,23,22,23,25,28,33,39,45,51,58,
    64,72,80,88,96,103,110,118,127,138,151,165,178,191,203,216,229,243,
    259,275,293,313,335,359,384,411,442,478,519,564,612,664,721,781,843,
    910,979,1054,1142,1247,1363,1482,1602,1727,1859,1997,2142,2291,2442,
    2601,2774,2965,3177,3411,3656,3904,4152,4412,4705,5061,5492,6008,6581,
    7179,7744,8264,8797,9428,10129,10933,11860,12942,14138,15318,16333,
    17232,18244,19556,20946,22414,23758,25184,26695,28297,29994,31794,
    33702,35724,37867,40139,42548)

temp <- array(temp/100000, dim=c(113,2,3,6))
 
#
# At this point temp is an array of age, sex, race, year.  
#  The first 3 ages are the q's for days 0-1, 1-7, and 7-28, the fourth
#  is the q for the entire first year.
# Change the array to one of daily hazard rates
#  For the 4th row, make it so the sum of the first year's hazard is correct,
#  i.e., 1*row1 + 6*row2 + 21*row3 + 337.24* row4 = -log(1-q)
 
temp2 <- -log(1- temp)/365.24    
temp2[1,,,] <- -log(1-temp[1,,,]) /1
temp2[2,,,] <- -log(1-temp[2,,,]) /6     #days 1-7
temp2[3,,,] <- -log(1-temp[3,,,]) /21    #days 7-28
temp2[4,,,] <- (-log(1-temp[4,,,]) -
		 (temp2[1,,,] + 6*temp2[2,,,] + 21*temp2[3,,,]))/ 337.24

#
# Now, add in the year 2000 extrapolation
#
survexp.usr <- array(0, dim=c(113,2,3,7))
survexp.usr[,,,1:6] <- temp2
data.restore('survexp2000.sdump')
for (i in 1:4) {
    survexp.usr[i,,,7] <- exp(log(temp2[i,,,6]) + survexp.2000[1,,1:3])
    }
survexp.usr[5:113,,,7] <- exp(log(temp2[5:113,,,6]) + survexp.2000[-1,,1:3])

attributes(survexp.usr) <- list (
	dim      =c(113,2,3,7),
	dimnames =list(c('0-1d','1-7d', '7-28d', '28-365d', 
	                                         as.character(1:109)),
	      c("male", "female"), c('white', 'nonwhite', 'black'),
	      10*(194:200)),
	dimid    =c("age", "sex", "race", "year"),
	factor   =c(0,1,1,10),
	cutpoints=list(c(0,1,7,28,1:109 * 365.24), NULL, NULL,
	               mdy.date(1,1, (194:200)*10)),
	summary = function(R) {
		     x <- c(format(round(min(R[,1]) /365.24, 1)),
			    format(round(max(R[,1]) /355.24, 1)),
			    sum(R[,2]==1), sum(R[,2]==2),
			    sum(R[,3]==1), sum(R[,3]==2), sum(R[,3]==3))
		     x2<-as.character(as.date(c(min(R[,4]), max(R[,4]))))

		     paste("  age ranges from", x[1], "to", x[2], "years\n",
			   " male:", x[3], " female:", x[4], "\n",
			   " date of entry from", x2[1], "to", x2[2], "\n",
			   " white:",x[7], " nonwhite:", x[8]," black:", x[9],
			   "\n")
		     })
rm(temp, temp2, survexp.2000)
oldClass(survexp.usr) <- 'ratetable'
