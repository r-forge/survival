# SCCS $Id: survexp.uswhite.s,v 4.7 1995-02-06 22:27:55 therneau Exp $
#
# Create the US total hazards table, whites only
#   The raw numbers below are q* 10^5.  Note that there are 24 leap years/100
#
survexp.uswhite  <- {
    temp1 <- c(
     3069,211,139,106,90,81,75,68,63,60,60,62,67,76,90,104,120,132,144,153,162,
     170,174,175,174,171,168,169,173,176,182,190,201,215,229,249,269,294,322,
     355,392,432,477,526,579,638,700,771,847,925,1012,1107,1211,1328,1453,1588,
     1731,1882,2041,2207,2381,2569,2772,2985,3207,3445,3706,4000,4320,4662,
     5028,5434,5886,6385,6920,7501,8120,8788,9489,10211,10994,11851,12801,
     13877,15052,16305,17591,18897,20205,21518,22903,24284,25763,27303,28959,
     30553,32227,33911,35617,37357,39112,40866,42600,44312,46014,47709,49403,
     51100,52810,54520,2355,189,112,87,69,61,54,47,44,40,40,39,40,43,48,54,59,
     63,67,69,74,76,80,82,85,88,92,97,101,109,114,123,131,141,150,161,173,188,
     204,222,241,264,287,313,343,372,407,441,480,519,561,608,663,721,784,854,
     930,1019,1117,1224,1340,1468,1608,1751,1900,2063,2254,2490,2764,3068,3410,
     3786,4203,4650,5129,5651,6221,6850,7537,8269,9062,9910,10830,11815,12855,
     13969,15142,16407,17746,19168,20653,22232,23835,25571,27393,29261,31159,
     33050,34954,36895,38839,40752,42600,44367,46076,47750,49417,51100,52810,
     54529)
    temp2 <- c(
     2592,153,101,81,69,62,57,53,49,45,42,42,47,59,75,93,111,126,139,149,159,
     169,174,172,165,156,149,145,145,149,156,163,171,181,193,207,225,246,270,
     299,332,368,409,454,504,558,617,686,766,856,955,1058,1162,1264,1368,1475,
     1593,1730,1891,2074,2271,2476,2690,2912,3143,3389,3652,3930,4225,4538,
     4871,5230,5623,6060,6542,7066,7636,8271,8986,9788,10732,11799,12895,13920,
     14861,16039,17303,18665,20194,21877,23601,25289,26973,28612,30128,31416,
     32915,34450,36018,37616,39242,40891,42562,44250,45951,47662,49378,51095,
     52810,54519,1964,135,81,63,55,47,41,37,33,30,28,28,29,32,36,41,47,51,54,
     55,56,58,60,62,63,65,68,71,74,79,85,91,97,105,113,122,133,145,158,174,190,
     209,229,252,276,303,331,362,396,432,473,517,560,601,642,687,740,805,886,
     981,1088,1203,1325,1454,1592,1742,1909,2100,2319,2567,2836,3129,3462,3845,
     4278,4742,5245,5827,6509,7294,8213,9231,10264,11235,12151,13625,15237,
     16936,18731,20611,22560,24536,26481,28322,29988,31416,32915,34450,36018,
     37616,39242,40891,42562,44250,45951,47662,49378,51095,52810,54519,2006,
     116,83,72,59,54,51,48,44,39,34,32,39,55,80,107,134,156,172,181,190,201,
     205,203,195,184,173,165,162,165,170,176,183,192,203,217,235,256,281,310,
     340,372,409,452,501,555,612,673,739,812,892,980,1081,1194,1318,1452,1594,
     1745,1906,2077,2258,2451,2657,2879,3120,3386,3674,3977,4284,4597,4916,
     5262,5655,6118,6647,7231,7843,8472,9103,9749,10466,11273,12127,13012,
     13942,15033,16321,17666,18947,20145,21344,22684,24152,25767,27426,29014,
     30431,31784,33085,34324,35479,36553,37550,38471,39320,40101,40818,41475,
     42075,42624,1532,101,67,54,47,40,36,32,29,27,25,24,26,31,38,46,55,61,64,
     64,64,65,65,66,67,68,70,72,75,79,84,90,97,104,113,122,133,146,161,177,193,
     211,230,254,280,308,336,366,397,430,466,507,550,596,646,699,758,819,884,
     953,1027,1110,1203,1309,1429,1563,1713,1883,2075,2288,2513,2759,3048,3396,
     3803,4255,4740,5264,5829,6440,7128,7893,8702,9539,10427,11465,12685,13944,
     15144,16303,17570,19046,20617,22206,23758,25298,26762,28133,29413,30615,
     31742,32794,33772,34679,35517,36289,36999,37651,38248,38793,1231,92,66,53,
     43,39,37,34,30,24,19,19,28,46,71,96,118,137,151,163,175,186,193,193,189,
     183,177,172,168,167,166,165,166,169,175,184,196,209,224,240,261,287,316,
     348,382,420,463,514,573,639,706,775,850,934,1027,1125,1227,1338,1464,1605,
     1762,1933,2119,2316,2523,2738,2968,3218,3495,3805,4148,4516,4901,5295,
     5703,6146,6642,7180,7762,8394,9099,9886,10733,11613,12523,13507,14592,
     15691,16774,17875,19058,20389,21864,23453,25061,26617,28001,29311,30545,
     31703,32784,33791,34724,35588,36384,37117,37790,38407,38971,39486,965,77,
     51,37,30,28,26,23,21,18,17,16,19,25,32,40,47,52,54,55,56,57,57,58,58,58,
     58,59,60,63,65,68,72,77,83,90,99,109,119,130,143,158,174,192,211,231,254,
     280,310,343,376,410,447,488,532,579,628,681,742,811,889,975,1067,1162,
     1259,1359,1470,1595,1740,1907,2092,2294,2517,2760,3027,3315,3637,4015,
     4467,4995,5589,6239,6949,7713,8539,9463,10491,11534,12559,13617,14831,
     16231,17709,19198,20690,22228,23729,25173,26551,27859,29094,30255,31342,
     32355,33297,34168,34973,35715,36397,37022)

    temp2 <- -log(1- c(temp1,temp2)/100000)/365.24    #daily hazard rate

    #Add in the extrapolated data
    temp <- array(0, c(110,2,6))
    temp[,,1:4] <- temp2
    fix  <- c( .00061*(0:109) - .1271, .00041*(0:109) - .1770)
    temp[,,5]   <- exp(log(temp[,,4]) + fix)
    temp[,,6]   <- exp(log(temp[,,5]) + fix)

    attributes(temp) <- list (
	dim      =c(110,2,6),
	dimnames =list(0:109, c("male", "female"), 10 * 195:200),
	dimid    =c("age", "sex", "year"),
	factor   =c(0,1,10),
	cutpoints=list(0:109 * 365.24, NULL, mdy.date(1,1, (195:200)*10)),
	summary = function(R) {
		     x <- c(format(round(min(R[,1]) /365.24, 1)),
			    format(round(max(R[,1]) /355.24, 1)),
			    sum(R[,2]==1), sum(R[,2]==2),
			    min(R[,3]), max(R[,3]))

		     paste("age ranges from", x[1], "to", x[2], "years,\n",
			   "male:", x[3], "female:", x[4], "\n",
			   "year of entry from", x[5], "to", x[6], "\n")
		     },
	class='ratetable')
    temp
    }
