#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library(RUnit)
library(Anaquin)

.test.data <- function()
{
    input <- c(15.1062,
               15.1062,
               966.797,
               241.699,
               30.2124,
               7734.38,
               483.398,
               30937.5,
               483.398,
               15.1062,
               483.398,
               241.699,
               60.4248,
               0.118017,
               7734.37,
               7734.38,
               120.85,
               1933.59,
               0.944138,
               120.85,
               7.5531,
               3.77655,
               3867.19,
               15468.8,
               1.88828,
               1933.59,
               120.85,
               3867.19,
               30.2124,
               0.472069,
               241.699,
               60.4248,
               0.944138,
               120.85,
               1.88828,
               60.4248,
               15468.8,
               60.4248,
               1933.59,
               1.88828,
               30.2124,
               0.944138,
               3867.19,
               15468.8,
               3867.19,
               30.2124,
               15.1062,
               1933.59,
               7.5531,
               0.472069,
               0.472069,
               0.236034,
               0.0590086,
               241.699,
               0.236034,
               0.472069,
               0.118017,
               120.85,
               0.118017,
               483.398,
               30937.5,
               7.5531,
               0.472069,
               483.398,
               7.5531,
               966.797,
               3.77655,
               30937.5,
               3.77655,
               3.77655,
               966.797,
               3.77655,
               0.236034,
               1.88828,
               0.0590086)
    
    o1 <- c(4.29354,
            6.41774,
            271.696,
            67.3212,
            3.1943,
            2628.81,
            245.034,
            14727.4,
            38.2906,
            8.41748,
            140.367,
            90.1342,
            5.03172,
            0.282232,
            2340.3,
            5418.16,
            124.648,
            157.443,
            0.437077,
            2.42797,
            3.88121,
            NA,
            2010.24,
            7762.18,
            2.05478,
            526.579,
            29.9481,
            1369.19,
            9.21093,
            0.305795,
            50.4639,
            11.4579,
            0.291139,
            76.4986,
            3.2669,
            21.1347,
            8064.47,
            37.4146,
            438.561,
            0.38055,
            11.2435,
            0.88877,
            1491.72,
            4786.46,
            867.185,
            7.30535,
            2.78194,
            559.819,
            3.69124,
            0.523401,
            0.474327,
            0.929315,
            0.722466,
            111.985,
            1.1903,
            1.47662,
            0.285894,
            35.2424,
            0.700978,
            419.301,
            25760.2,
            2.5859,
            1.31373,
            24.5161,
            4.35937,
            457.339,
            2.59279,
            15205,
            0.64947,
            1.1555,
            261.917,
            0.252222,
            0.443106,
            0.376111,
            1.0948)
    
    o2 <- c(3.09828,
            8.68183,
            272.248,
            80.5532,
            3.17143,
            2813.75,
            252.489,
            15776.6,
            55.5334,
            7.92421,
            134.694,
            92.8906,
            8.70729,
            0.88829,
            2313.62,
            5647.12,
            121.991,
            157.827,
            0.411342,
            2.28061,
            2.82815,
            2.02179,
            2007.42,
            8399.07,
            1.32179,
            479.284,
            30.1469,
            1354.44,
            7.04809,
            0.688492,
            45.6092,
            16.2209,
            0.655503,
            64.66,
            2.49032,
            23.6616,
            8801.15,
            37.5889,
            402.614,
            0.792792,
            13.1838,
            1.30589,
            1546.73,
            5024.49,
            870.876,
            6.57147,
            3.05968,
            568.894,
            3.72397,
            0.650601,
            NA,
            0.970407,
            0.902928,
            107.726,
            NA,
            1.74165,
            0.15758,
            42.8825,
            0.632109,
            436.377,
            27119.4,
            4.50098,
            0.297246,
            18.5351,
            2.64522,
            490.524,
            1.45904,
            15831,
            0.729047,
            0.449076,
            277.903,
            2.34995,
            0.385306,
            0.864491,
            0.944744)
    
    o3 <- c(5.78108,
            8.41477,
            261.53,
            71.0515,
            2.22543,
            2554.46,
            235.358,
            14074.1,
            52.9666,
            5.93463,
            130.117,
            89.6172,
            13.3467,
            2.71212,
            2186.4,
            5142.73,
            122.094,
            202.948,
            3.26205,
            3.628,
            2.82539,
            1.14586,
            1921.63,
            7459.34,
            0.987163,
            469.363,
            27.2883,
            1303.06,
            10.1981,
            1.0284,
            44.2452,
            14.9583,
            1.95817,
            62.0998,
            4.66003,
            24.9048,
            7705.65,
            39.1692,
            377.697,
            1.70661,
            12.2224,
            2.68944,
            1432.93,
            4672.95,
            821.054,
            6.00481,
            2.91698,
            549.193,
            6.69091,
            2.97898,
            1.91727,
            3.08212,
            2.09839,
            113.618,
            2.05095,
            4.40774,
            1.82189,
            36.1417,
            4.30722,
            410.602,
            23843.5,
            6.58454,
            1.3233,
            22.6614,
            6.2324,
            443.329,
            1.81664,
            14447.8,
            2.18542,
            0.148822,
            237.506,
            2.7184,
            2.3283,
            2.68241,
            3.31801)
    
    return (data.frame(Input=input, Observed1=o1, Observed2=o2, Observed3=o3))
}

test.PlotLinear_1 <- function()
{
    x <- .test.data()
    r <- plotLinear(row.names(x),
                    log2(x$Input),
                    log2(x$Observed1))
    
    checkEqualsNumeric(r$LOQ$model$breaks, 0.9170727)
}

test.PlotLinear_2 <- function()
{
    x <- .test.data()
    r <- plotLinear(row.names(x),
                    log2(x$Input),
                    log2(x[,c(2:4)]))

    checkEqualsNumeric(r$LOQ$model$breaks, 1.917068888)
}

test.PlotLinear_3 <- function()
{
    x <- "15.1062,
          15.1062,
          241.699,
          0.472069,
          120.85,
          483.398,
          0.0590086,
          483.398,
          15.1062,
          241.699,
          60.4248,
          120.85,
          120.85,
          0.472069,
          0.118017,
          0.0590086,
          0.944138,
          0.118017,
          7.5531,
          3.77655,
          0.236034,
          60.4248,
         1.88828,
          241.699,
          0.118017,
          0.236034,
          0.472069,
          30.2124,
          241.699,
          60.4248,
          0.944138,
          0.944138,
          0.118017,
          1.88828,
          60.4248,
          60.4248,
          0.944138,
          1.88828,
          0.472069,
          0.944138,
          15.1062,
          7.5531,
          15.1062,
          0.472069,
          15.1062,
          60.4248,
          7.5531,
          30.2124,
          30.2124,
          0.0590086,
          0.236034,
          0.0590086,
          241.699,
          0.236034,
          30.2124,
          120.85,
          0.118017,
          120.85,
          1.88828,
          7.5531,
          30.2124,
          483.398,
          7.5531,
          966.797,
          3.77655,
          241.699,
          3.77655,
          3.77655,
          3.77655,
          0.236034,
          1.88828,
          0.0590086"
    
    y <- "2292.02,
          1728.11,
          33095.9,
          51.557,
          13616.2,
          54615.3,
          5.75489,
          60652.2,
          1205.27,
          30255,
          7645.25,
          11118.4,
          16529.2,
          48.5283,
          12.6669,
          8.29854,
          96.6785,
          14.3095,
          1151.03,
          481.086,
          22.5602,
          9128.84,
          203.783,
          24093.9,
          15.2307,
          22.7496,
          50.3159,
          4104.9,
          19207.5,
          6752.21,
          123.603,
          187.043,
          14.0419,
          263.462,
          3310.3,
          8270.03,
          119.862,
          280.256,
          1e-05,
          115.157,
          2042.86,
          1084.78,
          1773.62,
          54.5923,
          1733.97,
          6365.94,
          901.141,
          4224.51,
          3612.83,
          14.9924,
          29.9387,
          5.04516,
          30960.8,
          35.6649,
          3664.05,
          17094.3,
          13.5061,
          13986.2,
          248.781,
          918.314,
          3502.82,
          64966.4,
          860.103,
          3426.51,
          515.013,
          32569.1,
          609.002,
          427.186,
          677.121,
          25.5264,
          199.22,
          8.10828"    
    
    z <- "R1_101,
          R1_102,
          R1_11,
          R1_12,
          R1_13,
          R1_14,
          R1_21,
          R1_22,
          R1_23,
          R1_31,
          R1_32,
          R1_33,
          R1_41,
          R1_42,
          R1_43,
          R1_51,
          R1_52,
          R1_53,
          R1_61,
          R1_62,
          R1_63,
          R1_71,
          R1_72,
          R1_73,
          R1_81,
          R1_82,
          R1_83,
          R1_91,
          R1_92,
          R1_93,
          R2_1,
          R2_105,
          R2_115,
          R2_116,
          R2_117,
          R2_140,
          R2_143,
          R2_151,
          R2_152,
          R2_153,
          R2_154,
          R2_18,
          R2_19,
          R2_20,
          R2_24,
          R2_26,
          R2_27,
          R2_28,
          R2_32,
          R2_33,
          R2_37,
          R2_38,
          R2_41,
          R2_42,
          R2_45,
          R2_46,
          R2_47,
          R2_53,
          R2_55,
          R2_57,
          R2_59,
          R2_6,
          R2_60,
          R2_63,
          R2_65,
          R2_66,
          R2_67,
          R2_68,
          R2_71,
          R2_72,
          R2_73,
          R2_76"

    r <- plotLinear(unlist(strsplit(z, split=',')),
                    log2(as.numeric(unlist(strsplit(x, split=',')))),
                    log2(as.numeric(unlist(strsplit(y, split=',')))),
                    showLOQ=TRUE)
    
    checkEquals(NULL, r$LOQ)
}

test.PlotLinear_4 <- function()
{
    x  <- "10.0708,
           5.0354,
    0.8886,
    14.2176,
    107.422,
    859.375,
    161.133,
    80.5664,
    1.7772,
    28.4352,
    5156.25,
    2578.12,
    483.398,
    3437.5,
    27500,
    53.7109,
    429.688,
    0.8886,
    14.2176,
    483.398,
    26.8555,
    214.844,
    3.5544,
    56.8704,
    0.0786781,
    0.0393391,
    454.963,
    7279.41,
    5156.25,
    2578.12,
    13.4277,
    107.422,
    966.797,
    966.797,
    0.629425,
    0.314713,
    120.85,
    0.4443,
    7.1088,
    2.5177,
    1.25885,
    227.482,
    3639.71,
    7734.38,
    7734.38,
    1.25885,
    0.629425,
    1289.06,
    644.531,
    60.4248,
    60.4248,
    3750,
    117.188,
    20.1416,
    10.0708,
    0.0143051,
    0.457764,
    120.85,
    120.85,
    30.2124,
    30.2124,
    80.5664,
    40.2832,
    0.269754,
    0.539507,
    1.07901,
    30.2124,
    30.2124,
    60.4248,
    2209.82,
    4419.64,
    8839.29,
    58.5938,
    1875,
    0.269754,
    0.539507,
    1.07901,
    15.1062,
    15.1062,
    0.134877,
    0.269754,
    0.539507,
    1933.59,
    1933.59,
    7734.38,
    7734.38,
    2578.12,
    1289.06,
    0.944138,
    20.1416,
    10.0708,
    1.67847,
    13.4277,
    1875,
    58.5938,
    0.839233,
    6.71387,
    0.0674384,
    0.134877,
    0.269754,
    0.0674384,
    0.134877,
    0.269754,
    0.0674384,
    0.134877,
    161.133,
    80.5664,
    0.0337192,
    0.134877,
    0.0314713,
    0.0629425,
    0.125885,
    0.25177,
    0.0168596,
    0.0674384,
    80.5664,
    40.2832,
    0.0168596,
    0.0337192,
    0.0674384,
    322.266,
    161.133,
    937.5,
    30000,
    0.4443,
    7.1088,
    0.134877,
    0.269754,
    0.4443,
    7.1088,
    29.2969,
    937.5,
    3.77655,
    937.5,
    30000,
    3.77655,
    0.114441,
    3.66211,
    28.4352,
    454.963,
    0.114441,
    3.66211,
    0.0157356,
    0.0314713,
    0.0629425,
    0.125885,
    0.0572205,
    1.83105,
    29.2969,
    937.5"

    y1 <- "0.837485,
       0.121353,
    0.0379273,
    0.768669,
    1.90631,
    0.744164,
    1.0572e-06,
    3.87601,
    8.51e-08,
    0.779118,
    214.214,
    1091.49,
    54.3986,
    380.045,
    0.15194,
    2.09719,
    17.9671,
    0.0795291,
    0.478201,
    3.56595,
    2.38232,
    15.8246,
    0.394103,
    0.895861,
    8.8e-09,
    0.0712489,
    35.6487,
    614.35,
    988.259,
    178.257,
    0.532468,
    3.26565,
    44.7566,
    41.5003,
    0.0639191,
    NA,
    4.2057,
    4.25e-08,
    0.596072,
    0.232307,
    0.0675099,
    2.28788,
    423.966,
    3907.65,
    3567.12,
    0.111393,
    3.9771e-06,
    133.993,
    154.723,
    4.34212,
    0.996274,
    85.0451,
    4.61353,
    0.163091,
    0.379198,
    NA,
    0.038774,
    3.21766,
    15.0352,
    1.38579,
    1.93683,
    9.35673,
    5.451,
    5e-08,
    0.0649624,
    0.384407,
    3.77595,
    7.44824,
    16.0071,
    270.303,
    3.494e-07,
    1671.87,
    6.812,
    325.095,
    3.935e-07,
    0.0616513,
    0.0305549,
    0.149523,
    0.726862,
    0.103768,
    1.596e-06,
    0.0463446,
    199.256,
    322.444,
    1093.72,
    1688.88,
    279.584,
    63.5678,
    0.110435,
    1.93013,
    0.374583,
    0.218183,
    1.74819,
    85.852,
    1.21e-08,
    0.377331,
    1.96912,
    NA,
    NA,
    NA,
    0.130965,
    1.18e-08,
    2.91e-08,
    NA,
    NA,
    15.1834,
    9.38777,
    NA,
    0.0799122,
    2.10385e-05,
    NA,
    0.083846,
    NA,
    0.0171778,
    NA,
    3.35384,
    5.73421,
    5.8264e-06,
    6.246e-07,
    0.0200194,
    17.5116,
    11.0122,
    260.919,
    3964.2,
    0.100965,
    1.11491,
    0.0259428,
    NA,
    0.288957,
    0.4821,
    1.20828,
    193.558,
    0.536125,
    187.495,
    7064.26,
    0.264505,
    0.0221402,
    0.396147,
    3.55531,
    73.5915,
    NA,
    NA,
    9.9647e-06,
    4.63e-08,
    4.7e-09,
    0.0237112,
    0.0392847,
    0.0456551,
    1.70831,
    91.1479"
    
    y2 <- "0.722778,
0.733873,
    5.219e-07,
    0.604539,
    1.73717,
    1.1534,
    1.9165e-06,
    3.91994,
    4.754e-07,
    0.898643,
    207.593,
    1121.36,
    54.4085,
    373.338,
    0.169796,
    2.09949,
    16.8435,
    0.0616389,
    0.619896,
    3.89005,
    2.5968,
    16.6549,
    0.182859,
    0.948602,
    1.016e-07,
    0.0403375,
    39.1606,
    609.005,
    1008.73,
    157.116,
    0.234875,
    3.39922,
    42.6545,
    40.196,
    0.181165,
    2.847e-07,
    4.69883,
    NA,
    0.189675,
    0.254871,
    0.391466,
    2.41451,
    426.951,
    3810.71,
    3577.81,
    0.0420185,
    4.103e-06,
    120.872,
    150.976,
    4.00902,
    1.09308,
    88.7421,
    4.36902,
    0.325112,
    0.487459,
    0.0794166,
    7.003e-07,
    3.75052,
    13.7295,
    1.36716,
    2.148,
    8.88517,
    5.73348,
    0.00058588,
    0.351373,
    0.329472,
    4.29572,
    7.05352,
    18.2358,
    266.712,
    3.379e-07,
    1650,
    6.13767,
    312.299,
    NA,
    0.0414129,
    NA,
    0.158325,
    1.03642,
    0.0698399,
    2.7166e-06,
    0.0904658,
    195.812,
    309.716,
    1161.1,
    1716.5,
    275.204,
    61.8612,
    NA,
    1.05298,
    0.427896,
    3.555e-07,
    2.10418,
    80.0865,
    1.96e-08,
    0.512922,
    1.58903,
    6.00372e-05,
    0.0432046,
    NA,
    2.8e-06,
    1e-10,
    0.0702308,
    NA,
    NA,
    15.5818,
    8.15287,
    NA,
    0.089151,
    0.0283771,
    1.102e-07,
    1.27e-08,
    3.0816e-06,
    NA,
    0.0412629,
    3.26266,
    7.27648,
    NA,
    NA,
    NA,
    17.6655,
    10.9249,
    249.933,
    3641.64,
    0.0509423,
    0.640191,
    0.0444007,
    9.7502e-06,
    0.319625,
    0.436569,
    1.04639,
    191.152,
    0.279646,
    167.877,
    7014.58,
    0.392405,
    2.306e-07,
    0.334322,
    2.50338,
    65.9721,
    0.0358907,
    0.350269,
    3.09e-08,
    2.8e-09,
    0.0405143,
    9.942e-07,
    NA,
    0.010433,
    1.02221,
    95.2976"
    
    y3 <- "0.490311,
0.469879,
    0.16225,
    0.490533,
    2.13129,
    1.0798,
    2.2654e-06,
    4.24638,
    0.110956,
    0.62222,
    262.925,
    1096.04,
    48.0564,
    399.27,
    0.14716,
    3.38397,
    18.2888,
    3.547e-07,
    0.761565,
    4.13325,
    1.7191,
    17.8298,
    0.088615,
    1.34844,
    NA,
    NA,
    43.9125,
    605.221,
    1007.06,
    177.07,
    0.554892,
    3.52786,
    40.612,
    46.9536,
    2.0762e-06,
    0.0886553,
    4.04549,
    0.060803,
    0.413782,
    0.11363,
    NA,
    2.72187,
    429.049,
    3861.91,
    3759.93,
    0.108981,
    3.9713e-06,
    121.845,
    146.047,
    4.52467,
    1.6103,
    93.0103,
    5.10038,
    0.51797,
    0.195218,
    2.26379e-05,
    0.0758392,
    3.94007,
    11.767,
    1.5896,
    2.04017,
    9.66727,
    6.0005,
    0.000148475,
    0.165744,
    0.127973,
    3.72053,
    6.65285,
    16.6536,
    279.149,
    2.128e-07,
    1707.2,
    7.94687,
    298.357,
    NA,
    0.0357967,
    NA,
    0.156696,
    0.911175,
    NA,
    0.0614663,
    0.0476663,
    197.844,
    317.099,
    1223.97,
    1618.46,
    277.503,
    66.3776,
    0.072033,
    1.72695,
    0.216182,
    0.172833,
    1.6133,
    89.5909,
    2.87e-08,
    0.496378,
    1.83563,
    1.506e-07,
    NA,
    0.137899,
    NA,
    NA,
    NA,
    0.0146193,
    0.0325043,
    14.8062,
    7.7748,
    0.0247236,
    NA,
    4.48882e-05,
    NA,
    0.0305807,
    0.0408539,
    0.0168537,
    NA,
    3.43573,
    7.86324,
    NA,
    NA,
    0.0196412,
    17.5211,
    11.9901,
    253.366,
    3609.29,
    9.812e-07,
    0.907263,
    0.0254687,
    NA,
    0.208575,
    0.40283,
    0.836432,
    191.248,
    0.483889,
    197.458,
    6854.51,
    0.171308,
    4.1e-08,
    0.303812,
    1.92317,
    67.3608,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    0.0269071,
    2.66745,
    97.6445"
    
    z  <- "R1_101_1,
R1_101_2,
    R1_102_1,
    R1_102_2,
    R1_103_1,
    R1_103_2,
    R1_11_1,
    R1_11_2,
    R1_12_1,
    R1_12_2,
    R1_13_1,
    R1_13_2,
    R1_14_1,
    R1_21_1,
    R1_21_2,
    R1_22_1,
    R1_22_2,
    R1_23_1,
    R1_23_2,
    R1_24_1,
    R1_31_1,
    R1_31_2,
    R1_32_1,
    R1_32_2,
    R1_33_1,
    R1_33_2,
    R1_41_1,
    R1_41_2,
    R1_42_1,
    R1_42_2,
    R1_43_1,
    R1_43_2,
    R1_51_1,
    R1_51_2,
    R1_52_1,
    R1_52_2,
    R1_53_1,
    R1_61_1,
    R1_61_2,
    R1_62_1,
    R1_62_2,
    R1_63_1,
    R1_63_2,
    R1_71_1,
    R1_71_2,
    R1_72_1,
    R1_72_2,
    R1_73_1,
    R1_73_2,
    R1_81_1,
    R1_81_2,
    R1_82_1,
    R1_82_2,
    R1_83_1,
    R1_83_2,
    R1_91_1,
    R1_91_2,
    R1_92_1,
    R1_92_2,
    R1_93_1,
    R1_93_2,
    R2_115_1,
    R2_115_2,
    R2_116_1,
    R2_116_2,
    R2_116_3,
    R2_117_1,
    R2_117_3,
    R2_140_1,
    R2_14_1,
    R2_14_2,
    R2_14_3,
    R2_150_1,
    R2_150_2,
    R2_151_1,
    R2_151_2,
    R2_151_3,
    R2_152_1,
    R2_152_2,
    R2_153_1,
    R2_153_2,
    R2_153_3,
    R2_154_1,
    R2_154_2,
    R2_18_1,
    R2_18_2,
    R2_19_1,
    R2_19_2,
    R2_1_1,
    R2_20_1,
    R2_20_2,
    R2_24_1,
    R2_24_2,
    R2_26_1,
    R2_26_2,
    R2_27_1,
    R2_27_2,
    R2_28_1,
    R2_28_2,
    R2_28_3,
    R2_32_1,
    R2_32_2,
    R2_32_3,
    R2_37_2,
    R2_37_3,
    R2_41_1,
    R2_41_2,
    R2_42_1,
    R2_42_3,
    R2_45_1,
    R2_45_2,
    R2_45_3,
    R2_45_4,
    R2_46_1,
    R2_46_3,
    R2_47_1,
    R2_47_2,
    R2_53_1,
    R2_53_2,
    R2_53_3,
    R2_54_1,
    R2_54_2,
    R2_55_2,
    R2_55_3,
    R2_57_1,
    R2_57_2,
    R2_59_2,
    R2_59_3,
    R2_60_1,
    R2_60_2,
    R2_63_1,
    R2_63_3,
    R2_65_1,
    R2_66_1,
    R2_66_2,
    R2_67_1,
    R2_68_1,
    R2_68_2,
    R2_6_1,
    R2_6_2,
    R2_71_1,
    R2_71_2,
    R2_72_1,
    R2_72_2,
    R2_72_3,
    R2_72_4,
    R2_73_1,
    R2_73_2,
    R2_7_1,
    R2_7_2"

    x  <- log2(suppressWarnings(as.numeric(unlist(strsplit(x,  split=',')))))
    y1 <- log2(suppressWarnings(as.numeric(unlist(strsplit(y1, split=',')))))
    y2 <- log2(suppressWarnings(as.numeric(unlist(strsplit(y2, split=',')))))
    y3 <- log2(suppressWarnings(as.numeric(unlist(strsplit(y3, split=',')))))
    z  <- unlist(strsplit(z, split=','))

    r <- plotLinear(c(z),
                    x,
                    data.frame(y1=y1, y2=y2, y3=y3),
                    showLOQ=TRUE, showSD=TRUE)
    
    checkEqualsNumeric(0.8296060465, r$LOQ$model$breaks)
}

test.PlotLinear_1()
test.PlotLinear_2()
test.PlotLinear_3()
test.PlotLinear_4()