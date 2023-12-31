#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library(RUnit)
library(Anaquin)

test.PlotLogistic_1 <- function()
{
    data(UserGuideData_5.4.5.1)
    x <- UserGuideData_5.4.5.1

    r <- plotLogistic(row.names(x),
                      log2(x$Input),
                      x$Sn)

    # Expected fitted values by the sigmoid function
    fitted <- c(0.505755856395065,
                0.222643528893781,
                0.010300016219899,
                0.648997187338457,
                0.92026304251071,
                0.928213869071803,
                0.924632245089082,
                0.914366418394931,
                0.0378973207098846,
                0.833752363429975,
                0.928358782240202,
                0.928345459227316,
                0.927911002760357,
                0.928353148175993,
                0.928363361215984,
                0.89836631575432,
                0.927796005796745,
                0.010300016219899,
                0.648997187338457,
                0.927911002760357,
                0.823993501178706,
                0.926214306703906,
                0.129043885951807,
                0.901398794826717,
                9.82865701336804e-05,
                2.59121362851877e-05,
                0.927855061830411,
                0.928361094753171,
                0.928358782240202,
                0.928345459227316,
                0.627066585253743,
                0.92026304251071,
                0.928244209392258,
                0.928244209392258,
                0.00533476431533362,
                0.00141231744552197,
                0.921893786460962,
                0.00273763028622125,
                0.352602355170263,
                0.0712810947952379,
                0.0199170597244443,
                0.926437613006538,
                0.928354231359681,
                0.928361365200111,
                0.928361365200111,
                0.0199170597244443,
                0.00533476431533362,
                0.928294923575079,
                0.928103273636578,
                0.904289792326221,
                0.904289792326221,
                0.928354751455318,
                0.921502218658016,
                0.760780539182587,
                0.505755856395065,
                3.70222949553412e-06,
                0.00289892966464432,
                0.921893786460962,
                0.921893786460962,
                0.843210863415579,
                0.843210863415579,
                0.0115580895205083,
                0.914366418394931,
                0.877413012253192,
                0.00105034244261124,
                0.00397177716642234,
                0.0148884100362838,
                0.843210863415579,
                0.843210863415579,
                0.904289792326221,
                0.0115580895205083,
                0.928339214782058,
                0.92835713602241,
                0.928361860492618,
                0.902862250061517,
                0.928330169636822,
                0.00105034244261124,
                0.00397177716642234,
                0.0148884100362838,
                0.671229716373391,
                0.671229716373391,
                0.000277119451357618,
                0.00105034244261124,
                0.00397177716642234,
                0.928332087966142,
                0.928332087966142,
                0.928361365200111,
                0.928361365200111,
                0.928345459227316,
                0.928294923575079,
                0.0115580895205083,
                0.760780539182587,
                0.505755856395065,
                0.0340966826710067,
                0.627066585253743,
                0.928330169636822,
                0.902862250061517,
                0.00923829474485884,
                0.328897268974982,
                7.30693701509513e-05,
                0.000277119451357618,
                0.00105034244261124,
                7.30693701509513e-05,
                0.000277119451357618,
                0.00105034244261124,
                5.65192799587101e-05,
                1.92634758951033e-05,
                7.30693701509513e-05,
                0.000277119451357618,
                3.09052920252929e-07,
                1.17234994662851e-06,
                4.44713748735444e-06,
                1.68695415881578e-05,
                0.924632245089082,
                0.914366418394931,
                1.92634758951033e-05,
                7.30693701509513e-05,
                0.000277119451357618,
                1.68695415881578e-05,
                6.39890691148908e-05,
                0.000242687866644135,
                0.000919934543895209,
                5.0782655595262e-06,
                1.92634758951033e-05,
                7.30693701509513e-05,
                0.914366418394931,
                0.877413012253192,
                5.0782655595262e-06,
                1.92634758951033e-05,
                7.30693701509513e-05,
                0.927376994426196,
                0.924632245089082,
                0.928236933392859,
                0.928363390588447,
                0.00273763028622125,
                0.352602355170263,
                7.30693701509513e-05,
                0.000277119451357618,
                0.00105034244261124,
                0.00273763028622125,
                0.352602355170263,
                0.838521261207461,
                0.928236933392859,
                0.142553466244861,
                0.928236933392859,
                0.928363390588447,
                0.142553466244861,
                0.000202045715636617,
                0.135557087471656,
                0.833752363429975,
                0.927855061830411,
                0.000202045715636617,
                0.135557087471656,
                4.44713748735444e-06,
                1.68695415881578e-05,
                6.39890691148908e-05,
                0.000242687866644135,
                5.32713092910686e-05,
                0.0400403521055245,
                1.33872469655316e-06,
                5.0782655595262e-06,
                1.92634758951033e-05,
                0.838521261207461,
                0.928236933392859)
    
    checkEquals(r$LOA, 3.33)
    #all(checkEqualsNumeric(fitted, r$fitted))
}

test.PlotLogistic_2 <- function()
{
    data(UserGuideData_5.4.6.3)
    x <- UserGuideData_5.4.6.3

    checkException(plotLogistic(row.names(data),
                                log2(data$Input),
                                log2(data$Observed1)))
}

test.PlotLogistic_1()
test.PlotLogistic_2()