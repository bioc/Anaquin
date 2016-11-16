#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library(RUnit)
library(Anaquin)

test.RnaQuin.iso2gen_1 <- function()
{
    x <- c('R1_1_1', 'R1_1_2', 'R1_1_3', 'R1_2_A')
    r <- RnaQuin.iso2gen(x)

    checkEquals(4, length(r))
    checkEquals('R1_1', r[[1]])
    checkEquals('R1_1', r[[2]])
    checkEquals('R1_1', r[[3]])
    checkEquals('R1_2', r[[4]])    
}

test.RnaQuin.genes_1 <- function()
{
    x <- c('R1_1_1', 'R1_11', 'R1_1_3', 'R2_54')
    r <- RnaQuin.genes(x)
    
    checkEquals(2, length(r))
    checkEquals('R1_11', r[[1]])
    checkEquals('R2_54', r[[2]])
}

test.RnaQuin.isoforms_1 <- function()
{
    x <- c('R1_11_1', 'This_Is_Not_Sequin', 'R1_12_2')
    r <- RnaQuin.isoforms(x)
    
    checkEquals(2, length(r))
    checkEquals('R1_11_1', r[[1]])
    checkEquals('R1_12_2', r[[2]])
}

test.RnaQuin.aggregate_1 <- function()
{
    data('RnaQuinMixtureA')
    r <- RnaQuin.aggregate(RnaQuinMixtureA$ID, RnaQuinMixtureA$MixA)
    
    checkEquals(78, nrow(r))
    checkEqualsNumeric(15.10620117, r[r$ID=='R1_101',]$Mix)
    checkEqualsNumeric(30.21240235, r[r$ID=='R1_12',]$Mix)
}

test.RnaQuin.genes_1()
test.RnaQuin.iso2gen_1()
test.RnaQuin.isoforms_1()
test.RnaQuin.aggregate_1()