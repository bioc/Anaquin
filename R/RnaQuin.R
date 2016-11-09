#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Garvan Institute of Medical Research
#

#
# Filter RnaQuin genes from a list of genes names.
#
RnaQuin.genes <- function(genes)
{
    data('MixtureA')
    genes[genes %in% RnaQuin.iso2gen(row.names(MixtureA))]
}

#
# Convert RnaQuin sequin isoforms to sequin genes.
#
#    Eg: R1_1_1 to R1_1
#
RnaQuin.iso2gen <- function(names)
{
    names <- strsplit(as.character(names), '_')
    
    f <- function(x)
    {
        paste(x[1:length(x)-1], collapse='_')
    }
    
    unlist(lapply(names, f))
}