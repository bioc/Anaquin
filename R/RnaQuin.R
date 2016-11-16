#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Garvan Institute of Medical Research
#

#
# Implements functionality related to RnaQuin
#

#
# Aggregate isoform concentration for genes
#
RnaQuin.aggregate <- function(isos, concent)
{
    data <- data.frame(ID=isos, Mix=concent)
    
    # Convert isoforms to genes
    data$ID <- RnaQuin.iso2gen(data$ID)
    
    # Aggregate isoforms into genes    
    return (aggregate(Mix~ID, data=data, FUN=sum))
}

#
# Filter RnaQuin genes from a list of genes names.
#
RnaQuin.genes <- function(genes)
{
    data('RnaQuinMixtureA')
    genes[genes %in% RnaQuin.iso2gen(RnaQuinMixtureA$ID)]
}

#
# Filter RnaQuin isoforms from a list of isoform names.
#
RnaQuin.isoforms <- function(isos)
{
    data('RnaQuinMixtureA')
    isos[isos %in% RnaQuinMixtureA$ID]
}

#
# Convert RnaQuin sequin isoforms to sequin genes.
#
#    Eg: R1_1_1 to R1_1
#
RnaQuin.iso2gen <- function(isos)
{
    unlist(lapply(strsplit(as.character(isos), '_'),  function(x)
    {
        paste(x[1:length(x)-1], collapse='_')
    }))
}