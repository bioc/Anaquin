#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Garvan Institute pf Medical Research
#

.RnaQuin.isoforms <- function(object)
{
    return (data.frame(row.names=seqs(object), Mix=input(object)))
}

.RnaQuin.genes <- function(object)
{
    data <- data.frame(ID=seqs(object), Mix=input(object))

    # Convert isoforms to genes
    data$ID <- RnaQuin.gen2iso(data$ID)
    
    # Aggregate isoforms into genes    
    return (aggregate(Mix~ID, data=data, FUN=sum))
}

#
# Convert RnaQuin sequin isoforms to sequin genes. The inputs assumed be valid.
#
#    Eg: R1_1_1 to R1_1
#
RnaQuin.gen2iso <- function(names)
{
    names <- as.character(names)
    names <- strsplit(names, '_')
    
    f <- function(x)
    {
        paste(x[1:length(x)-1], collapse='_')
    }
    
    unlist(lapply(names, f))
}