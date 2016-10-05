#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Garvan Institute pf Medical Research
#

#
# Anaquin is a flexible framework for Bioinformatics, thus the data
# representation vary by analysis. The class representation lists all
# the possible slots, however, only some of them will be needed for an
# analysis. It's important to check R-man and vigenette on what inputs
# are required.
#
#   Supported analysis:
#
#     - Mixture 
#     - PlotROC
#     - PlotLODR
#     - PlotLinear
#     - PlotLogistic
#

.validate = function(object)
{
    aly <- analysis(object)

    if (is.null(aly))
    {
        return ('Type of analyis is required')
    }

    # Expected analysis types
    expAly = c('PlotLODR', 'PlotROC', 'PlotLinear', 'PlotLogistic', 'Mixture',
               'plotLODR', 'plotROC', 'plotLinear', 'plotLogistic', 'mixture')

    if (!(aly %in% expAly))
    {
        choices <- paste(expAly, collapse=', ')
        str <- c('Invalid analysis: ', aly, '. Should be {', choices, '}.')
        return (paste(str, sep = ""))
    }
    
    if (aly == 'PlotLODR' | aly == 'plotLODR')
    {
        if (is.null(ratio(object)))
        {
            return ('PlotLODR requires expected ratio.
                     Please specifiy it with "ratio"')
        }
      
        if (is.null(pval(object)))
        { 
            return ('PlotLODR requires p-value probability.
                     Please specifiy it with "pval"')
        }
      
        if (is.null(measured(object)))
        { 
            return ('PlotLODR requires measurement.
                     Please specifiy it with "measured"')
        }
    }
    
    if (aly == 'PlotROC' | aly == 'plotROC')
    {
        if (is.null(ratio(object)))
        {
            return ('PlotROC requires expected ratio.
                     Please specifiy it with "ratio"')
        }

        if (is.null(label(object)))
        {
            return ('PlotROC requires label (TP/FP).
                     Please specifiy it with "label"')
        }
        
        if (all(is.na(score(object))))
        {
            return ('PlotROC requires scoring for ranking.
                     Please specifiy it with "score"') }
    }

    if (aly == 'PlotLinear' | aly == 'plotLinear')
    {
        if (is.null(input(object)))
        {
            return ('PlotLinear requires input concentration.
                     Please specifiy it with "input"')
        }
        
        if (is.null(measured(object)))
        {
            return ('PlotLinear requires measurement (eg: FPKM).
                     Please specifiy it with "measured"')
        }
    }

    if (aly == 'PlotLogistic' | aly == 'plotLogistic')
    {
        if (is.null(input(object)))
        {
            return ('PlotLogistic requires input concentration.
                     Please specifiy it with "input"')
        }
        
        if (is.null(measured(object)))
        {
            return ('PlotLogistic requires measurement (eg: FPKM).
                     Please specifiy it with "measured"')
        }
      
        if (!is.vector(measured(object)))
        {
            return ('PlotLogistic requires a vector of measurement.')
        }
    }

    return (TRUE)
}

setClassUnion("factorOrNULL",members=c("factor", "character", "NULL")) 
setClassUnion("numericOrNULL",members=c("numeric", "NULL")) 
setClassUnion("characterOrNULL",members=c("character", "NULL"))
setClassUnion("data.frameORvectorOrNULL", c("data.frame", "vector", "NULL"))

#' An S4 class to represent an Anaquin data set.
#'
#' @slot analysis Analysis type
#' @slot seqs     Sequin names
#' @slot std      Standard deviation
#' @slot pval     P-value probability
#' @slot qval     Q-value probability
#' @slot ratio    Sequin ratio
#' @slot input    Input concentration (attomol/ul)
#' @slot measured Measured abundance (eg: FPKM)
#' @slot label    Classified label (eg: 'TP'/'FP')
#' @slot score    Scoring for ROC

setClass("AnaquinData", representation(analysis = 'character',
                                       seqs     = 'factorOrNULL',
                                       std      = 'numericOrNULL',
                                       pval     = 'numericOrNULL',
                                       qval     = 'numericOrNULL',
                                       ratio    = 'numericOrNULL',
                                       input    = 'numericOrNULL',
                                       measured = 'data.frameORvectorOrNULL',
                                       label    = 'factorOrNULL',
                                       score    = 'numericOrNULL'),
                          prototype(std         = NULL,
                                    pval        = NULL,
                                    qval        = NULL,
                                    ratio       = NULL,
                                    input       = NULL,
                                    measured    = NULL,
                                    label       = NULL,
                                    score       = NULL),
                        validity = .validate)

AnaquinData <- function(analysis, ...)
{
    new("AnaquinData", analysis=analysis, ...)
}

setMethod("show",
           signature = "AnaquinData",
           definition = function(object)
           {
               aly <- analysis(object)

               cat('\nAnaquin data for ', aly, ' analysis', '\n', sep='')
               cat('Number of sequins: ', length(seqs(object)), '\n\n', sep='')
               cat('- Use seqs(x) to get sequins \n', sep='')
               
               if (aly == 'PlotLODR')
               {
                     cat('- Use ratio(x) to get expected ratios \n', sep='')
                     cat('- Use score(x) to get measurement \n', sep='')
                     cat('- Use label(x) to get p-value \n', sep='')
               }

               if (aly == 'PlotROC')
               {
                     cat('- Use ratio(x) to get expected ratios \n', sep='')
                     cat('- Use score(x) to get scores \n', sep='')
                     cat('- Use label(x) to get labels \n', sep='')
               }
               
               if (aly == 'PlotLinear')
               {
                   cat('- Use input(x) to get sequins \n', sep='')
                   cat('- Use measured(x) to get measurements \n', sep='')
               }

               if (aly == 'PlotLogistic')
               {
                   cat('- Use input(x) to get input concentration \n', sep='')
                   cat('- Use measured(x) to get measurements \n', sep='')
               }
               
               cat('\nwhere x is your object name\n\n', sep='')
               invisible(NULL)
           })

setGeneric('std',      function(object, ...) standardGeneric('std'))
setGeneric('pval',     function(object, ...) standardGeneric('pval'))
setGeneric('qval',     function(object, ...) standardGeneric('qval'))
setGeneric('ratio',    function(object, ...) standardGeneric('ratio'))
setGeneric('input',    function(object, ...) standardGeneric('input'))
setGeneric('label',    function(object, ...) standardGeneric('label'))
setGeneric('score',    function(object, ...) standardGeneric('score'))
setGeneric('seqs',     function(object, ...) standardGeneric('seqs'))
setGeneric('measured', function(object, ...) standardGeneric('measured'))
setGeneric('analysis', function(object, ...) standardGeneric('analysis'))
setGeneric('RnaQuin.genes',
                       function(object, ...) standardGeneric('RnaQuin.genes'))
setGeneric('RnaQuin.isoforms',
                       function(object, ...) standardGeneric('RnaQuin.isoforms'))

setMethod('std',      'AnaquinData', function(object) object@std)
setMethod('pval',     'AnaquinData', function(object) object@pval)
setMethod('qval',     'AnaquinData', function(object) object@qval)
setMethod('ratio',    'AnaquinData', function(object) object@ratio)
setMethod('input',    'AnaquinData', function(object) object@input)
setMethod('label',    'AnaquinData', function(object) object@label)
setMethod('score',    'AnaquinData', function(object) object@score)
setMethod('seqs',     'AnaquinData', function(object) object@seqs)
setMethod('analysis', 'AnaquinData', function(object) object@analysis)
setMethod('measured', 'AnaquinData', function(object) object@measured)

setMethod('RnaQuin.genes',
                      'AnaquinData', function(object) .RnaQuin.genes(object))
setMethod('RnaQuin.isoforms',
                      'AnaquinData', function(object) .RnaQuin.isoforms(object))

.m2str <- function(m)
{
    eq <- substitute(italic(y) == a + b * italic(x)*','~~italic(r)^2~'='~r2, 
                     list(a  = format(coef(m)[1], digits = 2), 
                          b  = format(coef(m)[2], digits = 2), 
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}

.lm2str <- function(data)
{
    return (.m2str(lm(y~x, data)))
}