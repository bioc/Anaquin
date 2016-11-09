#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Garvan Institute of Medical Research
#

#
# Anaquin is a flexible framework for bioinformatics, thus the data
# representation vary by analysis. The class representation shows all
# possible slots, however, only some of them will be needed for any specific
# analysis. It's important to check the documentation on what inputs are
# required.
#
#   - PlotMA
#   - PlotROC
#   - PlotLODR
#   - PlotLinear
#   - PlotLogistic
#

# Whether the classified labels are expected
.verifyLabel <- function(x)
{
    return (all(levels(x) %in% c("FP", "TP"))) 
}

.validate <- function(object)
{
    aly <- analysis(object)

    if (is.null(aly))
    {
        return ('Type of analyis is required')
    }

    # Expected analysis types
    expAly = c('PlotLODR', 'PlotROC', 'PlotLinear', 'PlotLogistic', 'plotMA',
               'plotLODR', 'plotROC', 'plotLinear', 'plotLogistic', 'PlotMA')

    if (!(aly %in% expAly))
    {
        choices <- paste(expAly, collapse=', ')
        str <- c('Invalid analysis: ', aly, '. Should be {', choices, '}.')
        return (paste(str, sep = ""))
    }
    
    if (aly == 'PlotMA' | aly == 'plotMA')
    {
        if (is.null(mean(object)))
        {
            return ('PlotMA requires replicate means (x-axis).
                    Please specifiy it with "mean"')
        }
        
        if (is.null(lfc(object)))
        { 
            return ('PlotMA requires log-fold change (y-axis).
                    Please specifiy it with "lfc"')
        }
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
                     Please specifiy it with "score"')
        }
        
        if (!.verifyLabel(label(object)))
        {
            return ('PlotROC requires input labels are only \'TP\'
                     or \'FP\'')
        }
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
#' @slot ratio    Sequin ratio
#' @slot input    Sequin concentration (attomol/ul)
#' @slot std      Standard deviation
#' @slot pval     P-value probability
#' @slot qval     Q-value probability
#' @slot measured Measured abundance (eg: FPKM)
#' @slot label    Classified label ('TP'/'FP')
#' @slot score    Scoring for ROC
#' @slot mean     Mean (eg: baseMean in DESeq2)
#' @slot lfc      Log-fold change (eg: differential analysis)

setClass("AnaquinData", representation(analysis = 'character',
                                       seqs     = 'factorOrNULL',
                                       std      = 'numericOrNULL',
                                       pval     = 'numericOrNULL',
                                       qval     = 'numericOrNULL',
                                       ratio    = 'numericOrNULL',
                                       input    = 'numericOrNULL',
                                       measured = 'data.frameORvectorOrNULL',
                                       label    = 'factorOrNULL',
                                       score    = 'numericOrNULL',
                                       mean     = 'numericOrNULL',
                                       lfc      = 'numericOrNULL'),
                          prototype(std         = NULL,
                                    pval        = NULL,
                                    qval        = NULL,
                                    ratio       = NULL,
                                    input       = NULL,
                                    measured    = NULL,
                                    label       = NULL,
                                    score       = NULL,
                                    mean        = NULL,
                                    lfc         = NULL),
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
               
               if (aly == 'PlotMA')
               {
                   cat('- Use mean(x) to get means \n', sep='')
                   cat('- Use lfc(x) to get log-fold ratio \n', sep='')
               }
               else if (aly == 'PlotLODR')
               {
                     cat('- Use ratio(x) to get expected ratios \n', sep='')
                     cat('- Use score(x) to get measurement \n', sep='')
                     cat('- Use label(x) to get p-value \n', sep='')
               }
               else if (aly == 'PlotROC')
               {
                     cat('- Use ratio(x) to get expected ratios \n', sep='')
                     cat('- Use score(x) to get scores \n', sep='')
                     cat('- Use label(x) to get labels \n', sep='')
               }
               else if (aly == 'PlotLinear')
               {
                   cat('- Use input(x) to get sequins \n', sep='')
                   cat('- Use measured(x) to get measurements \n', sep='')
               }
               else if (aly == 'PlotLogistic')
               {
                   cat('- Use input(x) to get input concentration \n', sep='')
                   cat('- Use measured(x) to get measurements \n', sep='')
               }
               
               cat('\nwhere x is your object name\n\n', sep='')
               invisible(NULL)
           })

setGeneric('lfc',      function(object, ...) standardGeneric('lfc'))
setGeneric('std',      function(object, ...) standardGeneric('std'))
setGeneric('pval',     function(object, ...) standardGeneric('pval'))
setGeneric('qval',     function(object, ...) standardGeneric('qval'))
setGeneric('seqs',     function(object, ...) standardGeneric('seqs'))
setGeneric('mean',     function(object, ...) standardGeneric('mean'))
setGeneric('ratio',    function(object, ...) standardGeneric('ratio'))
setGeneric('input',    function(object, ...) standardGeneric('input'))
setGeneric('label',    function(object, ...) standardGeneric('label'))
setGeneric('score',    function(object, ...) standardGeneric('score'))
setGeneric('measured', function(object, ...) standardGeneric('measured'))
setGeneric('analysis', function(object, ...) standardGeneric('analysis'))

setMethod('lfc',      'AnaquinData', function(object) object@lfc)
setMethod('std',      'AnaquinData', function(object) object@std)
setMethod('pval',     'AnaquinData', function(object) object@pval)
setMethod('qval',     'AnaquinData', function(object) object@qval)
setMethod('seqs',     'AnaquinData', function(object) object@seqs)
setMethod('mean',     'AnaquinData', function(object) object@mean)
setMethod('ratio',    'AnaquinData', function(object) object@ratio)
setMethod('input',    'AnaquinData', function(object) object@input)
setMethod('label',    'AnaquinData', function(object) object@label)
setMethod('score',    'AnaquinData', function(object) object@score)
setMethod('analysis', 'AnaquinData', function(object) object@analysis)
setMethod('measured', 'AnaquinData', function(object) object@measured)

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