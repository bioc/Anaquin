\name{AnaquinData}
\alias{AnaquinData}
\title{Create Anaquin dataset}
\description{
Create an Anaquin dataset. The resulting S4 object is required for every
Anaquin analysis.
}
\usage{
AnaquinData(analysis, ...)
}
\arguments{
  \item{analysis}{Type of analysis to be performed
                  (eg. \code{'PlotLinear'},
                       \code{'PlotLogistic'},
                       \code{'PlotROC'},
                       \code{'PlotLODR'}.)}
  \item{...}{Analysis specific data inputs.}

The function requires the following mandatory data inputs:

\tabular{ll}{
   \code{analysis} \tab Type of analysis \cr
   \code{seqs}     \tab Sequin names  \cr
   \code{...}      \tab Analysis specific data inputs  \cr   
}

\code{analysis} indicates the type of the analysis. The current release
supports the following analysis:

\tabular{ll}{
   \code{PlotLinear}   \tab Linear model with sequins \cr
   \code{PlotLogistic} \tab GLM logistic model with sequins  \cr
   \code{PlotROC}      \tab ROC analysis with sequins  \cr
   \code{PlotLODR}     \tab LODR (LOESS) analysis with sequins  \cr
}

\code{seqs} gives the sequin names. The function will give
an error message unless both \code{analysis} and \code{seqs} are given.

The function accepts the following optional data inputs:

\tabular{ll}{
   \code{std}      \tab Standard deviation \cr
   \code{pval}     \tab P-value probability  \cr
   \code{qval}     \tab Q-value probability  \cr
   \code{ratio}    \tab Expected sequin ratio \cr
   \code{input}    \tab Input concentration (attomol/ul)  \cr
   \code{measured} \tab Measured variable (eg: FPKM)  \cr
   \code{label}    \tab Classified labels (eg: 'TP', 'FP')  \cr
   \code{score}    \tab Value used for ranking sequins \cr
}  
}
\details{
Create an Anaquin dataset. The resulting R-object is required for every Anaquin analysis.
}
\value{
An S4 object of class \code{AnaquinData}.
}
\examples{
#
# Example 1: Create Anaquin dataset for PlotLogistic
#

data(UserGuideData_5.4.5.1)

title <- 'Assembly Plot'
xlab  <- 'Input Concentration (log2)'
ylab  <- 'Sensitivity'

# Sequin names
seqs <- row.names(UserGuideData_5.4.5.1)

# Input concentration
input <- log2(UserGuideData_5.4.5.1$InputConcent)

# Measured sensitivity
measured <- UserGuideData_5.4.5.1$Sn

anaquin <- AnaquinData(analysis='PlotLogistic',
                           seqs=seqs,
                          input=input,
                       measured=measured)

plotLogistic(anaquin, title=title, xlab=xlab, ylab=ylab, showLOA=TRUE)

#
# Example 2: Create Anaquin dataset for PlotLinear
#

data(UserGuideData_5.4.6.3)

title <- 'Gene Expression'
xlab  <- 'Input Concentration (log2)'
ylab  <- 'FPKM (log2)'

# Sequin names
seqs <- row.names(UserGuideData_5.4.6.3)

# Input concentration
input <- log2(UserGuideData_5.4.6.3$InputConcent)

# Measured FPKKM
measured <- log2(UserGuideData_5.4.6.3[,2:4])

anaquin <- AnaquinData(analysis='PlotLinear',
                           seqs=seqs,
                          input=input,
                       measured=measured)

plotLinear(anaquin, title=title, xlab=xlab, ylab=ylab, showLOQ=TRUE)

#
# Example 3: Create Anaquin dataset for plotROC
#

data(UserGuideData_5.6.3)

# Sequin names
seqs <- row.names(UserGuideData_5.6.3)

# Expected log-fold
ratio <- UserGuideData_5.6.3$ExpLFC

# How ROC points are ranked (scoring function)
score <- 1-UserGuideData_5.6.3$Pval

# Classified labels (TP/FP)
label <- UserGuideData_5.6.3$Label

anaquin <- AnaquinData(analysis='PlotROC',
                           seqs=seqs,
                          ratio=ratio,
                          score=score,
                          label=label)

plotROC(anaquin, title='ROC Plot', refRats=0)

#
# Example 4: Create Anaquin dataset for plotLODR
#

data(UserGuideData_5.6.3)

xlab  <- 'Average Counts'
ylab  <- 'P-value'
title <- 'LODR Curves'

# Sequin names
seqs <- row.names(UserGuideData_5.6.3)

# Measured average mean
measured <- UserGuideData_5.6.3$Mean

# Expected log-fold
ratio <- UserGuideData_5.6.3$ExpLFC

# P-value
pval <- UserGuideData_5.6.3$Pval

# Q-value
qval <- UserGuideData_5.6.3$Qval

anaquin <- AnaquinData(analysis='PlotLODR',
                           seqs=seqs,
                       measured=measured,
                          ratio=ratio,
                           pval=pval,
                           qval=qval)

plotLODR(anaquin, xlab=xlab, ylab=ylab, title=title, FDR=0.1)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}