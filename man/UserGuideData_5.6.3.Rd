\docType{data}
\name{UserGuideData_5.6.3}
\alias{UserGuideData_5.6.3}
\title{Differential expression (RnaQuin)}
\format{Data frame:
\itemize{
  \item ExpLFC: Expected log-fold change 
  \item ObsLFC: Observed log-fold change
  \item SD: Standard deviation of the measurment
  \item Pval: P-value probability
  \item Qval: Q-value probability
  \item Mean: Average counts across the samples
  \item Label: Average counts across the samples  
}}
\usage{
data(UserGuideData_5.6.3)
}
\description{
Differential gene expression estimated by DESeq2. \code{Section 5.6.3}
has details on the data set.
}
\value{
Data frame with columns defined in \code{Format}.
}
\source{
S.A Hardwick. Spliced synthetic genes as internal controls in RNA sequencing
experiments. \code{Nature Methods}, 2016.
}