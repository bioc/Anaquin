\docType{data}
\name{UserGuideData_5.4.6.3}
\alias{UserGuideData_5.4.6.3}
\title{Gene expression (RnaQuin)}
\format{Data frame:
\itemize{
  \item InputConcent: Input concentration in attomol/ul
  \item Observed1: Measured FPKM for the first replicate
  \item Observed2: Measured FPKM for the second replicate
  \item Observed3: Measured FPKM for the third replicate  
}}
\usage{
data(UserGuideData_5.4.6.3)
}
\description{
Gene expression estimated by Cufflinks. \code{Section 5.4.6.3} of the Anaquin
user guide has details on the data set.
}
\value{
Data frame with columns defined in \code{Format}.
}
\source{
S.A Hardwick. Spliced synthetic genes as internal controls in RNA sequencing
experiments. \code{Nature Methods}, 2016.
}