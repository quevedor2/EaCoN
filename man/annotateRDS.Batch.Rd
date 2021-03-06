% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output_builder.R
\name{annotateRDS.Batch}
\alias{annotateRDS.Batch}
\title{Batch wrapper for annotateRDS()}
\usage{
annotateRDS.Batch(
  all.fits,
  segmenter,
  nthread = 1,
  cluster.type = "PSOCK",
  bin.size = 50000,
  ...
)
}
\arguments{
\item{all.fits}{Sample named-list of all samples containing two elements:
fit[data.frame]: data.frame of the *.gammaEval.txt file
sample[character]: sample name}

\item{segmenter}{Segmenter used (E.g. ASCAT) (Only ASCAT works currently)}

\item{nthread}{Max number of threads to use [Default=1]}

\item{cluster.type}{Cluster type [Default=PSOCK]}

\item{...}{}
}
\value{
Returns annotateRDS() objects
}
\description{
Batch wrapper for annotateRDS()
}
\examples{
    gr.cnv <- annotateRDS.Batch(all.fits, toupper(segmenter), nthread=3)
}
