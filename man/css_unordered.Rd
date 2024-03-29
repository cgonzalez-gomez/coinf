% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CSS.R
\name{css_unordered}
\alias{css_unordered}
\title{CSS score for unordered disease query}
\usage{
css_unordered(R, S)
}
\arguments{
\item{R}{named vector with drug perturbation data representing the drug gene
list of interest with expression information.}

\item{S}{unordered named vector with the sign of the disease expression
of each gene (-1 or +1).}
}
\value{
CSS score
(https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-258)
}
\description{
CSS score for unordered disease query
}
