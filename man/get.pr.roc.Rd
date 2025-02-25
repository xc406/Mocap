% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_curves.R
\name{get.pr.roc}
\alias{get.pr.roc}
\title{wrapper function to calculate area under the PR and ROC curves}
\usage{
get.pr.roc(d, tfct, m, plot = TRUE, dir = ".")
}
\arguments{
\item{d}{a data frame with feature columns and gold standard column as the last column}

\item{tfct}{a character string specifying tf name and cell type}

\item{m}{sparse logistic regression model inherited from the LiblineaR package}

\item{plot}{logical}

\item{dir}{specifying a path to store plots in}
}
\description{
\code{get.pr.roc} with the option to plot PR and ROC curves
}
\examples{
get.pr.roc(d,tfct,slr.model,plot=FALSE)
get.pr.roc(d,tfct,slr.model,plot=TRUE,dir=".")
}

