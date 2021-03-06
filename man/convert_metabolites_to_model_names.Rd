% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MetabolicModelFunctions.R
\name{convert_metabolites_to_model_names}
\alias{convert_metabolites_to_model_names}
\title{Convert metabolites name to the model equivalent}
\usage{
convert_metabolites_to_model_names(metabolites, model,
  exchange_met = build_exchange_met(model))
}
\arguments{
\item{metabolites}{A data.frame containing the metabolites and their
concentrations}

\item{model}{A genome-scale metabolic model as a modelorg object}

\item{exchange_met}{A data.frame as build by the function build_exchange_met}
}
\value{
A data.frame containing the exchange metabolite model id and the
equivalent name
}
\description{
Convert metabolites name to the model equivalent
}
\keyword{internal}
