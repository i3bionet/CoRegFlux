#' SC_GRN_1 data
#'
#' A coregnet object infered from the m3d dataset describing the gene
#' regulatory network for S. cerevisiae as described in Banos, D. T., Trébulle,
#' P., & Elati, M. (2017). Integrating transcriptional activity in genome-scale
#' models of metabolism. BMC systems biology, 11(7), 134.
#' @format a coregnet object inferred using the package _CoRegNet_
#' \describe{
#'   \item{Number of transcription factor }{200}
#'   \item{Number of targets genes}{3748}
#'   \item{Evidences}{TRUE}
#'   ...
#' }
"SC_GRN_1"

#' aliases_SC data
#'
#' A data.frame containing the gene ID used in the metabolic model and their
#' common name, used in the gene regulatory network
#' @format a two colums data.frame which first columns correspond to the name
#' used in the model and the second to the ID used in the GRN (common name).
#' Those columns should be named geneName_model and geneName_GRN respectively.
#' \describe{
#' \item{geneName_model}{Aliases or gene names used in the gene-association
#' field in the genome-scale metabolic model}
#' \item{geneName_GRN}{Aliases or gene names used in the gene regulatory network}}
"aliases_SC"

#' SC_EXP_DATA data
#'
#' A matrix of S. cerevisiae gene expression in various experimental designs,
#' derived from the m3d dataset to infer _S. cerevisiae_ gene regulatory network.
#' The dataset was shorten to 3600 genes in order to limit the size of the
#' object
#' @format a matrix of 3600 genes by 247 samples
#' @source subset of m3d dataset available at <http://m3d.mssm.edu/>
"SC_EXP_DATA"

#' SC_Test_data data
#'
#' A matrix of S. cerevisiae gene expression during diauxic shift (Brauer and
#' al.)
#' @format a matrix of 6028 genes by 13 samples during diauxic shift
#' @source E-GEOD-4398 (Brauer MJ and al.)
"SC_Test_data"

#' SC_experiment_influence data
#'
#' A vector of influence computed from the first sample of SC_Test_data
#' @format a named numerical vector
"SC_experiment_influence"

#' iMM904
#'
#' A _S. cerevisiae_ genome-scale metabolic model as a modelOrg object
#' @format a modelOrg object as required by _sybil_. See _sybilSBML_ for more
#' information on how to load other model.
"iMM904"

#' PredictedGeneState data
#'
#' Predicted gene states as obtained by the function
#' \code{predict_linear_model_influence}
#' @format a named vector containing the gene name and its associated predicted
#' gene state.
"PredictedGeneState"

#' ODcurveToMetCurve data
#'
#' List as obtained by the function \code{ODCurveToMetabolicGeneCurves}
#' @format List as obtained by the function \code{ODCurveToMetabolicGeneCurves}
"ODcurveToMetCurve"

#' ODtoflux data
#'
#' List as obtained by the function \code{ODCurveToFluxCurves}
#' @format List as obtained by the function \code{ODCurveToFluxCurves}
"ODtoflux"
