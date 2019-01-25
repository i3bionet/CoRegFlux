#' SC_GRN_1
#' A coregnet object infered from the m3d dataset describing the gene
#' regulatory network for S. cerevisiae as described in Banos, D. T., Trébulle,
#' P., & Elati, M. (2017). Integrating transcriptional activity in genome-scale
#' models of metabolism. BMC systems biology, 11(7), 134.
#' \describe{
#'   \item{Number of transcription factor }{200}
#'   \item{Number of targets genes}{3748}
#'   \item{Evidences}{TRUE}
#'   ...
#' }
"SC_GRN_1"

#' aliases_SC
#' A data.frame containing the gene ID used in the metabolic model and their
#' common name, used in the gene regulatory network
#' @format a two colums data.frame which first columns correspond to the name
#' used in the model and the second to the ID used in the GRN (common name).
#' Those columns should be named geneName_model and geneName_GRN respectively.
"aliases_SC"

#' SC_EXP_DATA
#' A matrix of S. cerevisiae gene expression in various experimental designs,
#' derived from the m3d dataset to infer S. cerevisiae gene regulatory network.
#' The dataset was shorten to 3600 genes in order to limit the size of the
#' object
#' @format a matrix of 3600 genes by 247 samples
"SC_EXP_DATA"

#' SC_Test_data
#' A matrix of S. cerevisiae gene expression during diauxic shift (Brauer and
#' al.)
#' @format a matrix of 6028 genes by 13 samples during diauxic shift
"SC_Test_data"

#' SC_experiment_influence
#' A vector of influence computed from the first sample of SC_Test_data
#' @format a named numerical vector
"SC_experiment_influence"

#' iMM904
#' A S. cerevisiae genome-scale metabolic model as a modelOrg object
#' @format a modelOrg object as required by sybil. See sybilSBML for more
#' information on how to load other model.
"iMM904"

#' PredictedGeneState
#' Predicted gene state as obtained by the function
#' predict_linear_model_influence
#' @format a named vector
"PredictedGeneState"

#' ODcurveToMetCurve
#' List as obtained by the function ODcurveToMetCurve
#' @format List as obtained by the function ODcurveToMetCurve
"ODcurveToMetCurve"

#' ODtoflux
#' List as obtained by the function ODCurveToFluxCurves
#' @format List as obtained by the function ODcurveToMetCurve
"ODtoflux"
