#' Get intervals of flux variability (FVA) from an observed growth rate
#' @param model a genome-scale metabolic model of class modelorg
#' @param observed_growth_rate a numerical value for the observed growth rate
#' @param metabolites_rates Optional. a data.frame consisting of the name of the
#' metabolites, their concentrations and rates in mmol/gDW/h to adjust the model
#' uptake rates. The column name must be "name", "concentrations","rates"
#' @param biomass_flux_index Optional. Index of the biomass flux as returned by
#' \code{get_biomass_flux_position()}
#' @return Return the interval of fluxes values compatible with the observed
#' growth rate through flux variability analysis
#'
#' @examples
#' data("iMM904")
#' metabolites_rates<-data.frame(name=c("D-Glucose","Glycerol"),
#'                          concentrations=c(16,0),rates=c(-2.81,-8.01))
#'
#' FluxesVarFromObs<-get_fva_intervals_from_observations(iMM904,0.205,
#' metabolites_rates=metabolites_rates)
#' @export

get_fva_intervals_from_observations <- function(model,
                                                observed_growth_rate,
                                                metabolites_rates = NULL,
                                                biomass_flux_index =
                                            get_biomass_flux_position(model)) {

    sybil::uppbnd(model)[biomass_flux_index] <- observed_growth_rate

    if (!is.null(metabolites_rates)) {
        model <- adjust_constraints_to_observed_rates(model, metabolites_rates,
                                        exchange_met=build_exchange_met(model))
    }

    ts <- sybil::fluxVar(model)

    matrix(data = cbind(ts@lp_obj[seq_len(model@react_num)],
                        ts@lp_obj[(1 + model@react_num):(model@react_num * 2)]),
           ncol = 2,
           dimnames = list(model@react_id, c("min", "max")))
}

#' Get fluxes balance from an observed growth rate
#'
#' @param model a genome-scale metabolic model of class modelorg
#' @param observed_growth_rate a numerical value for the observed growth rate
#' @param metabolites_rates Optional, a data.frame consisting of the name of the
#' metabolites, their concentrations and rates in mmol/gDW/h to adjust the model
#' uptake rates. The column name must be "name", "concentrations","rates"
#' @param biomass_flux_index Optional. Index of the biomass flux as returned by
#' \code{get_biomass_flux_position()}
#' @param backward_fluxes Optional, only relevant for irreversible model
#' @param forward_fluxes Optional, only relevant for irreversible model
#' @return Return fluxes values compatible with the observed
#' growth rate through flux balance analysis
#' @export
#'
#' @examples
#' data("iMM904")
#' metabolites_rates<-data.frame(name=c("D-Glucose","Glycerol"),
#'                          concentrations=c(16,0),rates=c(-2.81,-8.01))
#'
#' FluxesFromObs<-get_fba_fluxes_from_observations(iMM904,0.205,
#' metabolites_rates = metabolites_rates)
get_fba_fluxes_from_observations <- function(model,
                                             observed_growth_rate,
                                             metabolites_rates = NULL,
                                             biomass_flux_index =
                                             get_biomass_flux_position(model),
                                             backward_fluxes="_b",
                                             forward_fluxes="_f"){

    sybil::uppbnd(model)[biomass_flux_index] <- observed_growth_rate
    if (!is.null(metabolites_rates)) {
        model <- adjust_constraints_to_observed_rates(model, metabolites_rates,
                                       exchange_met = build_exchange_met(model),
                                       backward_fluxes,forward_fluxes)
    }

    ts <- sybil::optimizeProb(model)
    matrix(data = ts@fluxdist@fluxes, dimnames = list(sybil::react_id(model)))
}

#' Adjust the constraint of the model to observed rates
#'
#' @param model a genome-scale metabolic model of class modelorg
#' @param metabolites_with_rates is a data.frame consisting of the name of the
#' metabolites, their concentrations and rates in mmol/gDW/h.
#' The column name must be "name", "concentrations","rates"
#' @param exchange_met Optional. a data.frame as given by build_exchange_met
#' @param backward_fluxes Optional. Useful for irreversible model
#' @param forward_fluxes Optional. Useful for irreversible model
#' @return Return the model with updated bounds corresponding to the observed
#' rates provided
#'
#' @examples
#' data("iMM904")
#' metabolites_rates<-data.frame(name=c("D-Glucose","Glycerol"),
#'                          concentrations=c(16,0),rates=c(-2.81,-8.01))
#' iMM904_adjusted<-adjust_constraints_to_observed_rates(iMM904,
#' metabolites_rates)
#' @export
adjust_constraints_to_observed_rates <- function(model,
                                                 metabolites_with_rates,
                                                 exchange_met=
                                                     build_exchange_met(model),
                                                 backward_fluxes="_b",
                                                 forward_fluxes="_f" ) {

    metabolites <- get_metabolites_exchange_fluxes(model = model,
                                                metabolites =
                                                    metabolites_with_rates,
                                            exchange_met = exchange_met,
                                            backward_fluxes = backward_fluxes,
                                            forward_fluxes = forward_fluxes)
    metabolites_fluxes_indexes <- metabolites$flux_position

    uptake_bounds <- metabolites$rates
    if (is.null(uptake_bounds)){
        stop("No rates were provided. Please, make sure your metabolite
                data.frame fits the requirement")
    }
    if (!methods::.hasSlot(model, "irrev")) {
        for (i in seq_along(uptake_bounds)) {
                if (uptake_bounds[i] <= 0 ){
                    sybil::lowbnd(model)[metabolites_fluxes_indexes[i]] <-
                    uptake_bounds[i]
                }else{
                    sybil::uppbnd(model)[metabolites_fluxes_indexes[i]] <-
                        uptake_bounds[i]
                }
        }
    } else {
        for (i in seq_along(uptake_bounds)) {
            sybil::uppbnd(model)[metabolites_fluxes_indexes[i]] <-
                abs(uptake_bounds[i])
        }
    }
    return(model)
}
