#' Get biomass flux position
#'
#' @param model An object of class modelOrg, the genome scale metabolic model
#' @param biomass_reaction_id Default value "biomass"
#' @param biomass_reaction_name Optional, the react_name in the modelOrg under
#' which the biomass function can be found, such as "growth"
#' @return the position of the biomass generating reaction according the the
#' objective
#'in our case we had the biomass reactions for models iMM904 and iTO977
#' @export
#'
#' @examples data("iMM904")
#' get_biomass_flux_position(iMM904)
get_biomass_flux_position <- function(model, biomass_reaction_id = "biomass",
                                      biomass_reaction_name = NULL){
    if(!is.null(biomass_reaction_name)){
        result <- grep(biomass_reaction_name, sybil::react_name(model),
                       ignore.case = TRUE)
        if (length(result) != 1 ){
            if (length(result) > 1){
                stop(paste('More than one reaction match the react_name "',
                           biomass_reaction_name,'"',sep=""))
            }else{
                stop("Biomass reaction not found in the model")
            }
            return(result)
            }
        }
    result <- grep(biomass_reaction_id, sybil::react_id(model),
                  ignore.case = TRUE)
    if (length(result) > 1){
        stop(paste('More than one reaction match the react_id "',
                   biomass_reaction_id,'"',sep=""))
    }
    if (length(result) == 0){
        if (missing(biomass_reaction_name)){
        stop("Biomass reaction not found in the model")
    } else{
        result <- grep(biomass_reaction_name, sybil::react_name(model),
                      ignore.case = TRUE)
        print(result)
        if (length(result) != 1 ){
            if (length(result) > 1){
                stop(paste('More than one reaction match the react_name "',
                           biomass_reaction_name,'"',sep=""))
            }
            else{
                stop("Biomass reaction not found in the model")
                }
        }else{
            return(result)
            }
    }
    }
    return(result)
}

#' Convert metabolites name to the model equivalent
#' @param model A genome-scale metabolic model as a modelorg object
#' @param metabolites A data.frame containing the metabolites and their
#' concentrations
#' @param exchange_met A data.frame as build by the function build_exchange_met
#' @return A data.frame containing the exchange metabolite model id and the
#' equivalent name
#' @keywords internal

convert_metabolites_to_model_names <- function(metabolites, model,
                                               exchange_met=
                                           build_exchange_met(model)) {
    result <- merge(x = metabolites, y = exchange_met,
                    by.x = colnames(metabolites)[1],
                    by.y = colnames(exchange_met)[2])
    result <- result[stats::complete.cases(result), ]
    return(result)
    }

#' Build the exchange metabolite data.frame
#'
#' @param model An object of class modelOrg, the genome scale metabolic model
#'
#' @return a data.frame containing the exchange metabolite model id and the
#' equivalent name
#' @examples
#' data("iMM904")
#' exchanged_met<-build_exchange_met(iMM904)
#' head(exchanged_met)
#' @export


build_exchange_met <- function(model){
        exchanged_met <- sybil::findExchReact(model)
        df <- data.frame("model_name"=exchanged_met@met_id,
                         "exchange_metabolite_name" =
                             sybil::met_name(model)[exchanged_met@met_pos],
                         stringsAsFactors = FALSE )
    return(df)
}

#' Get metabolites exchange fluxes
#' @param metabolites A data.frame containing the names and concentrations of
#' metabolites
#' @param model An object of class modelOrg, the genome scale metabolic model
#' @param exchange_met A data.frame as build by the function build_exchange_met
#' @param backward_fluxes Optional parameter for irreversible model to indicate
#' backward fluxes
#' @param forward_fluxes Optional parameter for irreversible model to indicate
#' forward fluxes
#' @return a data.frame containing the exchange metabolite model id and the
#' equivalent name
#' @export
#'
#' @examples
#' data("iMM904")
#' metabolites<-data.frame("name"=c("D-Glucose","Glycerol"),
#'                         "concentrations"=c(16,0))
#' get_metabolites_exchange_fluxes(iMM904,metabolites)
get_metabolites_exchange_fluxes <- function(model,metabolites,
                                            exchange_met =
                                              build_exchange_met(model),
                                            backward_fluxes = "_b",
                                            forward_fluxes = "_f") {
    metabolites_model_names <- convert_metabolites_to_model_names(metabolites,
                                                    model,
                                                    exchange_met = exchange_met)

    exch_reactions <- sybil::findExchReact(model)
    flux_position <- exch_reactions@react_pos
    reaction_id <- exch_reactions@react_id
    metabolites_id <- exch_reactions@met_id
    metabolites_model_names$model_name <- as.factor(
        metabolites_model_names$model_name)
    model_exchange_fluxes <- data.frame(flux_position, reaction_id,
                                        metabolites_id)
    metabolites_model_names$metabolites_id <- metabolites_model_names$model_name
    full_met <- plyr::join(metabolites_model_names, model_exchange_fluxes)
    if (methods::.hasSlot(model,"irrev")){
        full_met<- unique(full_met)
        if (!is.null(full_met$rates)){
            met_positive_rate <- subset(full_met, full_met$rates > 0)
            met_negative_rate <- subset(full_met, full_met$rates <= 0)

            met_pos_f<- met_positive_rate[-grep(backward_fluxes,
                                                met_positive_rate$reaction_id),]
            met_neg_b<- met_negative_rate[-grep(forward_fluxes,
                                                met_negative_rate$reaction_id),]
            full_met<-rbind(met_pos_f,met_neg_b)
        }
    }
    return(full_met)
}
