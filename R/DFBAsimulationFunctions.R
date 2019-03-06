#' Euler step biomass
#'
#' This function updates the biomass given the biomass at the
#' previous time step and the size of the time step.
#'
#' @param biomass_t0  biomass in the current time step
#' @param rate  the growth rate for exponential growth
#' @param time_step time step given by t1-t0
#'
#' @return biomass after time.step=t1-t0: biomass.t1 = biomass.t0 *
#' exp(rate * time.step)
#'
#' @seealso euler_step_metabolites, update_system_state
#' @keywords internal

euler_step_biomass <- function(biomass_t0, rate, time_step){
    return(biomass_t0 * exp( rate * time_step ))
}


#' Compute the metabolite concentrations at the next time point
#'
#' @param met_concentrations_t0 metabolites concentration in the current time
#' step
#' @param fluxes fluxes states in the current time step for the exchange
#'  metabolites
#' @param biomass_t0  biomass in the current time step
#' @param rate the growth rate for exponential growth
#' @param time_step time step given by t1-t0
#'
#' @return metabolites concentrations at the next time point
#' @seealso euler_step_biomass, update_system_state
#' @keywords internal

euler_step_metabolites <- function(met_concentrations_t0,
                                   fluxes, biomass_t0, rate, time_step){
    metabolites_concentrations_t1 <- met_concentrations_t0 -
        fluxes / rate * biomass_t0 * ( 1 - exp ( rate * time_step ) )
    metabolites_concentrations_t1[metabolites_concentrations_t1 < 0] <- 0
    return(metabolites_concentrations_t1)
}

#' Flux balance analysis solution for the model
#'
#' Here we just perform a flux balance analysis solution and store the objective
#'  and the fluxes
#'
#' @param model An object of class modelOrg, the metabolic model.
#'
#' @return  list with the element
#' fluxes - fluxes of the flux balance analysis solution
#' obj- value of the objective function for the flux balance analyis solution.
#' @keywords internal
FBA_step <- function(model){

    Esol <- sybil::optimizeProb(model, algorithm = "fba")

    temp <- methods::slot(Esol, "fluxdist")
    list(fluxes = as.vector(methods::slot(temp, "fluxes")),
         obj = methods::slot(Esol, "lp_obj"))
}

#' Update the fluxes constraints given the metabolite concentrations
#'
#' @param model An object of class modelOrg, the metabolic model.
#' @param met_fluxes_indexes Indexes of the metabolites fluxes
#' @param met_concentrations_t0 Metabolites concentrations at t0
#' @param biomass_t0 Biomasss at t0
#' @param time_step time_step studied
#'
#' @return Return the updated model

update_uptake_fluxes_constraints_metabolites <- function(model,
                                                         met_fluxes_indexes,
                                                         met_concentrations_t0,
                                                         biomass_t0, time_step){
    #we check if the model is in irreversible form
    if ( !methods::.hasSlot(model, "irrev") ){
        #here the tricky part is dealing with negative bounds
        uptake_bounds <- met_concentrations_t0 / ( biomass_t0 * time_step )
        uptake_bounds <- ifelse( uptake_bounds < 1e-9, 0, uptake_bounds )
        sybil::lowbnd(model)[met_fluxes_indexes] <- ifelse(
          sybil::lowbnd(model)[met_fluxes_indexes] == 0,
          no = ifelse(
             abs(sybil::lowbnd(model)[met_fluxes_indexes]) < abs(uptake_bounds),
             yes = sybil::lowbnd(model)[met_fluxes_indexes],
             no = -abs(uptake_bounds)),
          yes = -uptake_bounds)
        return(model)
    }else{
        #we find the uptake reactions belonging to the extrogenous metabolites
        tmp <- sybil::findExchReact(model)
        uptake_reactions <- intersect(tmp@react_pos[tmp@uptake],
                                      met_fluxes_indexes)
        uptake_bounds <- met_concentrations_t0 / ( biomass_t0 * time_step )

        uptake_bounds <- ifelse( uptake_bounds < 1e-9, 0, uptake_bounds )
        sybil::uppbnd(model)[met_fluxes_indexes[uptake_reactions]] <- ifelse(
            sybil::uppbnd(model)[met_fluxes_indexes[uptake_reactions]] == 0,
            no = apply(cbind(
                abs(sybil::uppbnd(model)[met_fluxes_indexes[uptake_reactions]]),
                abs(uptake_bounds[uptake_reactions])), 1, min),
            yes = uptake_bounds[uptake_reactions])

    }
    return(model)
}

#' Update the fluxes constraints to simulate TF KO or overexpression
#'
#' Update the constraints according to the influence & regulatory network for a
#' single KO or over-expression
#'
#' @param model An object of class modelOrg, the metabolic model.
#' @param coregnet Object of class CoRegNet, containing the regulatory and
#' coregulatory interactions.
#' @param regulator_table A data.frame containing 3 columns: "regulator",
#' "influence","expression" containing respectively the name of a TF present in
#' the CoRegNet object as string, its influence in the condition of interest as
#' a numerical and an expression factor of 0 for a KO, or an integer >1 for an
#'  overexpression
#' @param aliases Optional, a data.frame containing the gene names used in the
#' metabolic model and the aliases to use to match the regulatory network
#' @return Return the model with updated bounds
#' @examples data("SC_GRN_1")
#' data("iMM904")
#' data("aliases_SC")
#' regulator_table <- data.frame("regulator" = "MET32",
#'                               "influence" = -1.20322 ,
#'                               "expression" = 3,
#'                               stringsAsFactors = FALSE)

#' model_TF_KO_OV_constraints <- update_fluxes_constraints_influence(
#' model= iMM904,coregnet = SC_GRN_1,regulator_table =  regulator_table,
#' aliases = aliases_SC)
#' @export



update_fluxes_constraints_influence <- function(model, coregnet,
                                                regulator_table,
                                                aliases) {

    regulator_table<-PrepareRegulatorTable(regulator_table = regulator_table,
                                               coregnet= coregnet,
                                               aliases= aliases)
    for(i in seq_len(dim(regulator_table)[1])){

    TF <- regulator_table$regulator[i]
    influence <- regulator_table$influence[i]
    expression_factor <- regulator_table$expression[i]

    target_genes_activating <- CoRegNet::targets(coregnet,
                                                 regulator = TF,
                                                 type = ("activating"))
    target_genes_repressing <- CoRegNet::targets(coregnet,
                                                 regulator = TF,
                                                 type = ("repressing"))

    if (!is.null(aliases)){
        target_genes_activatingDF <- merge(aliases,
                                         data.frame(target_genes_activating),
                                         by.x = "geneName_GRN",
                                         by.y = "target_genes_activating")
        target_genes_activating <- target_genes_activatingDF$geneName_model
        target_genes_repressingDF <- merge(aliases,
                                         data.frame(target_genes_repressing),
                                         by.x = "geneName_GRN",
                                         by.y = "target_genes_repressing")
        target_genes_repressingDF <- target_genes_repressingDF$geneName_model
        }

    if (is.null(expression_factor)|expression_factor==1) {
        # the overexpression factor N will take its default value of 1 and KO is
        # carried out
            model <- update_fluxes_constraints_GRegulation(model,
                       as.character(target_genes_activating),
                       function(x) perturbation_function(old_bound = x,
                                                     tf_influence = influence))
            # only if its two sided
            model <- update_fluxes_constraints_GRegulation(model,
                       as.character(target_genes_repressing),
                       function(x) perturbation_function_twosided(old_bound = x,
                                                     tf_influence = influence))

    } else {
            # The activated target will have a new value of
            # old_bound*(logistic(tf_influence, expression_factor))
            model <- update_fluxes_constraints_GRegulation(model,
                       as.character(target_genes_activating),
                       function(x) perturbation_function_twosided(old_bound = x,
                                                     tf_influence = influence,
                                                     expression_factor))
            # only if its two sided.
            model <- update_fluxes_constraints_GRegulation(model,
                       as.character(target_genes_repressing),
                       function(x) perturbation_function(old_bound = x,
                                                     tf_influence = influence,
                                                     expression_factor))
        }
    rm(TF,influence,expression_factor,target_genes_activating,
       target_genes_activatingDF,target_genes_repressingDF,
       target_genes_repressing)}
    return(model)
}

#' Update the fluxes constraints to simulate gene KO or overexpression
#'
#' Update the constraints of the reactions associated with the knock-out or
#' overexpressed gene
#'
#' @param model An object of class modelOrg, the metabolic model.
#' @param gene_table A data.frame containing 2 columns: "gene",
#' "expression" containing respectively the name of a gene present in
#' the CoRegNet object as string and an expression factor of 0 for a KO,
#' or an integer >1 for an overexpression
#' @param aliases Optional. A data.frame containing the gene names used in the
#' metabolic model and the aliases to use to match the regulatory network
#' @return Return the model with updated bounds
#' @examples
#' data("iMM904")
#' data("aliases_SC")
#' gene_table <- data.frame("gene" = c("YGL202W","YIL162W"),
#' "expression" =c(2,0), stringsAsFactors = FALSE)
#'
#' model_gene_KO_OV_constraints <- update_fluxes_constraints_geneKOOV(
#' model= iMM904,
#' gene_table =  gene_table,
#' aliases = aliases_SC)
#' @export

update_fluxes_constraints_geneKOOV <- function(model,
                                               gene_table,
                                               aliases=NULL) {
    gene_table<-PrepareGeneTable(gene_table = gene_table,
                                 model = model,
                                 aliases= aliases)
    for(i in seq_len(dim(gene_table)[1])){
        gene <- gene_table$gene[i]
        expression_factor <- gene_table$expression[i]
        model <- update_fluxes_constraints_GRegulation(model,as.character(gene),
                                        function(x) return(x*expression_factor))
    }
    return(model)

}
#' Update fluxes constraints
#'
#' Update the fluxes constraints to simulate a knockout/over-expression using a
#' penalty function (logistic) over the targets of a given TF
#' @param model An object of class modelOrg, the metabolic model.
#' @param target_genes_names Targets which will be affected by the knock-out
#' @param p_function Perturbation function
#'
#' @return Return the model with updated bounds
#' @keywords internal

update_fluxes_constraints_GRegulation <- function(model,
                                                  target_genes_names,
                                                  p_function){

    if (!any(sybil::allGenes(model) %in% target_genes_names)){
        return(model)
        }
    target_indexes <- which( sybil::allGenes(model) %in% target_genes_names)
    indexes <- sybil::geneDel(model, target_indexes)
    sybil::lowbnd(model)[indexes] <- p_function(sybil::lowbnd(model)[indexes])
    sybil::uppbnd(model)[indexes] <- p_function(sybil::uppbnd(model)[indexes])
    return(model)
}

#' Perturbation function
#'
#' @param old_bound Current value of the flux bound
#' @param tf_influence Influence of the chosen TF in the studied condition
#' @param expression_factor Numerical value by which the current value of
#' the bounds can be multiply with
#' @return New bound value
#' @keywords internal


perturbation_function <- function(old_bound,
                                  tf_influence,
                                  expression_factor=1){
    old_bound * (expression_factor - logistic(tf_influence,
                                                  expression_factor))
}


perturbation_function_twosided <- function(old_bound,
                                           tf_influence,
                                           expression_factor=1){
    old_bound * (logistic(tf_influence, expression_factor))
}


logistic <- function(x, expression_factor){
    expression_factor / ( 1 + exp(-x))
}



#' Update the fluxes bounds given the metabolites concentrations
#'
#' Updates the fluxes values by performing flux balance analysis, it is
#' necessary to indicate the indexes of the extraneous metabolites and the
#' biomass reaction. It performs an fba step and selects the given fluxes
#'
#' @param model An object of class modelOrg, the metabolic model.
#' @param met_fluxes_indexes index of the fluxes corresponding to uptake of
#' extrogenous metabolites
#' @param biomass_flux_index index of the flux corresponding to the biomass
#' reaction.
#'
#' @return a list of elements:
#' fluxes : array containing the fluxes values for the fba solution.
#' biomass_yield:  flux value for the corresponding biomass reaction
#' met_fluxes:  fluxes corresponding to the uptake reactions of the extraneous
#' metabolites
#' obj: value of the objective function obtained by the flux balance analysis
#' model.
#' @keywords internal
update_fluxes_state <- function(model, met_fluxes_indexes,
                                biomass_flux_index){
    FBA_solution <- FBA_step(model)
    if (methods::.hasSlot(model, "irrev")){
        tmp <- sybil::findExchReact(model)
        FBA_solution$fluxes[intersect(met_fluxes_indexes,
                                      tmp@react_pos[tmp@uptake])] <-
            -1 * FBA_solution$fluxes[intersect(met_fluxes_indexes,
                                               tmp@react_pos[tmp@uptake])]
    }
    met_fluxes <- FBA_solution$fluxes[met_fluxes_indexes]
    fluxes <- FBA_solution$fluxes
    biomass_yield <- FBA_solution$fluxes[biomass_flux_index]

    return(list(fluxes = fluxes,
                biomass_yield = biomass_yield,
                met_fluxes = met_fluxes,
                objective = FBA_solution$obj))
}


#' Update the fluxes, metabolites concentration and biomass in the system
#'
#' @param flux_state Current fluxes states
#' @param biomass_t0 Initial biomass
#' @param met_concentrations_t0 Metabolites concentration at t0
#' @param time_step Time step
#' @param Biomass_Yield_To_Rate Function to convert biomass to rate.
#' Default is f(x)=1*x
#' @return Return updated biomass and metabolites concentrations
#' @keywords internal

update_system_state <- function(flux_state,
                                biomass_t0,
                                met_concentrations_t0,
                                time_step,
                                Biomass_Yield_To_Rate){
    #in case we have another rate to biomass function
    rate <- ifelse(missing(Biomass_Yield_To_Rate),
                   flux_state$biomass_yield,
                   Biomass_Yield_To_Rate(flux_state$biomass_yield))

    #in case of rate equals zero, we may be able to remove this later
    rate <- ifelse(rate == 0, 1e-12, rate)

    biomass_t1 <- euler_step_biomass(biomass_t0, rate, time_step)

    met_concentrations_t1 <- euler_step_metabolites(met_concentrations_t0,
                                                    flux_state$met_fluxes,
                                                    biomass_t0,
                                                    rate,
                                                    time_step)
    list(rate = rate,
         biomass_t1 = biomass_t1,
         met_concentrations_t1 = met_concentrations_t1)
}


PrepareRegulatorTable<-function(regulator_table,coregnet,aliases){
    #In case of expression_factor = 1 remove the corresponding line
    regulator_table$expression<-round(regulator_table$expression)
    if(any(grepl('^1$',regulator_table$expression))){
       regulator_table <- regulator_table[!grepl(1,regulator_table$expression),]
    }
    # Replace 0 by 1 to carry out a KO in update_fluxes_constraints_influence
    regulator_table$expression<-as.numeric(gsub(0,1,regulator_table$expression))
    res <- regulator_table
    if ( dim(regulator_table)[1]>0 ){
        if ( !is.null(coregnet) ){
            if (any(regulator_table$regulator %in% names(CoRegNet::
                                                        regulators(coregnet)))){
                InGRN<- regulator_table$regulator  %in% names(
                    CoRegNet::regulators(coregnet))
                myTF <- as.character(regulator_table$regulator[InGRN])
                myInfluence <- regulator_table$influence[InGRN]
                myExpression <- regulator_table$expression[InGRN]
                rm(InGRN)
                res<-data.frame(regulator = myTF,
                                influence = myInfluence,
                                expression = myExpression,
                                stringsAsFactors = FALSE)
                return(res)
            }
        }
    }
    else{
        return(NULL)
    }
    return(res)
}


PrepareGeneTable<-function(gene_table,model,aliases=NULL){
    gene_table$expression<-round(gene_table$expression)
    if(!is.null(aliases)){
        myGenes<-gene_table$gene
        if(all(myGenes %in% aliases$geneName_model)){
            gene_table<-gene_table[gene_table$gene %in% allGenes(model),]
            return(gene_table)
        }
        if(all(myGenes %in% aliases$geneName_GRN)){
            myGenes_DF<-merge(aliases,
                              data.frame(myGenes),
                              by.x = "geneName_GRN",
                              by.y = "myGenes")
            myGenes<-myGenes_DF$geneName_model
            gene_table$gene<-myGenes
            gene_table<-gene_table[gene_table$gene %in% allGenes(model),]
            return(gene_table)
        }
        else{
            aliases_vector<-as.vector(aliases$geneName_model)
            names(aliases_vector)<-as.vector(aliases$geneName_GRN)

            Rename<-function(ToRenamed,AssociationVector){
                if (is.na(AssociationVector[ToRenamed])== TRUE){
                    return(ToRenamed)}
                else{return(AssociationVector[ToRenamed])}
            }
            gene_table$gene<-unname(unlist(lapply(gene_table$gene,Rename,
                                                  aliases_vector)))
            gene_table<-gene_table[gene_table$gene %in% allGenes(model),]
            return(gene_table)}
    }
    gene_table<-gene_table[gene_table$gene %in% allGenes(model),]
    return(gene_table)
}

#' Single simulation step
#'
#' Single simulation step in which fluxes are reconstrained according to
#' metabolite concentrations, then given the continuous evaluation of the
#'  gpr rules and the softplus function of the gene regulatory state.
#' @param model An object of class modelOrg, the metabolic model.
#' @param coregnet  Optional, object of class CoRegNet object containing
#' the information about regulatory and coregulatory relationships
#' @param metabolites data frame of metabolites names
#' @param met_concentrations_t0 data frame of metabolites concentrations at t0,
#' before performing a time step
#' @param biomass_t0 biomass at t0 before performing a step
#' @param regulator_table A data.frame containing 3 columns: "regulator",
#' "influence","expression" containing respectively the name of a TF present in
#' the CoRegNet object as string, its influence in the condition of interest as
#' a numerical and an expression factor of 0 for a KO, or an integer >1 for an
#'  overexpression
#' @param gene_table A data.frame containing 2 columns: "gene",
#' "expression" containing respectively the name of a gene present in
#' the CoRegNet object as string and an expression factor of 0 for a KO,
#' or an integer >1 for an overexpression
#' @param time_step size of the time step to perform; that is t1-t0
#' @param gene_state  data frame with rows being gene names and columns being
#' the gene expression or any other continuous value representing metabolites
#'  activity to be evaluated using the gpr rules
#' @param softplus_parameter Softplus parameter identify through calibration.
#' Default to 0.
#' @param aliases Optional. A data.frame containing the gene names currently
#' used in the network under the colname "geneName" and the alias under the
#' colnames "alias"
#' @param biomass_flux_index index of the flux corresponding to the biomass
#' reaction.
#' @seealso Simulation
#' @return list of:
#' fluxes: fluxes of the resulting fba soulution to the metabolic and genetic
#' constraints
#'  biomass_yield: biomass yield that is used as proxy for the growth rate in
#'  the dynamic flux balance analysis solution. Corresponds to the flux of the
#'  biomass reaction of the model.


Simulation_Step <- function(model, coregnet, metabolites,
                            met_concentrations_t0,
                            biomass_t0,
                            regulator_table,
                            gene_table,
                            time_step,
                            gene_state,
                            softplus_parameter,
                            aliases,
                            biomass_flux_index){
    message("simulation step")
    #we update the fluxes constraints according to the amount of metabolites
    model <- update_uptake_fluxes_constraints_metabolites(model = model,
            met_fluxes_indexes = as.integer(metabolites$flux_position),
            met_concentrations_t0 = met_concentrations_t0,
            biomass_t0 = biomass_t0,
            time_step = time_step)

    soft_plus_positive <- NULL
    soft_plus_negative <- NULL
    #we update the gene constraints according to the softplus function of the
    #continuous gpr rules.
    if (!missing(gene_state) & !is.null(gene_state) ){
        predicted <- as.matrix(gene_state$State)
        rownames(predicted) <- gene_state$Name
        model_coreg <- coregflux_static(model,
                                        predicted_gene_expression = predicted,
                                        gene_parameter = softplus_parameter,
                                        aliases = aliases)
        model <- model_coreg$model
        soft_plus_positive <- model_coreg$softplus_positive
        soft_plus_negative <- model_coreg$softplus_negative
    }
    # we apply constraint to simulate the knock_out of influent TF

    if ( !is.null(regulator_table) ){
             model <- update_fluxes_constraints_influence(model = model,
                                    coregnet = coregnet,
                                    regulator_table = regulator_table,
                                    aliases = aliases)

    }
    if(!is.null(gene_table)){
        model<- update_fluxes_constraints_geneKOOV(model=model,
                                                   gene_table = gene_table,
                                                   aliases = aliases)
    }
    # we update the state of the fluxes by solving the fba problem given by the
    # previous constraints
    flux_state <- update_fluxes_state(model = model,
                                      met_fluxes_indexes =
                                          metabolites$flux_position,
                                      biomass_flux_index)
    # we update the system state by solving an euler step for metabolites and
    # biomass given the previously computed fluxes
    system_state <- update_system_state(flux_state = flux_state,
                                        biomass_t0 = biomass_t0,
                                        met_concentrations_t0 =
                                          met_concentrations_t0,
                                        time_step = time_step,
                                     Biomass_Yield_To_Rate = function(x){1 * x})
    # if the model is irreversible, we have to make sure that the metabolite
    # update for the forward and backward reactions gets summed
    if (methods::.hasSlot(model, "irrev")){

        tmp <- data.frame(name = names(system_state$met_concentrations_t1),
                        t1 = system_state$met_concentrations_t1,
                        t0 = met_concentrations_t0)

        tmp <- plyr::ddply(tmp, "name", transform, t1 = tmp$t0 + sum(
            tmp$t1 - tmp$t0) )
        system_state$met_concentrations_t1 <- tmp$t1
    }

    list(fluxes = flux_state$fluxes,
         biomass_yield = flux_state$biomass_yield,
         met_fluxes = flux_state$met_fluxes,
         rate = system_state$rate,
         biomass_t1 = system_state$biomass_t1,
         met_concentrations_t1 = system_state$met_concentrations_t1,
         objective = flux_state$objective,
         metabolites = metabolites$name,
         soft_plus_positive = soft_plus_positive,
         soft_plus_negative = soft_plus_negative

    )
}

#' Simulation using Dynamic Flux balance analysis over time as in \cite{varma}
#' @param model An object of class modelOrg, the genome-scale metabolic model
#' (GEM).
#' @param time Timepoints at which the flux balance analysis
#' solution will be evaluated.
#' @param metabolites A data.frame containing the extraneous metabolites and the
#' initial concentrations
#' @param initial_biomass The value of the biomass at the beginning of the
#' simulation
#' @param coregnet Object of class CoRegNet, containing the regulatory and
#' coregulatory interactions.
#' @param regulator_table A data.frame containing 3 columns: "regulator",
#' "influence","expression" containing respectively the name of a TF present in
#' the CoRegNet object as a string, its influence in the condition of interest
#' as a numerical and an expression factor of 0 for a KO, or an integer >1 for
#' an overexpression
#' @param gene_table A data.frame containing 2 columns: "gene" and "expression"
#' containing respectively the name of a gene present in the modelOrg as a
#' string and an expression factor of 0 for a KO, or an integer >1 for
#' an overexpression
#' @param gene_state_function Function to obtain the gene state for a given
#' subset of gene
#' @param time_step_fba_bounds Bounds for the fba problem at each time point,
#' overrides any other form of constraining for a given flux.
#' @param softplus_parameter the softplus parameter identify through calibration
#' @param aliases Optional. A data.frame containing the gene names currently
#' used in the network under the colname "geneName" and the alias under the
#' colnames "alias"
#' @param biomass_flux_index index of the flux corresponding to the biomass
#' reaction.
#' @import CoRegNet
#' @import sybil
#' @return Return a list containing the simulation information such as the
#' objective_history, fluxes_history, met_concentration_history, biomass_history
#' @export
#' @details The simulation function allows the user to run several kind of
#' simulations based on the provided arguments. When providing only the GEM,
#' time, initial biomass and the metabolites, a classical dFBA is carried out.
#' To integrate the gene expression to the GEM, the gene_state_function must
#' be provided while if the user wants to simulate a TF knock-out or
#' overexpression, then a coregnet object and the regulator table should also be
#'  provided. See the vignette and quick-user guide for more examples.
#' @examples
#' data("SC_GRN_1")
#' data("SC_EXP_DATA")
#' data("SC_experiment_influence")
#' data("iMM904")
#' data("aliases_SC")
#' data("PredictedGeneState")
#'
#' metabolites<-data.frame("name"=c("D-Glucose","Glycerol"),
#'                         "concentrations"=c(16,0))
#'
#' result_without_any_constraint<-Simulation(iMM904,time=seq(1,10,by=1),
#'                    metabolites,
#'                    initial_biomass=0.45,
#'                    aliases = aliases_SC)
#'
#' GeneState<-data.frame("Name"=names(PredictedGeneState),
#'                     "State"=unname(PredictedGeneState))
#'
#' result<-Simulation(iMM904,time=seq(1,10,by=1),
#'                    metabolites,
#'                    initial_biomass=0.45,
#'                    gene_state_function=function(a,b){GeneState},
#'                    aliases = aliases_SC)
#'
#' result$biomass_history
Simulation <- function(model, time = c(0,1), metabolites, initial_biomass,
                       biomass_flux_index =
                           CoRegFlux::get_biomass_flux_position(model),
                       coregnet = NULL,
                       regulator_table = NULL,
                       gene_table = NULL,
                       gene_state_function = NULL,
                       time_step_fba_bounds = NULL,
                       softplus_parameter = 0,
                       aliases = NULL){

    if ( missing(biomass_flux_index) ){
        message(paste("Default biomass flux index use is",biomass_flux_index,
                      "corresponding to ",
                      model@react_name[biomass_flux_index]))
    }
    if(!is.null(metabolites$rates)){
        model <- adjust_constraints_to_observed_rates(model,
                                          metabolites_with_rates = metabolites)
    }
    #we compute the time step sizes according to the vector of time points
    time_step <- diff(time)

    exchange_met = CoRegFlux::build_exchange_met(model)
    # We find the exchange fluxes for the extraneous metabolites, gets tricky
    #for irreversible models
    metabolites <- get_metabolites_exchange_fluxes(metabolites = metabolites,
                                                   model = model,
                                                   exchange_met)

    #we initialize the arrays containing the history of the values
    objective_history <- vector(length = length(time_step) )

    fluxes_history <- matrix( ncol = length(time_step),
                              nrow = sybil::react_num(model) )
    rownames(fluxes_history) <- sybil::react_id(model)

    met_concentration_history <- matrix(ncol = length(time_step) + 1,
                                      nrow = length(metabolites$concentrations))
    rownames(met_concentration_history) <- metabolites$name

    met_fluxes_history <- matrix(ncol = length(time_step),
                                      nrow = length(metabolites$concentrations))

    rate_history <- vector(length = length(time_step))

    biomass_history <- vector(length = length(time_step) + 1)


    met_concentration_history[, 1] <- metabolites$concentrations
    biomass_history[1] <- initial_biomass
    gene_state_history <- vector(length = length(time_step) + 1)

    #we initialize the gene state vector in case of gene_state function
    if ( !is.null(gene_state_function) ){
        gene_state <- gene_state_function()

        gene_state_history[1] <- list(gene_state)
    }else{
        gene_state <- NULL
    }
    soft_plus_positive <- list()
    soft_plus_negative <- list()
    #we perform a simulation step for each time point
    for (i in seq_len(length (time) - 1)){
        if ( is.matrix(time_step_fba_bounds) ){
            if ( dim(time_step_fba_bounds)[2] == (length(time) - 1)){
                sybil::lowbnd(model) <- ifelse(time_step_fba_bounds[, i] < 0,
                                              time_step_fba_bounds[, i],
                                              lowbnd(model))
                sybil::uppbnd(model) <- ifelse(time_step_fba_bounds[, i] > 0,
                                              time_step_fba_bounds[, i],
                                              uppbnd(model))
            }
        }
        step_results <- Simulation_Step(model = model,
                                        metabolites = metabolites,
                                        coregnet = coregnet,
                                        met_concentrations_t0 =
                                            met_concentration_history[, i],
                                        biomass_t0 = biomass_history[i],
                                        regulator_table =
                                            regulator_table,
                                        gene_table = gene_table,
                                        time_step = time_step[i],
                                        gene_state = gene_state,
                                        softplus_parameter = softplus_parameter,
                                        aliases, biomass_flux_index)

        if ( !is.null(gene_state_function) ){
            gene_state_history[i] <- list(gene_state)
            gene_state <- gene_state_function(gene_state,
                                    step_results$met_concentrations_t1)
        }
        objective_history[i] <- step_results$objective
        fluxes_history[, i] <- step_results$fluxes
        met_concentration_history[, i + 1] <- step_results$met_concentrations_t1
        met_fluxes_history[, i] <- step_results$met_fluxes
        rate_history[i] <- step_results$rate
        biomass_history[i + 1] <- step_results$biomass_t1
        soft_plus_positive[[i]] <- step_results$soft_plus_positive
        soft_plus_negative[[i]] <- step_results$soft_plus_negative
    }
    result <- list(objective_history = objective_history,
                 metabolites = step_results$metabolites,
                 fluxes_history = fluxes_history,
                 met_concentration_history = met_concentration_history,
                 met_fluxes_history = met_fluxes_history,
                 rate_history = rate_history, biomass_history = biomass_history,
                 time = time,
                 gene_state_history = gene_state_history,
                 soft_plus_positive = soft_plus_positive,
                 soft_plus_negative = soft_plus_negative
    )
    return(result)

}
