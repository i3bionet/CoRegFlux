#' Transform a logic gpr rule string to  a list and then evaluate the
#' continuous version of the rules
#'
#' @param gpr A gene-protein rules from the genome-scale metabolic model
#' @param expression A numerical matrix corresponding to the gene expression.
#' Rownames should contain gene names/ids while samples should be in columns.
#'
#' @return A numerical value corresponding to the evaluation of the flux based
#' on the expression of the metabolic genes present in the gpr
#' @keywords internal

gpr_expression <- function(gpr, expression) {
    gpr <- gsub("[()]", "", gpr)
    gpr <- gsub("[[:space:]]", "", gpr)
    complex <- lapply(gpr, function(gpr) {
        unlist(strsplit(gpr, "or"))
    })
    genes <- lapply(complex, function(complex) {
        strsplit(complex, "and")
    })
    genes[lengths(genes) == 0] <- NA
    min_complex <- lapply(genes, function(gene) {
        lapply(gene, function(gene) {
            gene <- unlist(gene)
            gene <- gene[gene %in% rownames(expression)]

            if (length(gene) == 0) {
                min_complex <- NA
            } else {
                # AND
                min_complex <- min(rowMeans(expression, na.rm = TRUE)[gene],
                                   na.rm = TRUE)
                }
            return(min_complex)
        })
    })
    exp <- unlist(lapply(min_complex, function(min_complex) {
        # OR
        if (all(is.na(unlist(min_complex)))) {
            0
        } else {
            max(unlist(min_complex), na.rm = TRUE)
        }
    }))
    return(exp)
}

#' Continuous evaluation of the gpr rules
#'
#' This function was adapted from code present in \cite{exp2flux}.
#'  Unlike exp2flux, we consider NA gene expression as uninformative, thus terms
#'  involving these items will not be evaluated. Also, unlike exp2flux, our
#'  equivalence is OR(A,B) <- sum(A,B); AND(A,B) <- min(A,B),
#'  more in accordance with \cite{tiger}
#'
#' @param model A genome-scale metabolic model of class modelorg
#' @param expression A numerical matrix corresponding to the gene expression.
#' Rownames should contain gene names/ids while samples should be in columns.
#' @param scale wether to scale the gene expression to unit variance
#'
#' @return A metabolic model with lower and upper bound corresponding to the
#' continuous version of the rule evaluation and zero on unaffected fluxes.
#' @keywords internal


continuous_gpr <- function(model, expression, scale = FALSE) {
    model@lowbnd <- 0
    model@uppbnd <- 0
    exp <- gpr_expression(gpr = model@gpr, expression = expression)
    if (scale == TRUE) {
        exp <- round( ( exp / max(exp, na.rm = TRUE)), 6) * 1000
    }
    lb <- model@lowbnd
    ub <- model@uppbnd
    model@lowbnd <- -1 * exp
    model@lowbnd[!model@react_rev] <- 0
    model@uppbnd <- exp
    model@lowbnd[model@react_id %in% sybil::findExchReact(model)@react_id] <-
        lb[model@react_id %in% sybil::findExchReact(model)@react_id]
    model@uppbnd[model@react_id %in% sybil::findExchReact(model)@react_id] <-
        ub[model@react_id %in% sybil::findExchReact(model)@react_id]
    return(model)
}

#' Update the model using the provided gene regulatory network and expression
#'
#' \code{coregflux_static()} uses the gene states to update the fluxes bounds
#' from the metabolic model.
#'
#' @param model A genome-scale metabolic model of class modelorg
#' @param predicted_gene_expression The vector of predicted gene expression for
#' the genes present in the metabolic model as given by
#' \code{predict_linear_model_influence()}
#' @param gene_parameter Parameter of the softplus function
#' @param tol  Fluxes values below this threshold will be ignored.
#' @param aliases a data.frame containing the gene names currently used in the
#' network under the colname "geneName" and the alias under the colnames
#' "alias"
#' @return list containing: \item{model}{the metabolic model with the coregflux
#' constraints added}
#' \item{softplus_positive}{the results of evaluating ln(1+exp(gpr(x+theta)))
#'  where gpr() are the continuous version of the gpr rules applied to a set of
#'  gene expression x}
#' \item{softplus_negative}{the results of evaluating ln(1+exp(gpr(x+theta)))
#'  where gpr() are the continuous version of the gpr rules applied to a set of
#'  gene expression x}
#' @export
#' @examples
#' data("SC_GRN_1")
#' data("SC_experiment_influence")
#' data("SC_EXP_DATA")
#' data("aliases_SC")
#' data(iMM904)
#' data(PredictedGeneState)
#' static_list<-coregflux_static(iMM904,PredictedGeneState)

coregflux_static <- function(model,
                            predicted_gene_expression,
                            gene_parameter = 0,
                            tol = 1e-10,
                            aliases=NULL) {
    # this is the model that will be readjusted by coregflux
    coregflux_model <- model
    # here we adjust the bounds based on the expression results
    message("Adjusting bounds of metabolic model")

    tmp <- as.matrix(predicted_gene_expression)


    exp <- data.frame(expression = tmp, geneName_GRN = rownames(tmp))
    if ( !is.null(aliases) ) {
        exp <- suppressMessages(plyr::join_all(list(exp, aliases)))
    } else {
        exp <- cbind(exp, data.frame(rownames(tmp)))
        colnames(exp) <- c("expression", "geneName_GRN", "geneName_model")
    }

    tmp <- as.matrix(exp$expression)
    rownames(tmp) <- exp$geneName_model

    expression_data <- tmp
    # we evaluate the continuous version of the gpr rules
    expression_coregflux_bounds <- continuous_gpr(model, expression_data,
                                                  scale = FALSE)
    # we find the exchange reactions
    ex_idx <- sybil::react_pos(sybil::findExchReact(model))
    not_ex <- setdiff(seq_len(length(sybil::react_id(model))), ex_idx)
    # we find the reactions which lower bound was changed by the rules
    changed_idx_lower <- which(abs(sybil::lowbnd(expression_coregflux_bounds)) >
                                   tol)
    # we find the reactions which upper bound was changed by the rules
    changed_idx_upper <- which(abs(sybil::uppbnd(expression_coregflux_bounds)) >
                                   tol)
    # we change the reactions whose lowerbound was changed and are not exchange
    # reactions
    to_change_lower <- intersect(not_ex, changed_idx_lower)
    # we change the reactions whose uppbnd changed & are not exchange reactions
    to_change_upper <- intersect(not_ex, changed_idx_upper)
    # we select the valid fluxes to change for the upper bounds, if the model is
    # irreversible, this will be the only change made
    valid_fluxes <- intersect(to_change_upper, which(sybil::uppbnd(model) >= 0))
    # we evaluate the softplus over the bounds computed using the continuous
    # version of the gpr rules
    softplus_eval <- evaluate_softplus(gene_parameter,
                     sybil::uppbnd(expression_coregflux_bounds)[valid_fluxes])
    # we adjust the upper bounds to the new values
    sybil::uppbnd(coregflux_model)[valid_fluxes] <- softplus_eval
    # we store the results of the softplus evaluation
    softplus_result_positive <- data.frame(fluxes = valid_fluxes,
                                           softplus = softplus_eval)
    softplus_result_negative <- NULL
    if (!methods::.hasSlot(model, "irrev")) {
        # the same but for negative fluxes
        valid_fluxes <- intersect(to_change_lower, which(abs(
            sybil::lowbnd(model)) >= 0))

        softplus_eval <- -1 * evaluate_softplus(gene_parameter, abs(
            sybil::lowbnd(expression_coregflux_bounds)[valid_fluxes]))

        sybil::lowbnd(coregflux_model)[valid_fluxes] <- softplus_eval

        softplus_result_negative <- data.frame(fluxes = valid_fluxes,
                                               softplus = softplus_eval)
    }

    return(list(model = coregflux_model,
                softplus_positive = softplus_result_positive,
                softplus_negative = softplus_result_negative))

}

evaluate_softplus <- function(gene.parameter, X) {
    return(log(1 + exp(gene.parameter + X)))
}


#' Train a linear model and predict the gene expression from an experiment
#' influence
#' @param train_expression Gene expression of the training data set, not
#'  necessary if train_influence is supplied. Should be numerical matrix
#'  corresponding to the gene expression. Rownames should contain gene names/ids
#'   while samples should be in columns.
#' @param train_influence Regulator influence scores of the train data set.
#' @param experiment_influence Regulator influence scores for the condition of
#' interest as a named vector with the TF as names.
#' @param minTarget The minimum number of targets for a regulator to be
#' considered for actvity prediction when computing the influence.
#' Default set to 10
#' @param network CoRegNet object to be interrogated for building the linear
#' models
#'
#' @return The predicted gene expression levels compute from the linear model
#' and the experiment influence
#' @keywords internal

train_continuous_model <- function(train_expression,
                                   train_influence,
                                   minTarget=10,
                                  experiment_influence,
                                  network) {
    number_of_regulators <- sum(CoRegNet::targets(network) %in%
                                    rownames(train_expression)) /
        length(CoRegNet::targets(network))
    if (is.null(rownames(train_expression)))
        stop("Training expression matrix must have rownames")
    if (ncol(train_expression) != ncol(train_influence))
        stop("Training set expression and influence do not have the same number
             of samples (columns)")
    if (number_of_regulators <= 0.05) {
        warning("Less than 5% of the network targets are present in the training
                expression data set, make sure the expression matrix rownames
                names match  the CoRegNet names")
    }

    # prediction results
    predict_continuous <- matrix(nrow = dim(train_expression)[1], ncol = 1)
    rownames(predict_continuous) <- rownames(train_expression)
    for (i in seq_len(dim(train_expression)[1])) {
        temp2 <- data.frame(expression = train_expression[i, ])
        # we find the set of regulators
        regulators_mgene <- CoRegNet::regulators(network,
                                                 rownames(train_expression)[i])

        if (!any(is.na(regulators_mgene))) {
            # we create an explanatory variables matrix from the m3d influence
            if (length(intersect(regulators_mgene, rownames(train_influence)))
                != 0) {
                if (length(intersect(regulators_mgene, rownames(train_influence)
                                     )) == 1) {
                    # Solve a bug in case of only one regulator by making sure
                    # it won't take the format of a vector
                    temporal_exp <- as.matrix(train_influence[
                        rownames(train_influence) %in% regulators_mgene, ])
                    colnames(temporal_exp) <- intersect(regulators_mgene,
                                                    rownames(train_influence))
                }
                if (length(intersect(regulators_mgene, rownames(train_influence)
                                     )) > 1) {

                    temporal_exp <- t(train_influence[rownames(train_influence)
                                                      %in% regulators_mgene, ])
                }

                # the train matrix will be of expression ~ Influences of
                # regulators

                temporal_set <- cbind(temp2, temporal_exp)
                # we perform linear regression
                temporal_regression <- stats::lm(expression ~ ., temporal_set)
                # we use the experimental influences as explanatory variables
                # for prediction
                temporal_prediction <- as.data.frame(cbind(t(
                    experiment_influence[names(experiment_influence) %in%
                                             regulators_mgene])))

                # until here there is no difference between training a discrete
                # and a continuous model we perform multinomial regression to
                # predict gene state in {-1, 0, 1}
                predict_continuous[i, ] <- as.numeric(
                    as.character(stats::predict(temporal_regression,
                                                temporal_prediction)))
                # we additionally recover the attached probabilities
                rm(temporal_exp, temporal_prediction, temporal_regression,
                   temporal_set)
            } else {
                predict_continuous[i, ] <- NA
            }
        } else {
            predict_continuous[i, ] <- NA
        }
    }
    return(list(predict_continuous = predict_continuous))
}


#' Train a linear model
#'
#' Here we train a linear regression model of the form x= alpha + beta*I
#' where x is the gene expression of the metabolic genes of the train data set
#' train_expression, alpha is an intercept, I is the influence of the regulators
#'  of the training data set and beta are the coefficients.
#'
#' train_expression Gene expression of the training data set, not
#'  necessary if train_influence is supplied. Should be numerical matrix
#'  corresponding to the gene expression. Rownames should contain gene names/ids
#'   while samples should be in columns.
#' @param train_influence Optional. Regulator influence scores computed using
#' the function CoRegNet::regulatorInfluence for the training data set,
#' default minTarg = 10
#' @param network CoRegNet object use to build the linear model and to compute
#' the influence.
#' @param train_expression Gene expression of the training data set, not
#'  necessary if train_influence is supplied. Should be numerical matrix
#'  corresponding to the gene expression. Rownames should contain gene names/ids
#'   while samples should be in columns.
#' @return A linear model
#' @seealso predict_linear_model_influence

get_linear_model <- function(train_expression,
                             train_influence=regulatorInfluence(network,
                                                train_expression, minTarg = 10),
                             network) {
    number_of_regulators <- sum(CoRegNet::targets(network) %in%
                rownames(train_expression)) / length(CoRegNet::targets(network))
    if (is.null(rownames(train_expression)))
        stop("Training expression matrix must have rownames")
    if (ncol(train_expression) != ncol(train_influence))
        stop("Training set expression and influence do not have the same number
             of samples (columns)")
    if (number_of_regulators <= 0.05) {
        warning("Less than 5% of the network targets are present in the training
                expression data set, make sure the expression matrix rownames
                names match  the CoRegNet names")
    }
    train <- as.matrix(train_expression)
    rownames(train) <- rownames(train_expression)
    # lets get rid of all zeros rows

    # prediction results
    linear_models <- vector(mode = "list", length = dim(train)[1])

    for (i in seq_len(dim(train)[1])) {
        # for each row of expression, that is for each metabolic gene
        temp2 <- data.frame(expression = train[i, ])
        # we find the set of regulators
        regulators_mgene <- regulators(network, rownames(train)[i])
        if (!any(is.na(regulators_mgene))) {
            # we create an explanatory variables matrix from the influence
            if (length(intersect(regulators_mgene, rownames(train_influence)))
                != 0) {
                if (length(intersect(regulators_mgene,
                                     rownames(train_influence))) == 1) {
                    # Solve a bug in case of only one regulator by making sure
                    # it won't take the format of a vector
                    temporal_exp <- as.matrix(train_influence[
                        rownames(train_influence) %in% regulators_mgene, ])
                    colnames(temporal_exp) <- intersect(regulators_mgene,
                                                    rownames(train_influence))
                }
                if (length(intersect(regulators_mgene,
                                     rownames(train_influence))) > 1) {

                    temporal_exp <- t(train_influence[rownames(train_influence)
                                                      %in% regulators_mgene, ])
                }

                # the train matrix will be expression ~ Influences of regulators
                temporal_set <- cbind(temp2, temporal_exp)
                linear_models[[i]] <- list(gene = rownames(train)[i],
                                           linear_model = stats::lm(expression ~
                                ., temporal_set), predictors = regulators_mgene)
            } else {
                linear_models[[i]] <- list(gene = rownames(train)[i],
                                           linear_model = NA, predictors = NA)
            }
        } else {
            linear_models[[i]] <- list(gene = rownames(train)[i],
                                       linear_model = NA, predictors = NA)
        }

    }
    return(linear_models)
}

predict_continuous_model <- function(linear_model, experiment_influence) {
    predict_continuous <- matrix(nrow = length(linear_model), ncol = 1)
    rownames(predict_continuous) <- unlist(lapply(linear_model,
                                                 function(x) x$gene))
    for (i in seq_along(linear_model)) {
        trained_model <- linear_model[[i]]
        model_formula <- trained_model$linear_model
        if (all(!is.na(trained_model))) {
            if (!any(is.na(model_formula$coefficients))) {

                temporal_prediction <- as.data.frame(cbind(t(
                    experiment_influence[names(experiment_influence) %in%
                                             trained_model$predictors])))
                predict_continuous[i, ] <- as.numeric(as.character(
                    stats::predict(model_formula, temporal_prediction)))

            } else predict_continuous[i, ] <- NA
        } else predict_continuous[i, ] <- NA
    }
    return(list(predict_continuous = predict_continuous))
}

PrepareTrainDataInf<-function(network,
                              train_expression,
                              experiment_influence,
                              minTarget,verbose=0){
    if(verbose==1){
        message("Computing influence for the training dataset")
    }
    train_influence <- CoRegNet::regulatorInfluence(object = network,
                                                    expData = train_expression,
                                                    minTarg = minTarget)
    TF_to_remove <- setdiff(rownames(train_influence),
                            names(experiment_influence))
    train_influence<-train_influence[-which(
        rownames(train_influence) %in% TF_to_remove),]
    return(train_influence)
}

#' Predict the gene expression level based on condition-specific influence
#'
#' Build a linear model and use it to predict the gene expression level from the
#' influence of an experiment
#' @param network a coregnet object
#' @param experiment_influence Regulator influence scores for the condition of
#' interest as a named vector with the TF as names.
#' @param train_expression Gene expression of the training data set, not
#'  necessary if train_influence is supplied. Should be numerical matrix
#'  corresponding to the gene expression. Rownames should contain gene names/ids
#'   while samples should be in columns.
#' @param tol Fluxes values below this threshold will be ignored. Default
#' @param train_influence Optional, if is train_expression is provided.
#' An influence matrix as computed by the function \code{regulatorInfluence()} from
#' CoRegNet
#' @param model  A genome-scale metabolic model from a class modelOrg.
#' @param min_Target Optional. Use in case train_influence is not provided.
#' Default value = 10. See regulatorInfluence for more information.
#' @param aliases Optional, A two columns data.frame containing the name used in
#' the gene regulatory network and their equivalent in the genome-scale
#' metabolic model to allow the mapping of the GRN onto the GEM.
#' The colnames should be geneName_model and geneName_GRN
#' @param verbose Default to 0. Give informations about the process status
#' @return The predicted genes expressions/states
#' @examples data("SC_GRN_1")
#' data("SC_experiment_influence")
#' data("SC_EXP_DATA")
#' data("iMM904")
#' data("aliases_SC")
#' PredictedGeneState <- predict_linear_model_influence(network = SC_GRN_1,
#'                     experiment_influence = SC_experiment_influence,
#'                     train_expression = SC_EXP_DATA,
#'                     min_Target = 4,
#'                     model = iMM904,
#'                     aliases= aliases_SC)
#'
#' GeneState<-data.frame("Name"=names(PredictedGeneState),
#'                     "State"=unname(PredictedGeneState))
#'
#' @export

predict_linear_model_influence <- function(network,
                                           model,
                                        train_influence = regulatorInfluence(
                                           network,train_expression,min_Target),

                                           experiment_influence,
                                           train_expression,
                                           min_Target = 10,
                                           tol = 1e-10,
                                           #linear_model = NULL,
                                           aliases = NULL,
                                           verbose = 0) {

        if (!is.null(aliases)) {
            if(verbose==1){
                       message("Processing aliases ...")}
            aliases_metabolic <- merge(aliases, data.frame(
                                                geneName_model =
                                                    sybil::allGenes(model)),
                                       by="geneName_model"
                                       )
            train_expression_metabolic <- train_expression[
            rownames(train_expression) %in% aliases_metabolic$geneName_GRN, ]
        } else {
            train_expression_metabolic <- train_expression[
                rownames(train_expression) %in% sybil::allGenes(model), ]
        }
        if(verbose==1){
            message("Training linear regression model ...")}

        TF_to_remove <- setdiff(rownames(train_influence),
                              names(experiment_influence))
        if(length(TF_to_remove) !=0){
            train_influence<-train_influence[-which(
                rownames(train_influence) %in% TF_to_remove),]
        }

        train_results <- train_continuous_model(train_expression =
                                                   train_expression_metabolic,
                                               #threshold = discrete_threshold,
                                              train_influence = train_influence,
                                              experiment_influence =
                                                   experiment_influence,
                                              network = network)

    return(train_results$predict_continuous[stats::complete.cases(
        train_results$predict_continuous), ])
}
