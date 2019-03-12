# Transforms gpr rules vector in sum(min(...)) form
# examples gpr_expression_tmp(gpr(iMM904)[[3]])

gpr_expression_tmp <- function(gpr) {
  gpr <- gsub("[()]", "", gpr)
  gpr <- gsub("[[:space:]]", "", gpr)
  complex <- lapply(gpr, function(gpr) {
    unlist(strsplit(gpr, "or"))
  })
  genes <- lapply(complex, function(complex) {
    strsplit(complex, "and")
  })
  genes[lengths(genes) == 0] <- ""
  min_complex <- lapply(genes, function(gene) {
    sub = paste(sort(unlist(lapply(gene, function(g){
      if(length(g) > 1)
        paste0("min(",paste(sort(g), collapse=","),")")
      else if(length(g) == 1)
        paste(sort(g), collapse=",")
    }))), collapse=",")
    if(length(gene) > 1)
      paste0("sum(", sort(sub),")")
    else if(length(gene) == 1)
      sort(sub)
  })
  unlist(min_complex)
}

#2 possible input shapes:
## - aggregator(sub1, sub2, sub3) each subs being a possible input for parseSubs
## - litteral without aggregation nor ","
# @example parseSubs("sum(YDR261C,YLR300W,YOR190W)")
parseSubs <- function(str){
  if(!grepl(",", str)) NULL
  else{
    chars = unlist(strsplit(str, ""))
    acc = ""
    level = 0
    blocks = list()
    for(i in seq_along(chars)){
      acc = paste0(acc, chars[i])
      if(chars[i] == ")"){
        level = level - 1
      }else if(chars[i] == "("){
        if(level == 0)
          acc = ""
        level = level + 1
      }else if(chars[i] == "," && level == 1){
        blocks[[length(blocks) + 1]] = gsub(",$", "", acc)
        acc = ""
      }
    }
    if(acc != "") blocks[[length(blocks) + 1]] = gsub("\\)?$", "", acc)
    blockSubs = lapply(blocks, parseSubs)
#print(data.frame(from = str, to = unlist(blocks), stringsAsFactors = FALSE))
    rbind(do.call("rbind", blockSubs), data.frame(from = str,
                                                  to = unlist(blocks),
                                                  stringsAsFactors = FALSE))
  }
}

fluxBoundsToGPR <- function(fluxBounds, softplusParam){

  toGPR <- function(vals){
    exps = exp(sapply(seq_along(vals), function(i){vals[i]}) - 1)
    exps = sapply(exps, function(x){if(is.infinite(x)) .Machine$integer.max
                                    else x})
    gprs = log(exps) - softplusParam
    sapply(gprs, function(g){if(is.infinite(g) || is.na(g)) 0.001
                            else max(g, 0.001)})
  }

  fluxBounds$absubound = abs(fluxBounds$ubound)
  fluxBounds$abslbound = abs(fluxBounds$lbound)
  fluxBounds$minbound = sapply(seq_along(fluxBounds$abslbound), function(i){
    #If bounds of different sign, then min is 0, since interval is bigger than
    # just abs value of max
    if(fluxBounds$ubound[i] != 0 &&
       (fluxBounds$lbound[i] / fluxBounds$ubound[i]) < 0) 0
    else min(fluxBounds$abslbound[i], fluxBounds$absubound[i])
  })
  fluxBounds$maxbound = sapply(seq_along(fluxBounds$abslbound), function(i){
    max(fluxBounds$abslbound[i], fluxBounds$absubound[i])
  })

  data.frame(name = fluxBounds$name, lbound = toGPR(fluxBounds$minbound),
             ubound = toGPR(fluxBounds$maxbound),
             pointEstimate = toGPR(abs(fluxBounds$pointEstimate)))
}

fluxCurvesToGPRCurves <- function(fluxCurves, softplusParam){
  times = sort(unique(fluxCurves$time))
  do.call("rbind", lapply(seq_along(times), function(t){
    timeFluxes = fluxCurves[which(fluxCurves$time == times[t]),]
    cbind(R.cache::evalWithMemoization(fluxBoundsToGPR(timeFluxes,
                                                       softplusParam),
                              key = list(softplus = softplusParam,
                                         fluxes = digest::digest(timeFluxes))),
          "time" = times[t])
  }))
}

rulesToGraph <- function(minMaxRules){
  completeEdgelist = as.matrix(do.call("rbind", lapply(minMaxRules, parseSubs)))
  res = igraph::graph_from_edgelist(completeEdgelist, directed = TRUE)
#Check whether singletons from gpr expressions are missing in the built edgelist
  remainingSingletons = unique(minMaxRules[-which(minMaxRules %in%
                                                      igraph::V(res)$name |
                                                      minMaxRules == "")])
  #Add these missing singletons as new nodes without edges
  res = res + igraph::vertices(remainingSingletons)
  res
}

setInitialBounds <- function(G, exprs, gprs){
  igraph::V(G)$lbound = sapply(igraph::V(G)$name, function(n){
    i = which(exprs == n)
    if(length(i) == 0) 0
    else min(gprs$lbound[i])
  })

  igraph::V(G)$ubound = sapply(igraph::V(G)$name, function(n){
    i = which(exprs == n)
    if(length(i) == 0) 1000
    else max(gprs$ubound[i])
  })

  igraph::V(G)$label = paste(paste("[", paste(paste(igraph::V(G)$lbound, ";"),
                                              igraph::V(G)$ubound)), "]")

  G
}

propagateBoundsDown <- function(G){
  #Propagate max / min equality as children constraints
  for(i in seq_along(igraph::V(G)$name)){
    if(grepl("^sum\\(", igraph::V(G)$name[i])){
      children = igraph::neighbors(G, igraph::V(G)[i], mode = "out")
      igraph::V(G)$ubound[children] = apply(cbind(igraph::V(G)$ubound[children],
                                          igraph::V(G)$ubound[i]), 1, min)
    }else if(grepl("^min\\(", igraph::V(G)$name[i])){
      children = igraph::neighbors(G, igraph::V(G)[i], mode = "out")
      igraph::V(G)$lbound[children] = apply(cbind(igraph::V(G)$lbound[children],
                                          igraph::V(G)$lbound[i]), 1, max)
    }
  }

  igraph::V(G)$label = paste(sapply(igraph::V(G)$name, function(n){
    if(grepl("^min", n)) "m" else if(grepl("^sum", n)) "M" else ""
  }), paste(paste("[", paste(paste(igraph::V(G)$lbound, ";"),
                             igraph::V(G)$ubound)), "]"))

  G
}

propagateBoundsUp <- function(G){
  intervals = igraph::V(G)$ubound - igraph::V(G)$lbound
  #Propagate certain values to parents with following rules :
  ## If parent is min(...) = alpha and node has certain value or interval values
  # > alpha, then remove node from parent min
  ## If parent is sum(...) = alpha and node has certain value or interval values
  # < alpha, then remove node from parent max
  for(i in seq_along(igraph::V(G)$name)){
    interval = intervals[i]
    parents = igraph::neighbors(G, igraph::V(G)[i], mode = "in")
    ub = igraph::V(G)[i]$ubound
    lb = igraph::V(G)[i]$lbound
    for(p in parents){
      pub = igraph::V(G)[p]$ubound
      plb = igraph::V(G)[p]$lbound
      pint = pub - plb
      if(!is.infinite(pint) && !is.na(pint) && pint < 0.05){
        if(grepl("^sum", igraph::V(G)[p]$name) && ub < plb){
            igraph::delete.edges(G, paste0(igraph::V(G)[p]$name,"|",
                                           igraph::V(G)[i]$name))
        }else if(grepl("^min", igraph::V(G)[p]$name) && lb > pub){
            igraph::delete.edges(G, paste0(igraph::V(G)[p]$name,"|",
                                           igraph::V(G)[i]$name))
        }
      }
    }
  }

  igraph::V(G)$label = paste(sapply(igraph::V(G)$name, function(n){
    if(grepl("^min", n)) "m" else if(grepl("^sum", n)) "M" else ""
  }), paste(paste("[", paste(paste(igraph::V(G)$lbound, ";"),
                             igraph::V(G)$ubound)), "]"))

  G
}

#' ODToFluxBounds
#'
#' @param model An object of class modelOrg, the metabolic model.
#' @param odRate The values of OD measured over time
#' @param metabolites_rates A data.frame containing the extraneous metabolites,
#' their initial concentrations and their uptake rates. Columns must be named
#' "names","concentrations" and "rates".
#' @param biomass_flux_index Optional. index of the flux corresponding to the
#' biomass reaction.
#' @return Flux bounds from OD
ODToFluxBounds <- function(odRate, model, metabolites_rates = NULL,
                         biomass_flux_index = get_biomass_flux_position(model)){
  fluxPoints = if(is.null(metabolites_rates)){
      get_fba_fluxes_from_observations(model,observed_growth_rate = odRate,
                                       biomass_flux_index= biomass_flux_index)
  }else {
      get_fba_fluxes_from_observations(model,
                                        odRate,
                                        metabolites_rates = metabolites_rates,
                                        biomass_flux_index= biomass_flux_index)}
  M = get_fva_intervals_from_observations(
    model = model,
    metabolites_rates = metabolites_rates,
    observed_growth_rate = odRate,
    biomass_flux_index= biomass_flux_index
  )
  fluxBounds = data.frame(lbound = M[, "min"], ubound = M[, "max"],
                          pointEstimate = fluxPoints, stringsAsFactors = FALSE)
  cbind("name" = rownames(fluxBounds), fluxBounds, stringsAsFactors = FALSE)
}
#' ODCurveToFluxCurves
#'
#'This function takes measured ODs and turn them into a FluxCurves object to be
#'visualize using \code{visFluxCurves()}. It relies on flux variability analysis
#' to highlight the flux value interval required to meet the specified OD.
#' @param model An object of class modelOrg, the metabolic model.
#' @param ODs A vector of measured ODs
#' @param times A vector of timepoints at which the flux balance analysis
#' solution will be evaluated.
#' @param metabolites_rates A data.frame containing the extraneous metabolites,
#' their initial concentrations and their uptake rates. Columns must be named
#' "names","concentrations" and "rates".
#' @param biomass_flux_index Optional. index of the flux corresponding to the
#' biomass reaction.
#' @examples
#' data(iMM904)
#' ODs<-seq.int(0.099,1.8,length.out = 5)
#' times = seq(0.5,2,by=0.5)
#'
#' metabolites_rates <- data.frame("name"=c("D-Glucose"),
#' "concentrations"=c(16.6),"rates"=c(-2.81))
#'
#' ODtoflux<-ODCurveToFluxCurves(model = iMM904,
#' ODs = ODs,times = times, metabolites_rates = metabolites_rates)
#'
#' visFluxCurves(ODtoflux)
#' @seealso  visFluxCurves, ODCurveToMetabolicGeneCurves,
#' visMetabolicGeneCurves
#' @export
#' @return An object FluxCurves to visualize using the function
#' \code{visFluxCurves}
ODCurveToFluxCurves <- function(model, ODs,times, metabolites_rates = NULL,
                        biomass_flux_index = get_biomass_flux_position(model)){

    lRates <- logRates(times, ODs)
    positiveDeltaIds <- which(lRates > 0)
    positiveDeltaTimes <- times[positiveDeltaIds]
    positiveDeltaCurvePoints <- ODs[positiveDeltaIds]
    odRates <- lRates[positiveDeltaIds]
  do.call("rbind", lapply(seq_along(odRates), function(i){
    metab = if(is.null(metabolites_rates)) NULL
    else {
      #metabIds = (max(i - 2, 1)) * 2 + 1
      #metabIds = c(metabIds, metabIds + 1)

        metab<-metabolites_rates#[metabIds,]
    }
    cbind(R.cache::evalWithMemoization(ODToFluxBounds(odRates[i], model, metab,
                                    biomass_flux_index = biomass_flux_index),
                              key = list(model = digest::digest(model),
                                         od = odRates[i],
                                         metabolites_rates =
                                             digest::digest(metab))),
                                         "time" = i)
  }))
}


ODToMetabolicGenes <- function(model, odRate, gprsGraph, exprs,
                               metabolites_rates, softplusParam,
                               singlePointFluxEstimate = FALSE,
                        biomass_flux_index = get_biomass_flux_position(model),
                        aliases=NULL){
  cat("odRate: ", odRate,", ", file = stderr())
  #OD rate -> flux bounds
  fluxBounds = ODToFluxBounds(odRate = odRate, model =  model,
                              metabolites_rates = metabolites_rates,
                              biomass_flux_index = biomass_flux_index)
  #flux bounds -> single point estimate (mean) -> GPR_i values
  if(singlePointFluxEstimate){
    fluxBounds$lbound = fluxBounds$pointEstimate
    fluxBounds$ubound = fluxBounds$pointEstimate
  }
  gprs = fluxBoundsToGPR(fluxBounds, softplusParam)
  G = setInitialBounds(gprsGraph, exprs, gprs)
  G = propagateBoundsDown(G)
  #G = propagateBoundsUp(G)
  if(is.null(aliases)){
      res = data.frame(
          name = igraph::V(G)$name,
          lbound = round(igraph::V(G)$lbound, 2),
          ubound = round(igraph::V(G)$ubound, 2),
          stringsAsFactors = FALSE
      )
  }else{
      res = data.frame(
      name = igraph::V(G)$name,
      commonName = sapply(igraph::V(G)$name, function(n){
          a = which(aliases$geneName_model == n)[1]
          aliases[a, "geneName_GRN"]
      }),
      lbound = round(igraph::V(G)$lbound, 2),
      ubound = round(igraph::V(G)$ubound, 2),
      stringsAsFactors = FALSE
  )}

  res[which(!grepl(",", res$name)),]
}

relativeCurve <- function(xs, ys){
  sapply(seq_along(xs), function(i){
    max(1, ys[i]/ys[1])
  })
}

logRates <- function(xs, ys){
  logCurve = log(relativeCurve(xs, ys))
  #sapply(seq_along(xs), function(i){
  #  log(max(1, ys[i]/ys[1]))
  #})
  sapply(seq_along(logCurve), function(i){
    if(i == 1) logCurve[1]
    else (logCurve[i] - logCurve[i - 1]) / (xs[i] - xs[i - 1])
  })
}


#' ODCurveToMetabolicGeneCurves
#'
#'This function takes measured ODs and turn them into a ODcurveToMetCurve object
#'to be visualize using \code{visMetabolicGeneCurves()}. It relies on flux
#'variability analysis to highlight the flux value interval required to meet the
#' specified OD and to map it on the metabolic genes.
#'
#' @param times A vector of timepoints at which the flux balance analysis
#' solution will be evaluated.
#' @param ODs  vector of measured ODs.
#' @param metabolites_rates A data.frame containing the extraneous metabolites,
#' their initial concentrations and their uptake rates. Columns must be named
#' "names","concentrations" and "rates".
#' @param model An object of class modelOrg, the metabolic model.
#' @param softplusParam  Softplus parameter identify through calibration.
#' @param singlePointFluxEstimate Optional, logical.
#' @param biomass_flux_index index of the flux corresponding to the biomass
#' reaction.
#' @param aliases Optional. A data.frame containing the gene names used in the
#' metabolic model and the aliases to use to match the regulatory network.
#' @examples ODs<-c(0.4500000,0.5322392,0.6295079,0.7445529)
#' data("aliases_SC","iMM904")
#' ODcurveToMetCurve<-ODCurveToMetabolicGeneCurves(times = seq(0.5,2,by=0.5),
#' ODs = ODs,model = iMM904,aliases = aliases_SC)
#' visMetabolicGeneCurves(ODcurveToMetCurve,genes="YJR077C")
#' @export
#' @return Metabolic genes curves to visualize using the function
#' \code{visMetabolicGeneCurves}

ODCurveToMetabolicGeneCurves <- function(times,
                                        ODs,
                                        metabolites_rates = NULL,
                                        model,
                                        softplusParam = 0,
                                        singlePointFluxEstimate = FALSE,
                        biomass_flux_index = get_biomass_flux_position(model),
                        aliases = NULL){
  lRates <- logRates(times, ODs)
  exprs <- gpr_expression_tmp(gpr(model))
  gprsGraph <- rulesToGraph(exprs)
  positiveDeltaIds <- which(lRates > 0)
  positiveDeltaTimes <- times[positiveDeltaIds]
  positiveDeltaCurvePoints <- ODs[positiveDeltaIds]
  lrates <- lRates[positiveDeltaIds]
  L <- lapply(seq_along(lrates), function(i){
    od <- lrates[i]
    metab <- if(is.null(metabolites_rates)) NULL
     else {metab <- metabolites_rates}
    R.cache::evalWithMemoization(
      ODToMetabolicGenes(model,
                         od,
                         gprsGraph,
                         exprs,
                         metabolites_rates,
                         softplusParam = softplusParam,
                         singlePointFluxEstimate,
                         biomass_flux_index=biomass_flux_index,
                         aliases=aliases),
      key <- list("od" = round(od, 3), "model" = digest::digest(model),
                 "metab" = digest::digest(metab), "softplus" = softplusParam,
                 "pointEstimate" = singlePointFluxEstimate)#,
      #force <- FALSE
    )
  })
  L <- lapply(seq_along(L), function(i){
    L[[i]]$time <- positiveDeltaTimes[i]; L[[i]]
  })
  cat(class(L), file = stderr())
  resL <- do.call("rbind", L)
  resL$pointEstimate <- (resL$ubound + resL$lbound) / 2
  resL
}

GPRToMetabolicGene <- function(gprs,
                              gprsGraph,
                              exprs,
                              singlePointEstimate = FALSE,
                              aliases = NULL){
  cat("GPRToMetabolicGene\n", file=stderr())
  #if(singlePointFluxEstimate){
  #  fluxBounds$lbound = fluxBounds$pointEstimate
  #  fluxBounds$ubound = fluxBounds$pointEstimate
  #}
  #gprs = fluxBoundsToGPR(fluxBounds, softplusParam)
  if(singlePointEstimate){
    gprs$lbound = gprs$pointEstimate
    gprs$ubound = gprs$pointEstimate
  }
  G = setInitialBounds(gprsGraph, exprs, gprs)
  G = propagateBoundsDown(G)
  if(is.null(aliases)){
      res = data.frame(
          name = igraph::V(G)$name,
          lbound = round(igraph::V(G)$lbound, 2),
          ubound = round(igraph::V(G)$ubound, 2),
          stringsAsFactors = FALSE
      )
  }else{
      res = data.frame(
          name = igraph::V(G)$name,
          commonName = sapply(igraph::V(G)$name, function(n){
              a = which(aliases$geneName_model == n)[1]
              aliases[a, "geneName_GRN"]
          }),
          lbound = round(igraph::V(G)$lbound, 2),
          ubound = round(igraph::V(G)$ubound, 2),
          stringsAsFactors = FALSE
      )}

  res[which(!grepl(",", res$name)),]
}

GPRCurvesToMetabolicGeneCurves <- function(model,
                                          gprs,
                                          singlePointEstimate = FALSE){
  exprs = gpr_expression_tmp(gpr(model))
  gprsGraph = rulesToGraph(exprs)
  times = sort(unique(gprs$time))
  resL = do.call("rbind", lapply(seq_along(times), function(t){
    timeGPRS = gprs[which(gprs$time == times[t]),]
    cbind(R.cache::evalWithMemoization(GPRToMetabolicGene(timeGPRS, gprsGraph,
                                                          exprs,
                                                 singlePointEstimate),
                              key = list(gprsGraph = digest::digest(gprsGraph),
                                         exprs = digest::digest(exprs),
                                         gprs = digest::digest(timeGPRS),
                                    singlePointEstimate = singlePointEstimate)),
          "time" = times[t])
  }))
  resL$pointEstimate = (resL$ubound + resL$lbound) / 2
  resL
}

#' Visualize Metabolic Gene Curves
#' @param metabCurves result table from ODCurveToMetabolicGeneCurves
#' @param genes a vector containing the names of the metabolic genes to plot.
#' Default select the first 50 genes
#' @param ... Optional, others curves
#' @export
#' @examples
#' data("ODcurveToMetCurve")
#'
#' visMetabolicGeneCurves(ODcurveToMetCurve,genes="YJR077C")
#' @seealso  ODCurveToMetabolicGeneCurves,ODCurveToFluxCurves, visFluxCurves
#' @return a plot of the curves of the chosen metabolic genes
visMetabolicGeneCurves <- function(metabCurves,
                                  genes = unique(metabCurves$name)[seq_len(50)],...){
  otherCurves = list(...)
  allCurves = cbind(metabCurves, "context" = 1)
  if(length(otherCurves) > 0){
    allCurves = rbind(allCurves, do.call("rbind", lapply(seq_along(otherCurves),
                                                         function(id){
      cbind(otherCurves[[id]], "context" = id + 1)
    })))
  }
  allCurves$context = as.factor(allCurves$context)
  allCurves$lbound = round(allCurves$lbound, 3)
  allCurves$ubound = round(allCurves$ubound, 3)
  allCurves$pointEstimate = round(allCurves$pointEstimate, 3)
  allCurves <- allCurves[which(allCurves$name %in% genes),]
  ggplot2::ggplot(allCurves,
                  ggplot2::aes(x = allCurves$time,fill = allCurves$context)) +
  ggplot2:: geom_ribbon(ggplot2::aes(ymin = allCurves$lbound,
                                     ymax = allCurves$ubound),
                        alpha = 0.4) +
  ggplot2::geom_line(ggplot2::aes(y = allCurves$pointEstimate,
                                  color = allCurves$context),
                     alpha = 0.8, size = 1.5) +
  ggplot2::facet_wrap(~ name, ncol = 8, scales = "free_y") +
      ggplot2::xlab("Time") + ggplot2::ylab("Bounds") +
      ggplot2::guides(fill=ggplot2::guide_legend(title="context"),
                      color = ggplot2::guide_legend(title="context"))
}

#' Visualize Fluxes Curves
#' @param fluxCurves result table from ODCurveToFluxCurves
#' @param genes a vector containing the names of the metabolic genes to plot.
#' Default select the first 50 genes
#' @param ... Optional others curves
#' @export
#' @examples
#' data("ODtoflux")
#' visFluxCurves(ODtoflux,genes ="ADK3")
#' @seealso  ODCurveToFluxCurves, ODCurveToMetabolicGeneCurves,
#' visMetabolicGeneCurves
#' @return a plot of the curves of the chosen fluxes
visFluxCurves <- function(fluxCurves,
                          genes = unique(fluxCurves$name)[seq_len(50)],...){
    otherCurves = list(...)
    allCurves = cbind(fluxCurves, "context" = 1)
    if(length(otherCurves) > 0){
        allCurves = rbind(allCurves,
                          do.call("rbind", lapply(seq_along(otherCurves),
                                                  function(id){
                              cbind(otherCurves[[id]], "context" = id + 1)})))
    }
    allCurves$context = as.factor(allCurves$context)
    allCurves$lbound = round(allCurves$lbound, 3)
    allCurves$ubound = round(allCurves$ubound, 3)
    allCurves$pointEstimate = round(allCurves$pointEstimate, 3)
    allCurves <- allCurves[which(allCurves$name %in% genes),]
    ggplot2::ggplot(allCurves,
                    ggplot2::aes(x = allCurves$time,fill = allCurves$context)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = allCurves$lbound,
                                          ymax = allCurves$ubound),
                             alpha = 0.4) +
        ggplot2::geom_line(ggplot2::aes(y = allCurves$pointEstimate,
                                        color = allCurves$context),
                           alpha = 0.8, size = 1.5)+
        ggplot2::facet_wrap(~ name, ncol = 8, scales = "free_y") +
        ggplot2::xlab("Time") + ggplot2::ylab("Bounds") +
        ggplot2::guides(fill=ggplot2::guide_legend(title="context"),
                        color = ggplot2::guide_legend(title="context"))
}
