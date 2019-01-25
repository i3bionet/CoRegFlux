testthat::context("DFBAsimulationFunctions")

library("CoRegNet")
library("sybil")
library("testthat")
library("CoRegFlux")

data("iMM904")
data("SC_experiment_influence")
data("SC_GRN_1")
data("SC_EXP_DATA")
data("aliases_SC")
data("PredictedGeneState")

metabolites<-data.frame("names" = c("D-Glucose","Ethanol"),
                        "concentrations" = c(16.6,0))

testthat::test_that("Simulation without CoRegNet finishes ", {
    expect_type( Simulation(model = iMM904,
                            metabolites = metabolites,
                            time = 1:5,
                            initial_biomass  = 0.4),type = "list" )

})

testthat::test_that("Simulation without CoRegNet throw errors for missing
                    arguments ", {
    expect_error( Simulation(model = iMM904,
                             metabolites = metabolites,
                             time = 1:5))
    expect_error( Simulation(initial_biomass  = 0.4,
                             metabolites = metabolites,
                             time = 1:5))
    expect_error( Simulation(model = iMM904,
                             time = 1:5,
                             initial_biomass  = 0.4))

})

 testthat::test_that("uptake fluxes diminish metabolite concentrations",{
     expect_lt(euler_step_metabolites(met_concentrations_t0 = 10,
                                      fluxes = -1,rate = 0.1,
                                      time_step = 0.1,
                                      biomass_t0 = 0.3),expected = 10)
 })


 testthat::test_that("positive fluxes increase metabolite concentrations",{
     expect_gt(euler_step_metabolites(met_concentrations_t0 = 10,
                                      fluxes = 1,
                                      rate = 0.1,
                                      time_step = 0.1,
                                      biomass_t0 = 0.3),expected = 10)
 })

 testthat::test_that("update_uptake_fluxes_constraints_metabolites check that
 bounds are changed",{

               model <- update_uptake_fluxes_constraints_metabolites(
                   model = iMM904,
                   met_fluxes_indexes = 550,
                   biomass_t0 = 0.3,
                   met_concentrations_t0 = 1e-6,
                   time_step = 1)
               expect_gt(lowbnd(model)[550],expected=-10)

               }
 )

testthat::test_that("coregflux_static: check bounds are changed",{
     model_gene_constraints <- coregflux_static(model= iMM904,
                                                predicted_gene_expression =
                                                    PredictedGeneState,
                                                aliases = aliases_SC)$model
     expect_false(all(model_gene_constraints@lowbnd == iMM904@lowbnd))
     expect_false(all(model_gene_constraints@uppbnd == iMM904@uppbnd))
})

regulator_table <- data.frame("regulator" = "MET32",
                              "influence" =  -1.20322,
                              "expression" = 0,
                              stringsAsFactors = FALSE)

testthat::test_that("model_TF_KO_OV_constraints: check bounds are changed for
the regulator targets",{
              model_TF_KO_OV_constraints <-
                  update_fluxes_constraints_influence(model= iMM904,
                                                      coregnet = SC_GRN_1,
                                                      regulator_table =
                                                          regulator_table,
                                                      aliases = aliases_SC )
  expect_false(all(model_TF_KO_OV_constraints@lowbnd == iMM904@lowbnd))
  expect_false(all(model_TF_KO_OV_constraints@uppbnd == iMM904@uppbnd))
  expect_length(which(model_TF_KO_OV_constraints@lowbnd != iMM904@lowbnd), 10)

              mygpr <- gpr(iMM904)[which(model_TF_KO_OV_constraints@lowbnd !=
                                             iMM904@lowbnd)]

              target_MET32 <- aliases_SC[which(aliases_SC$geneName_GRN %in%
                                                   targets(SC_GRN_1,"MET32")),1]
              expect_true(all(grepl(paste(target_MET32,collapse="|"),mygpr)))
              })

gene_table <- data.frame("gene" = c("YGL202W","YIL162W","YHR128W","YOR278W"),
                         "expression" =c(2,0,2,1),
                         stringsAsFactors = FALSE)

testthat::test_that("model_gene_KO_OV_constraints: check bounds are changed for
the given genes",{
             model_gene_KO_OV_constraints <- update_fluxes_constraints_geneKOOV(
                  model= iMM904, gene_table =  gene_table, aliases = aliases_SC)
              expect_true(model_gene_KO_OV_constraints@lowbnd[1263] ==
                              model_gene_KO_OV_constraints@lowbnd[1263])
              expect_true(model_gene_KO_OV_constraints@uppbnd[1460] !=
                              iMM904@uppbnd[1460])
              expect_equal(model_gene_KO_OV_constraints@uppbnd[1460],0)
              expect_true(model_gene_KO_OV_constraints@lowbnd[1460] ==
                              model_gene_KO_OV_constraints@lowbnd[1460])
              expect_true(model_gene_KO_OV_constraints@uppbnd[1547] !=
                              iMM904@uppbnd[1547])
              expect_equal(model_gene_KO_OV_constraints@uppbnd[1547], 2000)
              expect_true(model_gene_KO_OV_constraints@uppbnd[1545] ==
                              iMM904@uppbnd[1545])

          })
