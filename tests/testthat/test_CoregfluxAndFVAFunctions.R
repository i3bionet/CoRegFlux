testthat::context("CoregfluxAndFVAFonctions")

library("CoRegNet")
library("sybil")
library("testthat")
library("CoRegFlux")
data("iMM904")
data("SC_experiment_influence")
data("SC_GRN_1")
data("SC_EXP_DATA")
data("aliases_SC")

testthat::test_that("predict_linear_model_influence predicts gene expression",{
        predicted <- predict_linear_model_influence(model = iMM904,
                                                    min_Target = 4,
                                                    aliases = aliases_SC,
                                                train_expression = SC_EXP_DATA,
                                                experiment_influence =
                                                SC_experiment_influence,
                                                network = SC_GRN_1)
        expect_true(is.numeric(predicted))
        expect_true(!all(is.na(predicted)))
    })

testthat::test_that("coregflux static returns a FBA model for the test data
                    set",{
        predicted <- CoRegFlux::predict_linear_model_influence(model = iMM904,
                                                min_Target = 4,
                                                aliases = aliases_SC,
                                                train_expression = SC_EXP_DATA,
                                                experiment_influence =
                                                        SC_experiment_influence,
                                                network = SC_GRN_1)

        CoRegFlux_model<-suppressWarnings(coregflux_static(model =iMM904,
                                          gene_parameter = 0,
                                          predicted_gene_expression = predicted
        ))

        ts<-optimizeProb(CoRegFlux_model$model)
        expect_length(ts@lp_obj,1)
        expect_gte(ts@lp_obj,0)
        expect_equal(ts@lp_stat,5)
})

metabolites_rates <- data.frame("name"=c("D-Glucose"),
                                "concentrations"=c(16.6),"rates"=c(-2.81))

testthat::test_that("adjust_constraints_to_observed_rates: check bounds are
                    changed",{

    model <- adjust_constraints_to_observed_rates(model = iMM904,
                                                  metabolites_with_rates =
                                                      metabolites_rates)
    expect_equal(lowbnd(model)[550],expected=-2.81)
})
