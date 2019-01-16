testthat::context("MetabolicModelFunctions")
library("CoRegNet")
library("testthat")
library("CoRegFlux")

data("iMM904")
name=c("D-Glucose","Glycerol")
concentrations=c(13.8768699776138,0.01)
metabolites<-data.frame(name,concentrations)

 testthat::test_that("Metabolites names to model names works and return a
                     data.frame",{
     expect_true(is.data.frame(
         convert_metabolites_to_model_names(metabolites = metabolites,
                                            model = iMM904)))
     expect_equal(
         dim(convert_metabolites_to_model_names(metabolites = metabolites,
                                                model = iMM904)),c(2,3))
     result <- convert_metabolites_to_model_names(metabolites,iMM904)
     expect_equal(result$model_name,expected = c("glc__D[e]","glyc[e]"))
 })

testthat::test_that("get_metabolites_exchange_fluxes returns the proper fluxes",
                    {
    result <- get_metabolites_exchange_fluxes(iMM904,metabolites)
    expect_equal(result$flux_position,expected = c(550,555))
    expect_equal(dim(result),c(2,6))
})

testthat::test_that("build_exchange_met returns fluxes",{
    result <- build_exchange_met(iMM904)
    expect_equal(dim(result),expected = c(164,2))
    expect_true(typeof(result)=="list")
})

testthat::test_that("get_biomass_flux_position works",{
    expect_equal(get_biomass_flux_position(model = iMM904),expected = 1577)
    expect_error(get_biomass_flux_position(iMM904,
                                           biomass_reaction_id = "growth"))
})

