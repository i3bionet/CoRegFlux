% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DFBAsimulationFunctions.R
\name{Simulation_Step}
\alias{Simulation_Step}
\title{Single simulation step}
\usage{
Simulation_Step(model, coregnet, metabolites, met_concentrations_t0,
  biomass_t0, regulator_table, gene_table, time_step, gene_state,
  softplus_parameter, aliases, biomass_flux_index)
}
\arguments{
\item{model}{An object of class modelOrg, the metabolic model.}

\item{coregnet}{Optional, object of class CoRegNet object containing
the information about regulatory and coregulatory relationships}

\item{metabolites}{data frame of metabolites names}

\item{met_concentrations_t0}{data frame of metabolites concentrations at t0,
before performing a time step}

\item{biomass_t0}{biomass at t0 before performing a step}

\item{regulator_table}{A data.frame containing 3 columns: "regulator",
"influence","expression" containing respectively the name of a TF present in
the CoRegNet object as string, its influence in the condition of interest as
a numerical and an expression factor of 0 for a KO, or an integer >1 for an
 overexpression}

\item{gene_table}{A data.frame containing 2 columns: "gene",
"expression" containing respectively the name of a gene present in
the CoRegNet object as string and an expression factor of 0 for a KO,
or an integer >1 for an overexpression}

\item{time_step}{size of the time step to perform; that is t1-t0}

\item{gene_state}{data frame with rows being gene names and columns being
the gene expression or any other continuous value representing metabolites
 activity to be evaluated using the gpr rules}

\item{softplus_parameter}{Softplus parameter identify through calibration.
Default to 0.}

\item{aliases}{Optional. A data.frame containing the gene names currently
used in the network under the colname "geneName" and the alias under the
colnames "alias"}

\item{biomass_flux_index}{index of the flux corresponding to the biomass
reaction.}
}
\value{
list of:
fluxes: fluxes of the resulting fba soulution to the metabolic and genetic
constraints
 biomass_yield: biomass yield that is used as proxy for the growth rate in
 the dynamic flux balance analysis solution. Corresponds to the flux of the
 biomass reaction of the model.
}
\description{
Single simulation step in which fluxes are reconstrained according to
metabolite concentrations, then given the continuous evaluation of the
 gpr rules and the softplus function of the gene regulatory state.
}
\seealso{
Simulation
}
