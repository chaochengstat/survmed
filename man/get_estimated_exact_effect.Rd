% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_exact.R
\name{get_estimated_exact_effect}
\alias{get_estimated_exact_effect}
\title{Obtain the estimated mediation effect measures based on exact expressions}
\usage{
get_estimated_exact_effect(
  res,
  Mname,
  Aname,
  Wname,
  a0 = 0,
  a1 = 1,
  cw = 0,
  t = 0.5
)
}
\arguments{
\item{res}{results from the unadj.f, orc1.f, orc2.f, or rrc.f}

\item{Mname}{name of the mediator}

\item{Aname}{name of the true exposure}

\item{Wname}{names of the covariates}

\item{a0}{baseline exposure level in the mediation effect measures}

\item{a1}{active exposure level in the mediation effect measures}

\item{cw}{levels of the covariates}
}
\value{
estimation of mediation effect measures
}
\description{
Obtain the estimated mediation effect measures based on exact expressions
}
