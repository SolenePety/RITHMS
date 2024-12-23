% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formatting_data.R
\name{assign_taxa}
\alias{assign_taxa}
\alias{read_input_data}
\title{Assign all taxa to a cluster, eventually under genetic control using hclust}
\usage{
assign_taxa(founder_object, taxa_g = 0.1)

read_input_data(
  path_to_microbiome,
  path_to_pedmap,
  threshold = 0.05,
  ind_selected = NULL
)
}
\arguments{
\item{founder_object}{Output of read_input_data function}

\item{taxa_g}{Percentage of taxa under genetic control, DEFAULT = 0.1}

\item{path_to_microbiome}{String giving the path to count table file}

\item{path_to_pedmap}{String giving the path and prefix to ped and map file}

\item{threshold}{Threshold for rarefaction, DEFAULT = 0.05}

\item{ind_selected}{Vector of string values with individuals to keep, have to match rownames of count table file, DEFAULT = NULL}
}
\value{
A vector with a length matching the total number of taxa with values from 0 to the number of taxa, 0 corresponding to the non under genetic control cluster

A data.frame corresponding to (rarefied) microbiome with individuals in rows and taxa in columns. Genotypes data.frame is an attribute called "population" reachable using attr(output_name,"population")
}
\description{
Assign all taxa to a cluster, eventually under genetic control using hclust

Formatting data from file paths to base population object
}
\examples{
#exemples with dataset
}
