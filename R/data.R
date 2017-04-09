#' PP: Wu et al. PNAS 2016 principal patterns
#' @importFrom utils data
#' @format data.frame (405 x 21)
#' @note see template405
#' @source \url{http://www.pnas.org/content/113/16/4290.full} supplements
"PP"
#' dmMapTerms: a data frame with 27 landmark symbols and definitions for cell fate map, 1986
#' @format data.frame
#' @source \url{https://books.google.com/books?id=ol3\_CAAAQBAJ}
"dmMapTerms"
#' exNmf21: an NMFfit to the expressionPatterns of Dm blastoderm
#' @format \code{\link[NMF]{NMFfit-class}} instance
"exNmf21"
#' expressionPatterns: 405 x 1640 data.frame of expression patterns for 701 genes with some replication
#' @format data.frame
#' @source \url{http://www.pnas.org/content/113/16/4290.full} supplements
"expressionPatterns"
#' planiTerms: 17 landmark descriptions for the drosophila cell fate map
#' @format data.frame
#' @source \url{http://www.sdbonline.org/sites/FLY/atlas/00atlas.htm}
"planiTerms"
#' sPPcoefGenes: 21 x 701 data.frame with selected gene-specific coefficients 
#' @format data.frame
#' @source \url{http://www.pnas.org/content/113/16/4290.full} supplements
"sPPcoefGenes"
#' sPPcoefPatterns: 22 x 1640 data.frame with intercept and coefficients of NMF factorization of expressionPatterns
#' @format data.frame
#' @source \url{http://www.pnas.org/content/113/16/4290.full} supplements
"sPPcoefPatterns"
#' template405: 16 x 32 binary matrix, with zeroes outside the planimetric model (elliptical) for Drosophila blastoderm
#' @format matrix
"template405"
#' uniqueGenes: character vector of names of genes with spatially restricted expression in Drosophila
#' @format character vector
"uniqueGenes"
