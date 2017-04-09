#' use a barplot to depict gene-oriented coefficients of a principal pattern
#' @param ppind principal pattern number: row index for 'sPPcoefGenes' to display
#' @param depth number of coefficients (starting at largest) to display
#' @param coefdata a data.frame with coefficients for genes
#' @export
genebar = function(ppind=1, depth=30, coefdata) {
 co = coefdata[ppind,order(coefdata[ppind,],decreasing=TRUE)[1:depth]] 
 nco = as.numeric(co)
 names(nco) = names(co)
 barplot(nco, las=2)
}
