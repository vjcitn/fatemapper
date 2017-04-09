#' render a view of spatial expression pattern of a gene
#' @param sym gene symbol, a column name in edata 
#' @param edata a matrix or data frame with columns representing genes and rows representing spatially organized measures of expression
#' @param tmpl template as used in \code{\link{matit}}
#' @param threshold numeric value beneath which expression is regarded as absent
#' @examples
#' data(expressionPatterns)
#' data(template405)
#' vizRestriction(expressionPatterns, sym="salm")
#' @export
vizRestriction = function(edata, sym="salm", tmpl=template405, threshold=.1) {
 mpnr = matit(edata[,sym], tmpl = tmpl)[16:1,]
 zz = getXY(t(mpnr),threshold)
# define ellipse path components
 ellm = matit(1:405, tmpl = template405)[16:1,]
 ellmXY = getXY(t(ellm), thresh=0) # use zero
 xyh = chull(ellmXY[,1], ellmXY[,2])
# render
 ggplot(data.frame(zz), aes(x=x,y=y,z=val)) + 
    geom_contour(show.legend=TRUE) + xlim(0,33) +
    geom_path(data=data.frame(ellmXY)[c(xyh,xyh[1]),], aes(x=x,y=y)) +
    xlab("<-- anterior") + ylab("dorsal -->")
}

