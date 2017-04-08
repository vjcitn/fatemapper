#' take an estimated basis vector from NMF and format to template tmpl
#' @param x numerical vector, in the drosmap application, has length 405 representing points in the blastocyst model ellipse
#' @param nr number of rows of output matrix
#' @param nc number of columns of output matrix
#' @param tmpl matrix of nr x nc that has 0 where cells do not exist; in drosmap application, sum(tmpl) == 405
#' @examples
#' require(NMF)
#' data(exNmf21)
#' data(template405)
#' bas = basis(exNmf21)
#' mm = matit(bas[,1], tmpl=template405)
#' contour(t(mm)) # first principal pattern using naive NMF
#' @export
matit = function(x, nr=16, nc=32, tmpl) {
 tmp = matrix(0, nr=nr, nc=nc)
 tmp[tmpl>0] = x
 tmp
}

#' convert matrix representing values on a grid to coordinates and values
#' @param mat data matrix
#' @param threshold remove points with values less than threshold
#' @export
getXY = function(mat, threshold) {
 myinds = expand.grid(1:nrow(mat), 1:ncol(mat))
 dm = data.matrix(myinds)
 ans = cbind(x=dm[,1], y=dm[,2], val=as.numeric(mat))
 ans[ans[,3]>threshold,]
}

getChull = function(mat, threshold, peel=0) {
 myinds = expand.grid(1:nrow(mat), 1:ncol(mat))
 dm = data.matrix(myinds)
 big = cbind(x=dm[,1], y=dm[,2], val=as.numeric(mat))
 big = big[ which(big[,"val"]>threshold), ]
 if (is.null(big) || nrow(big)<1) return(NULL)
 cbig = chull(big)
 ans = big[cbig,1:2]
 while (peel >= 1) {
    ans = ans[-cbig,]
    cbig = chull(ans)
    ans = ans[cbig,]
    peel = peel-1
    }
 rbind(ans, ans[1,])
}
