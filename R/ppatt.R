#' use ggplot to render an image of a Drosophila blastocyst
#' @param basis matrix of basis vectors for a reexpression of the spatial expression matrix
#' @param threshold numerical lower bound on basis elements to be included when convex hulls of points forming principal patterns are computed
#' @param template matrix of 0 and 1, 1 indicating which parts of image rectangle are represented in the basis
#' @param \dots not used
#' @note Sharp thresholding can lead to warnings about patterns (vectors) that have no elelements larger than the threshold.
#' @export
ggBlast = function(basis, threshold=3, template, peel=0, ...) {
  if (missing(template)) {
    warning("no template provided, using template405")
    if (!exists("template405")) data(template405)
    template = template405
    }
  npatt = ncol(basis)
  flip = function(x) x[nrow(x):1,]
  suppressMessages({
  allpol = lapply(1:npatt, function(i)
     try(getChull(t(flip(matit(basis[,i], tmpl=template))), threshold, peel=peel)))
  })
  bad = which(sapply(allpol, function(x) inherits(x, "try-error")))
  if (length(bad)>0) {
      allpol = allpol[-bad]
      warning(paste0("dropping ", length(bad), " vectors not passing threshold"))
      npatt = npatt - length(bad)
      }
  isn = which(sapply(allpol, is.null))
  if ((k<-length(isn))>0) allpol=allpol[-isn]
  npatt = npatt-k
  nrs = which(sapply(allpol, nrow)==0)
  if ((k <- length(nrs))>0) allpol = allpol[-which(nrs==0)]
  npatt = npatt-k
  for (i in 1:npatt) allpol[[i]] = cbind(allpol[[i]], grp=i)
  dd = data.frame(do.call(rbind, allpol))
  dd$grp = factor(dd$grp)
  ggplot(data=dd) + geom_polygon(aes(x=x, y=y, group=grp, colour=grp,
     fill=grp), alpha=.8)
}


ppatt = function(x = bas, thres=.33, pal=terrain.colors, n=21, ...) {
plot(big, pch=" ", ...)
for (i in 1:ncol(bas))  try({
   ch = getChull(t(matit(x[,i], tmpl=template405)), thres)
   polygon(ch, col=pal(n, alpha=.8)[i])
   text(mean(ch[,1]), mean(ch[,2]), i)
   })
}


#' Shiny app for interactive visualization of expression patterns related to cell fate map
#' @param basis The basis component of a factorization of the spatial expression matrix
#' @param init initial value for thresholding slider
#' @param template matrix of 0 and 1, 1 indicating which parts of image rectangle are represented in the basis
#' @note A spatial expression matrix has rows corresponding to a mapping of a grid, 
#' and columns corresponding to a gene; see \code{data(expressionPatterns)} 
#' for an example.  A reexpression of such a matrix using non-negative matrix
#' factorization is described in Wu et al., 
#' PNAS, 2016 \url{http://www.pnas.org/content/113/16/4290.full}.  The basis
#' matrix of this factorization yields what Wu et al. call 'principal patterns'.
#' To visualize these patterns in the context of the cell fate map, we take
#' the convex hulls of points having pattern weights exceeding a given
#' threshold.
#' @export
CFMexplorer = function(basis, init=3, template) {
 basn = as.numeric(basis)
 low = round(min(as.numeric(basis), na.rm=TRUE),2)
 hi = round(max(as.numeric(basis), na.rm=TRUE),2)
 med = round(median(as.numeric(basis), na.rm=TRUE),2)
 if (is.null(init)) init=med
 ui = fluidPage(
        sidebarPanel(
         helpText("Interactive Cell Fate Model explorer"),
         helpText(" "), 
         helpText("Inspired by ",a(href='http://www.pnas.org/content/113/16/4290.full',"Wu et al. PNAS 2016")), 
         helpText(" "), 
         helpText("Choose threshold for NMF basis element inclusion, for construction of convex hulls"),
         sliderInput("inThresh", "threshold", min=low, max=hi, value=init),
         helpText(" "),
         checkboxInput("pickLand", "Add landmarks"),
         helpText(" "),
         helpText("Select depth of convex hull peeling"),
         numericInput("peel", " ", 0L, min=0L, max=2, step=1L)
        ),  
        mainPanel(
         tabsetPanel(
          tabPanel("chulls", plotOutput("ggpatt")),
          tabPanel("classic", plotOutput("ggpatt2"))
         )   
        )   
       )   
 server = function(input, output, session) {
   im = readPNG(system.file("pngs/springmap1.png", package="fatemapper"))
   output$ggpatt = renderPlot({
           view = ggBlast(basis, threshold=input$inThresh, 
                  template=template, peel=input$peel)
           if (input$pickLand) view = view + geom_text(data=DmLandmarks(),
                 aes(x=x,y=y,label=landm))
           print(view)
   })
   output$ggpatt2 = renderPlot({
     grid.raster(im)
     })
   }
 shinyApp(ui=ui, server=server)
}
# landmarks of the 405 point template for Dm blastoderm, using notation of Hartenstein 1985
#export
#DmLandmarks = function() {
#landmarks = c(am=19, as=270, C3d=111, cl=22,
#     dEpi=177, dr=113, es=9, ms=250, ol=56, ph=49,
#     pl=73, pm=390, pNR=58, pr=364, sg=119, tr=274, vNR=263,
#     C1=84, C2=97, C3=112, T1=128, T2=143, T3=159,
#     A1=175, A5 = 271, A10=357)
# LM_mat= matit(1:405, tmpl=template405)[16:1,]
# LM_XY = getXY(t(LM_mat), -1) 
# LM_XY = data.frame(LM_XY)
# LM_XY$landm = rep(" ", 512)
# LM_XY$landm[match(landmarks, LM_XY[,3])] = names(landmarks)
# LM_XY[which(LM_XY$landm != " "),]
#}

