
getdf = function(ind, thresh=2) {
 mm2 = matit(bas[,ind], tmpl=template405)
 hh = data.frame( getXY(t(mm2),thresh) )
 hh$ind = ind
 stat_density_2d(data=hh, aes(x=x,y=y,fill=..level.., colour=factor(ind)), geom="polygon")
}

