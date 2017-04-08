
DmLandmarks = function() {
landmarks = c(am=19, as=270, C3d=111, cl=22,
     dEpi=177, dr=113, es=9, ms=250, ol=56, ph=49,
     pl=73, pm=390, pNR=58, pr=364, sg=119, tr=274, vNR=263,
     C1=84, C2=97, C3=112, T1=128, T2=143, T3=159,
     A1=175, A5 = 271, A10=357)
 LM_mat= matit(1:405, tmpl=template405)[16:1,]
 LM_XY = getXY(t(LM_mat), -1) 
 LM_XY = data.frame(LM_XY)
 LM_XY$landm = rep(" ", 512)
 LM_XY$landm[match(landmarks, LM_XY[,3])] = names(landmarks)
 LM_XY[which(LM_XY$landm != " "),]
}

