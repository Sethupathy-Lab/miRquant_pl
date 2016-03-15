genPileUp<- function(fileName) {
   library(ggplot2)

   mylines=readLines(fileName)
   p<-ggplot()
   nm<-1
   rp<-1
   for (i in 1:length(mylines)){
      vec2<-as.vector(unlist(strsplit(mylines[i],"\t")))
      nm[i]<-paste(vec2[1],sprintf("%0.4f%%",as.numeric(vec2[2])),sep=" - ")
      rp[i]<-as.numeric(vec2[2])
      values2<-as.numeric(vec2[-(1:7)])
      mshape = rep(0,length(values2))
      for (j in 3:6){
	 idx = -as.numeric(vec2[7])+as.numeric(vec2[j])
	    mshape[idx]<-1
      }
      start2<-as.numeric(vec2[7])
      maxId2<-length(values2)+start2-1
      ID2<-c(start2:maxId2)
      ypos2<-rep(nm[i],length(ID2))

      d1<-data.frame(x=ID2,y=ypos2,v=values2,s=mshape)
      p<-p+geom_point(data=d1,aes(x=x,y=y,colour=v,size=v,shape=factor(s)),na.rm=T)
   }
   plotnmo<-nm[order(rp)]
   plotnm<-plotnmo[max(1,length(rp)-100):length(rp)]
   p<-p+scale_y_discrete('',limits=plotnm)+scale_x_continuous('Position relative to tRNA start') 
   return(p)
}

