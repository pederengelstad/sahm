###############################################################################
##
## Copyright (C) 2010-2012, USGS Fort Collins Science Center. 
## All rights reserved.
## Contact: talbertc@usgs.gov
##
## This file is part of the Software for Assisted Habitat Modeling package
## for VisTrails.
##
## "Redistribution and use in source and binary forms, with or without 
## modification, are permitted provided that the following conditions are met:
##
##  - Redistributions of source code must retain the above copyright notice, 
##    this list of conditions and the following disclaimer.
##  - Redistributions in binary form must reproduce the above copyright 
##    notice, this list of conditions and the following disclaimer in the 
##    documentation and/or other materials provided with the distribution.
##  - Neither the name of the University of Utah nor the names of its 
##    contributors may be used to endorse or promote products derived from 
##    this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
## THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
## PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
## CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
## EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
## PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
## OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
## WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
## OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
## ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
##
## Although this program has been used by the U.S. Geological Survey (USGS), 
## no warranty, expressed or implied, is made by the USGS or the 
## U.S. Government as to the accuracy and functioning of the program and 
## related program material nor shall the fact of distribution constitute 
## any such warranty, and no responsibility is assumed by the USGS 
## in connection therewith.
##
## Any use of trade, firm, or product names is for descriptive purposes only 
## and does not imply endorsement by the U.S. Government.
###############################################################################

confusion.matrix <- function(Stats, split.type){

  par(oma = c(4, 3, 5, 3), mar = c(20, 6, 5, 2))

  if(split.type=="none"){

    lo <- layout(matrix(data = c(1, 2), nrow = 1, ncol = 2), c(4.5, 1), 1)

  } else {

    lo <- layout(matrix(data = c(1, 2, 3), nrow = 1, ncol = 3), c(4.5, 4.5, 1), 1)

    if(split.type == "crossValidation"){

      a <- lapply(Stats[names(Stats) != "train"], function(lst){lst$Cmx})
      cmx <- a[[1]]
      for(i in 2:length(a)) cmx <- cmx + a[[i]]
      csv.stats <- apply(do.call("rbind", (lapply(Stats[names(Stats) != "train"], function(lst){return(c(lst$Sens, lst$Specf, lst$Pcc, lst$Kappa, lst$Tss))}))), 2, mean)
      Stats$crossValidation <- list(Cmx = cmx, Sens = csv.stats[1], Specf = csv.stats[2], Pcc = csv.stats[3], Kappa = csv.stats[4], Tss = csv.stats[5])
      Stats <- list("crossValidation" = Stats$crossValidation, "train" = Stats$train)
    
    }
  
  }
 
#zlim<-c(min(unlist(lapply(Stats,function(lst){100*lst$Cmx/sum(lst$Cmx)}))),max(unlist(lapply(Stats,function(lst){100*lst$Cmx/sum(lst$Cmx)}))))
#instead of basing the zlim on the acutal confusion matricies, base them on the maximum achievable value for a cell given the ratio of pres/abs
       extract.max<-function(lst){
         max(100*table(lst$auc.data$pres.abs)/length(lst$auc.dat$pres.abs))
         }
         options(warn=-1)
zlim=c(0,100)
        options(warn=0)
        
  for(i in length(Stats):1){
      image((1:2),c(2,4),matrix(data=c(100*Stats[[i]]$Cmx[2]/sum(Stats[[i]]$Cmx[1:2]),100*Stats[[i]]$Cmx[4]/sum(Stats[[i]]$Cmx[3:4]),
                                       100*Stats[[i]]$Cmx[1]/sum(Stats[[i]]$Cmx[1:2]),100*Stats[[i]]$Cmx[3]/sum(Stats[[i]]$Cmx[3:4])),nrow=2),
               zlim=zlim,xaxt="n",yaxt="n",xlab="",
               ylab="",main=paste("Confusion matrix for \n", names(Stats)[i], "data",sep=" "),col=heat.colors(100)[100:1],cex.lab=2,cex.main=2.5)
          mtext("Absence",side=2,at=2,cex=2,lwd=1.3)
          mtext("Presence",side=2,at=4,cex=2,lwd=1.3)
          mtext("Presence",side=1,at=1,cex=2,line=1,lwd=1.3)
          mtext("Absence",side=1,at=2,cex=2,line=1,lwd=1.3)
          text(x=c(1,1,2,2),y=c(2,4,2,4),
          labels=c(Stats[[i]]$Cmx[2],Stats[[i]]$Cmx[1],
                   Stats[[i]]$Cmx[4],
                   Stats[[i]]$Cmx[3]),cex=5)
              abline(h=3,lwd=5)
              abline(v=1.5,lwd=5)
         mtext(paste(
                 "Pct Correctly Classified : ",signif(Stats[[i]]$Pcc,digits=3),
               "\nSensitivity                      : ",signif(Stats[[i]]$Sens,digits=3),
               "\nSpecificity                      : ",signif(Stats[[i]]$Specf,digits=3),
               "\nTrue Skills Stat              : ",signif(Stats[[i]]$Tss,digits=3),
               "\nCohen's Kappa              : ",signif(Stats[[i]]$Kappa,digits=3),
              sep=""),
         side=1,line=13,cex=1.4,adj=0)
        box()
    }
  mtext("Observed",1,outer=TRUE,lwd=2,cex=2.5)
  mtext("Predicted",2,outer=TRUE,lwd=2,cex=2.5)
  
### color scale
 image(1,seq(from=zlim[1],to=zlim[2],length=50),
               matrix(data=seq(from=zlim[1],to=zlim[2],length=50), ncol=50,nrow=1),
              col=heat.colors(50)[50:1],
              xlab="",ylab="",zlim=zlim,
              xaxt="n",cex.lab=2,cex.axis=2)

}

