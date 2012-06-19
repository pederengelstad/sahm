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

make.auc.plot.jpg<-function(out=out){

  plotname<-paste(out$dat$bname,"_modelEvalPlot.jpg",sep="")
  calib.plot<-paste(out$dat$bname,"_CalibrationPlot.jpg",sep="")
  modelname<-toupper(out$input$model)
  input.list<-out$dat$ma
 
########################################################################
######################### Calc threshold on train split #################
 if(out$input$model.family!="poisson"){
            input.list$train$thresh<-out$dat$ma$train$thresh<- as.numeric(optimal.thresholds(data.frame(ID=1:nrow(input.list$train$dat),pres.abs=input.list$train$dat[,1],
                pred=input.list$train$pred),opt.methods=out$input$opt.methods))[2]
              if(out$dat$split.type%in%c("test","eval"))  input.list$test$thresh<-out$dat$ma$test$thresh<-input.list$train$thresh
            }
            else input.list$train$thresh=NULL

##################################################################
### Standard residual analysis plots for glm
    if(out$input$script.name%in%c("glm","mars") & out$dat$split.type!="eval"){
          jpeg(paste(out$dat$bname,"_stand.resid.plots.jpeg",sep=""),height=1000,width=1000)
          par(mfrow=c(2,2))
          if(out$input$script.name=="glm") plot(out$mods$final.mod,cex=1.5,lwd=1.5,cex.main=1.5,cex.lab=1.5)
          if(out$input$script.name=="mars") plot(out$mods$final.mod$model.glm,cex=1.5,lwd=1.5,cex.main=1.5,cex.lab=1.5)
          par(mfrow=c(1,1))
          graphics.off()
    }


################# Calculate all statistics on test\train or train\cv splits
  Stats<-lapply(input.list,calcStat,family=out$input$model.family)

 ##### lst doesn't contain the training portion of the data
   train.mask<-seq(1:length(Stats))[names(Stats)=="train"]

      ## breaking of the non-train split must be done separately because list structure is different for the test only and cv
    lst<-list()
    if(out$dat$split.type%in%c("test","eval"))
      lst$Test<-Stats[[-c(train.mask)]]
    if(out$dat$split.type=="crossValidation") lst<-Stats[-c(train.mask)]
    if(out$dat$split.type%in%c("none")) lst<-Stats

 #################################################
 ############### Confusion Matrix Plot ###########

  if(out$input$model.family!="poisson"){

   jpeg(file=paste(out$dat$bname,"confusion.matrix.jpg",sep="."),width=1000,height=1000,pointsize=13,quality=100)
    confusion.matrix(Stats,out$dat$split.type)
    graphics.off()
   }
########################## PLOTS ################################
#########  Residual surface of input data  ##########
  if(out$input$ResidMaps){
        if(out$dat$split.label!="eval"){
        residual.smooth.fct<-resid.image(calc.dev(input.list$train$dat$response, input.list$train$pred, input.list$train$weight, family=out$input$model.family)$dev.cont,input.list$train$pred,
                input.list$train$dat$response,input.list$train$XY$X,input.list$train$XY$Y,out$input$model.family,out$input$output.dir,label=out$dat$split.label,out)
            }
        else{
             residual.smooth.fct<-resid.image(calc.dev(input.list$test$dat$response, input.list$test$pred, input.list$test$weight, family=out$input$model.family)$dev.cont,input.list$test$pred,
                input.list$test$dat$response,input.list$test$XY$X,input.list$test$XY$Y,out$input$model.family,out$input$output.dir,label=out$dat$split.label,out)
             }
      }
########## AUC and Calibration plot for binomial data #######################

    if(out$input$model.family%in%c("binomial","bernoulli")){
            jpeg(file=plotname,height=1000,width=1000,pointsize=20,quality=100)
## ROC AUC plots
            TestTrainRocPlot(DATA=Stats$train$auc.data,opt.thresholds=input.list$train$thresh,add.legend=FALSE,lwd=2)
                 if(out$dat$split.type=="none") legend(x=.8,y=.15,paste("AUC=",round(Stats$train$auc.fit,digits=3),sep=""))
            if(out$dat$split.type!="none") {
            #so here we have to extract a sublist and apply a function to the sublist but if it has length 2 the structure of the list changes when the sublist is extracted
           if(out$dat$split.type%in%c("test","eval")){ TestTrainRocPlot(do.call("rbind",lapply(lst,function(lst){lst$auc.data})),add.roc=TRUE,line.type=2,color="red",add.legend=FALSE)
                legend(x=.55,y=.2,c(paste("Training Split (AUC=",round(Stats$train$auc.fit,digits=3), ")",sep=""),paste("Testing Split  (AUC=",round(Stats$test$auc.fit,digits=3), ")",sep="")),lty=2,col=c("black","red"),lwd=2)
                }
             if(out$dat$split.type=="crossValidation"){
             ROC.list<-list(predictions=lapply(lst,function(lst){lst$auc.data$pred}),labels=lapply(lst,function(lst){lst$auc.data$pres.abs}))
              pred <- prediction(ROC.list$predictions, ROC.list$labels)
              perf <- performance(pred,"tpr","fpr")
              plot(perf,col="grey82",lty=3,xlab="1-Specificity (False Positive)",ylab="Sensitivity (True Positive)",main="ROC Plot for Cross-Validation")
              plot(perf,lwd=1,avg="vertical",spread.estimate="boxplot",add=TRUE)
              TestTrainRocPlot(DATA=Stats$train$auc.data,opt.thresholds=input.list$train$thresh,add.legend=FALSE,lwd=1.5,add.roc=TRUE,line.type=1,col="red")
              points(1-Stats$train$Specf,Stats$train$Sens,pch=21,cex=2.5,bg="red")
               segments(x0=0,y0=0,x1=1,y1=1,col="blue")
              text(x=(.96-Stats$train$Specf),y=Stats$train$Sens+.03,label=round(Stats$train$thresh,digits=2))
                legend(x=.6,y=.22,c(paste("Training Split (AUC=",round(Stats$train$auc.fit,digits=3), ")",sep=""),
                     paste("Cross Validation Mean \n (AUC=",round(mean(unlist(lapply(lst,function(lst){lst$auc.fit}))),digits=3), ")",sep="")),lwd=c(4,1),lty=c(1,1),col=c("red","black"))
                }}
                graphics.off()

            #I'm pretty sure calibration plots should work for count data as well but I'm not quite ready to make a plot
            jpeg(file=calib.plot,height=1000,width=1000,pointsize=20,quality=100)
                cal.results<-switch(out$dat$split.type,
                            none = Stats$train$calibration.stats,
                             test = Stats$test$calibration.stats,
                             eval = Stats$test$calibration.stats,
                                crossValidation =  apply(do.call("rbind",lapply(lst,function(lst){lst$calibration.stats})),2,mean))
## Calibration plot
            a<-do.call("rbind",lapply(lst,function(lst){lst$auc.data}))
            if(out$input$PsdoAbs==TRUE) {
                pocplot(a$pred[a$pres.abs==1], a$pred[a$pres.abs==0], 
                title=paste("Presence Only Calibration Plot for \n",switch(out$dat$split.type,none="Training Data",test="Test Split",
                eval="Test Split",crossValidation="Cross Validation Split"),sep=""))
            } else{ 
            pacplot(a$pred,a$pres.abs,title=paste("Calibration Plot for ",switch(out$dat$split.type,none="Training Data",test="Test Split",eval="Test Split",crossValidation="Cross Validation Split"),sep=""))
                legend(x=0,y=1,c(as.expression(substitute(Int~~alpha==val, list(Int="Intercept:",val=signif(cal.results[1],digits=3)))),
                 as.expression(substitute(Slope~~beta==val, list(Slope="Slope:",val=signif(cal.results[2],digits=3)))),
                 as.expression(substitute(P(alpha==0, beta==1)==Prob,list(Prob=signif(cal.results[3],digits=3)))),
                 as.expression(substitute(P(alpha==0~a~beta==1)==Prob,list(Prob=signif(cal.results[4],digits=3),a="|"))),
                 as.expression(substitute(P(beta==1~a~alpha==0)==Prob,list(Prob=signif(cal.results[5],digits=3),a="|")))),bg="white")
             }
            dev.off()
            }
#Some residual plots for poisson data
    if(out$input$model.family%in%c("poisson")){
            jpeg(file=plotname)
            par(mfrow=c(2,2))
             plot(log(Stats$train$auc.data$pred[Stats$train$auc.data$pred!=0]),
                  (Stats$train$auc.data$pres.abs[Stats$train$auc.data$pred!=0]-Stats$train$auc.data$pred[Stats$train$auc.data$pred!=0]),
                  xlab="Predicted Values (log scale)",ylab="Residuals",main="Residuals vs Fitted",ylim=c(-3,3))
              abline(h=0,lty=2)
              panel.smooth(log(Stats$train$auc.data$pred[Stats$train$auc.data$pred!=0]),
              (Stats$train$auc.data$pres.abs[Stats$train$auc.data$pred!=0]-Stats$train$auc.data$pred[Stats$train$auc.data$pred!=0]))
               if(out$input$script.name!="rf"){
                    #this is the residual plot from glm but I don't think it will work for anything else
                    qqnorm(residuals(out$mods$final.mod),ylab="Std. deviance residuals")
                    qqline(residuals(out$mods$final.mod))
                     yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name("Std. Deviance Resid"))))
                    plot(log(Stats$train$auc.data$pred[Stats$train$auc.data$pred!=0]),
                       sqrt((abs(residuals(out$mods$final.mod,type="deviance")[Stats$train$auc.data$pred!=0]))),
                       xlab="Predicted Values (log Scale)",ylab=yl)
              }
            graphics.off()}

 ##################### CAPTURING TEXT OUTPUT #######################
    capture.output(cat("\n\n============================================================",
                        "\n\nEvaluation Statistics"),file=paste(out$dat$bname,"_output.txt",sep=""),append=TRUE)
      #this is kind of a pain but I have to keep everything in the same list format
      train.stats=list()
     if(out$dat$split.type=="none") train.stats<-Stats
      else train.stats$train=Stats[[train.mask]]

    capture.stats(train.stats,file.name=paste(out$dat$bname,"_output.txt",sep=""),label="train",family=out$input$model.family,opt.methods=out$input$opt.methods,out)
    if(out$dat$split.type!="none"){
    capture.output(cat("\n\n============================================================",
                        "\n\nEvaluation Statistics"),file=paste(out$dat$bname,"_output.txt",sep=""),append=TRUE)
        capture.stats(lst,file.name=paste(out$dat$bname,"_output.txt",sep=""),label=out$dat$split.label,family=out$input$model.family,opt.methods=out$input$opt.methods,out)
    }

        ############ getting statistics along with appropriate names into a data frame for creating the appended output
                        last.dir<-strsplit(out$input$output.dir,split="\\\\")
                        parent<-sub(paste("\\\\",last.dir[[1]][length(last.dir[[1]])],sep=""),"",out$input$output.dir)
                        
                       if(out$input$model.family%in%c("binomial","bernoulli")){
                           csv.stats<-lapply(Stats,function(lst){
                               return(c("","",lst$correlation,lst$pct.dev.exp,lst$Pcc,lst$Sens,lst$Specf))})
                            stat.names<-c("Correlation Coefficient","Percent Deviance Explained","Percent Correctly Classified","Sensitivity","Specificity")
                        } else{
                        csv.stats<-lapply(Stats,function(lst){
                            return(c("","",lst$correlation,lst$pct.dev.exp,lst$prediction.error/100))})
                                stat.names<-c("Correlation Coefficient","Percent Deviance Explained","Prediction Error")
                               }
                            csv.vect<-c(t(t(as.vector(unlist(csv.stats[train.mask])))),if(out$dat$split.type!="none") unlist(csv.stats[-c(train.mask)]))
                            csv.vect[seq(from=2,by=length(csv.vect)/length(Stats),length=length(Stats))]<-if(out$dat$split.type=="none"){"Train"}else{c("Train",names(lst))}
                           x=data.frame(cbind(rep(c("","",stat.names),times=length(Stats)),
                             csv.vect),row.names=NULL)

                        Header<-cbind(c("","Original Field Data","Field Data Template","PARC Output Folder","PARC Template","Covariate Selection Name",""),
                            c(last.dir[[1]][length(last.dir[[1]])],
                            out$dat$input$OrigFieldData,out$dat$input$FieldDataTemp,out$dat$input$ParcOutputFolder,
                            out$dat$input$ParcTemplate,ifelse(length(out$dat$input$CovSelectName)==0,"NONE",out$dat$input$CovSelectName),""))
                        assign("Evaluation.Metrics.List",list(x=x,Header=Header,Parm.Len=length(stat.names),parent=parent),envir=.GlobalEnv)
                      AppendOut(compile.out=out$input$Append.Dir,Header,x,out,Parm.Len=length(stat.names),parent=parent,split.type=out$dat$split.type)

    return(list(thresh=train.stats$train$thresh,residual.smooth.fct=residual.smooth.fct))
}

