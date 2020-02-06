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

chk.libs <- function(Model){
#Checks libraries for many functions I should probably pass just the list of libs to check but this helps me update 
#documentation on all libraries required by SAHM 
#Written by Marian Talbert 2/2012
#Updated 7/18/2018: Added PRROC to libs list for models [P. Engelstad]
#Updated 10/22/2019: Added ggplot2, dplyr to libs list for models [P. Engelstad]
#Updated 12/09/2019: Added xgboost for new model [P. Engelstad]

     if(Model=="PairsExplore") libs=list("gam")
     if(Model=="Pred.inspect") libs=list("raster","gam")
     if(Model%in%c("mars","glm","rf","gam","ann","brt","maxent","udc","xgb")) libs <- c("PresenceAbsence","rgdal","sp","survival","tools","raster","tcltk2","foreign","ade4","ROCR","ncf","splines","PRROC","gbm",'ggplot2','dplyr')
     if(Model=="udc")                libs<-as.list(c("rjson",libs))
     if(Model=="mars")               libs<-as.list(c("mda","earth","plotrix",libs))
     if(Model%in%c("glm","maxent"))  libs<-as.list(libs)
     if(Model=="rf")                 libs<-as.list(c("randomForest",libs))
     if(Model=="gam")                libs<-as.list(c("gam",libs))
     if(Model=="xgb")                libs<-as.list(c("xgboost",libs))
     if(Model=="brt")                libs<-as.list(c("lattice","gbm",libs))
    
     if(Model=="GenPsdAbs")   libs<-list("adehabitatHR","ks","raster","rgdal","sp","spatstat")
      lib.mssg <- unlist(suppressMessages(suppressWarnings(lapply(libs,require,quietly = T, warn.conflicts=F,character.only=T))))
      if(any(!lib.mssg)){
            install.packages(unlist(libs[!lib.mssg]), repos = "http://cran.r-project.org")
            lib.mssg <- unlist(suppressMessages(suppressWarnings(lapply(libs,require,quietly = T, warn.conflicts=F,character.only=T))))
            }
        if(any(!lib.mssg)) stop(paste(paste("\n\nthe following package(s) could not be loaded: ",paste(unlist(libs[!lib.mssg]),sep="")),
        "\n THIS IS OFTEN BECAUSE YOU DID NOT FOLLOW THE INSTALL INSTRUCTIONS ON PAGE 4 OF THE SAHM MANUAL\nINSTALL R IN A LOCATION WHERE YOU HAVE WRITE PERMISSION (Often NOT Program Files)!!!!",sep=""))

      }