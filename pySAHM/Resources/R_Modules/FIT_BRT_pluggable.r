###############################################################################
##
## Copyright (C) 20010-2012, USGS Fort Collins Science Center. 
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

make.p.tif=T
make.binary.tif=T

tc=NULL
n.folds=3
alpha=1
n.trees=NULL
learning.rate = NULL
bag.fraction = 0.75
prev.stratify = TRUE
max.trees = 10000
tolerance.method = "auto"
tolerance = 0.001
seed=NULL
opt.methods=2
MESS=FALSE
# multCore=TRUE
predSelect=TRUE

# Interpret command line argurments #
# Make Function Call #

Args <- commandArgs(trailingOnly=FALSE)

    for (i in 1:length(Args)){
     if(Args[i]=="-f") ScriptPath<-Args[i+1]
     argSplit <- strsplit(Args[i], "=")
     if(argSplit[[1]][1]=="--file") ScriptPath <- argSplit[[1]][2]
     }

    for (arg in Args) {
    	argSplit <- strsplit(arg, "=")
    	argSplit[[1]][1]
    	argSplit[[1]][2]
    	if(argSplit[[1]][1]=="c") csv <- argSplit[[1]][2]
    	if(argSplit[[1]][1]=="o") output <- argSplit[[1]][2]
    	if(argSplit[[1]][1]=="rc") responseCol <- argSplit[[1]][2]
   		if(argSplit[[1]][1]=="mpt") make.p.tif <- as.logical(argSplit[[1]][2])
 			if(argSplit[[1]][1]=="mbt")  make.binary.tif <- as.logical(argSplit[[1]][2])
      if(argSplit[[1]][1]=="tc")  tc <- as.numeric(argSplit[[1]][2])
 			if(argSplit[[1]][1]=="nf")  n.folds <- as.numeric(argSplit[[1]][2])
 			if(argSplit[[1]][1]=="alp")  alpha <- as.numeric(argSplit[[1]][2])
 			if(argSplit[[1]][1]=="ntr")  n.trees <- as.numeric(argSplit[[1]][2])
      if(argSplit[[1]][1]=="lr")  learning.rate <- as.numeric(argSplit[[1]][2])
 			if(argSplit[[1]][1]=="bf")  bag.fraction <- as.numeric(argSplit[[1]][2])
 			if(argSplit[[1]][1]=="ps")  prev.stratify <- as.logical(argSplit[[1]][2])
 			if(argSplit[[1]][1]=="mt")  max.trees <- as.numeric(argSplit[[1]][2])
 			if(argSplit[[1]][1]=="om")  opt.methods <- as.numeric(argSplit[[1]][2])
 			if(argSplit[[1]][1]=="seed")  seed <- as.numeric(argSplit[[1]][2])
 		  if(argSplit[[1]][1]=="tolm")  tolerance.method <- argSplit[[1]][2]
 		  if(argSplit[[1]][1]=="tol")  tolerance <- as.numeric(argSplit[[1]][2])
 		  if(argSplit[[1]][1]=="mes")  MESS <-as.logical(argSplit[[1]][2])
	    if(argSplit[[1]][1]=="multicore")  multCore <- as.logical(argSplit[[1]][2])
   	  if(argSplit[[1]][1]=="pst")  predSelect <- as.logical(argSplit[[1]][2])
 			
    }

ScriptPath<-dirname(ScriptPath)
source(file.path(ScriptPath,"LoadRequiredCode.r"))
source(file.path(ScriptPath,"BRT.helper.fcts.r"))

 

    FitModels(ma.name=csv,
		output.dir=output,
		response.col=responseCol,
		make.p.tif=make.p.tif,make.binary.tif=make.binary.tif,
		simp.method="cross-validation",debug.mode=F,tc=tc,n.trees=n.trees, n.folds=n.folds,alpha=alpha,script.name="brt",
		learning.rate =learning.rate, bag.fraction = bag.fraction,prev.stratify = prev.stratify,max.trees = max.trees,seed=seed,
    opt.methods=opt.methods,MESS=MESS,tolerance.method = tolerance.method,tolerance=tolerance,ScriptPath=ScriptPath,multCore=multCore,predSelect=predSelect)



