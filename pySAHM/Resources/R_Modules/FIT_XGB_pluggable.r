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


 # Interpret command line argurments #
# Make Function Call #
#Set defaults for optional commands
make.p.tif = T
make.binary.tif = T
responseCurveForm = "pdf"


#xgb params
eta = 1.0
xgb.gamma = 0.0
max_depth = 6
min_child_weight = 1.0
nrounds = 10
subsample = 0.632
make.r.curves = T
seed = NULL
opt.methods = 2
MESS = FALSE
xtest = NULL
ytest = NULL
# multCore=TRUE

Args = commandArgs(trailingOnly=FALSE)

    for (i in 1:length(Args)){
     if(Args[i]=="-f") ScriptPath=Args[i+1]
     argSplit = strsplit(Args[i], "=")
     if(argSplit[[1]][1] == "--file") ScriptPath = argSplit[[1]][2]
     }

    for (arg in Args) {
    	argSplit = strsplit(arg, "=")
    	argSplit[[1]][1]
    	argSplit[[1]][2]
    	if(argSplit[[1]][1] == "c") csv = argSplit[[1]][2]
    	if(argSplit[[1]][1] == "o") output = argSplit[[1]][2]
    	if(argSplit[[1]][1] == "rc") responseCol = argSplit[[1]][2]
    	if(argSplit[[1]][1] == "mpt")  make.p.tif = as.logical(argSplit[[1]][2])
 			if(argSplit[[1]][1] == "mbt")  make.binary.tif = as.logical(argSplit[[1]][2])
      if(argSplit[[1]][1] == "eta")  eta <- as.numeric(argSplit[[1]][2])
 			if(argSplit[[1]][1] == "gamma")  xgb.gamma = as.numeric(argSplit[[1]][2])
 			if(argSplit[[1]][1] == "md")  max_depth = as.integer(argSplit[[1]][2])
 			if(argSplit[[1]][1] == "mcw")  min_child_weight = as.numeric(argSplit[[1]][2])
      if(argSplit[[1]][1] == "nrounds")  nrounds = as.integer(argSplit[[1]][2])
 			if(argSplit[[1]][1] == "samp")  subsample = as.numeric(argSplit[[1]][2])
 		  if(argSplit[[1]][1] == "curves")  make.r.curves =  as.logical(argSplit[[1]][2])
 		  if(argSplit[[1]][1] == "om")  opt.methods = as.numeric(argSplit[[1]][2])
      if(argSplit[[1]][1] == "mes")  MESS = as.logical(argSplit[[1]][2])
      if(argSplit[[1]][1] == "seed")  seed = as.numeric(argSplit[[1]][2])
      if(argSplit[[1]][1] == "multicore")  multCore = as.logical(argSplit[[1]][2])

 		  
    }

ScriptPath = dirname(ScriptPath)

source(file.path(ScriptPath,"LoadRequiredCode.r"))

FitModels(ma.name = csv,
          tif.dir = NULL,
          output.dir = output,
          response.col = responseCol,
          make.p.tif = make.p.tif,
          make.binary.tif = make.binary.tif,
          debug.mode = F,
          xtest = xtest,
          ytest = ytest,
          eta = eta,
          xgb.gamma = xgb.gamma,
          max_depth = max_depth,
          min_child_weight = min_child_weight,
          nrounds = nrounds,
          subsample = subsample,
          make.r.curves=make.r.curves,
          seed = seed,
          script.name = "xgb",
          opt.methods = opt.methods, 
          MESS = MESS,
          ScriptPath = ScriptPath,
          multCore = multCore)