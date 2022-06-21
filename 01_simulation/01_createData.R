library(magrittr)
library(MplusAutomation)
library(dplyr)
library(doParallel)
library(parallel)
library(progress)

# ---- declare-functions ---------------------------------------------------------
create_nu_mod <- function(modPop=modPopulation(), n=500, nReps=50, seed=16, dataDir="./data/trainedModels/") {
  if(is.null(modPop)){
    errorCondition("Please provide population model from modPopulation function")
  }
  nu_mod <- paste0(
    "MONTECARLO:\n",
    ## Begin with the names statement
    paste0("NAMES ARE Y1-Y", modPop$outInd, " X1-X", modPop$outCau,";\n" ,collapse = ""),
    ## Create our binary variables here
    paste0("GENERATE = Y1-Y", modPop$outInd, " (1); \n", collapse = ""),
    paste0("CATEGORICAL = Y1-Y", modPop$outInd, "; \n", collapse = ""),
    ## Now do the number of observations
    paste0("NOBSERVATIONS = ", n,";\n" ,collapse = ""),
    ## Number of models
    paste0("NREPS = ", nReps,";\n" ,collapse = ""),
    ## SEED
    paste0("SEED = ", seed,";\n" ,collapse = ""),
    ## Save output
    paste0("REPSAVE = All",";\n" ,collapse = ""),
    ## data dir
    paste0("SAVE = ", "MIMICrep_", "n_", n, "difMag_", modPop$outMag, "numInd_", modPop$outInd, "numCause_", modPop$outCau, "*.dat;\n", collapse = "")
  )
  
  return(nu_mod)
}

modPopulation <- function(numCaus=1, numInd=20, magDif=.3, numDif=1, facLoad = .6, magEff=.2, minDif=-2, maxDif=2){
  ## This function will return the necassary text for the MODEL POPULATION statment in mplus
  ## This first section will return all of the scaled indicator varaibles
  ## Prep all of the cause varaibles
  sec.one <- paste0(
    paste0("[X", 1:numCaus, "@0];",sep="", collapse = "\n"),
    paste("\n", collapse = "\n"),
    paste0("X", 1:numCaus, "@1;",sep="", collapse = "\n")
  )
  ## Prep all of the indicator variables
  facLoad.sqr <- facLoad^2
  resid.var <- 1 - facLoad.sqr
  sec.two <- paste0(
    paste("\n", collapse = "\n"),
    paste("\n", collapse = "\n"),
    paste0("Y", 1:numInd, "@",resid.var, ";",sep="", collapse = "\n")
  )
  ## NOW IDENTIFY ALL OF THE DIFFICULTY VALUES HERE
  dif.vals <- round(runif(n = numInd, min = minDif, max = maxDif), 3)
  sec.two.dif <- paste0(
    paste("\n", collapse = "\n"),
    paste0("[Y", 1:numInd, "$1*", dif.vals, "];", sep='', collapse = "\n")    
  )
  ## Now create our factor variable -- only 1 for the moment
  sec.thr <- paste0("\nf BY Y1-", paste0("Y", numInd,collapse = ""),"*" , facLoad, ";" ,collapse = "")
  ## Assign varaince of the factor
  sec.fou <- paste0("\n f*1;\n",collapse = "")
  ##Assign the DIF here
  sec.fiv <- paste0("Y1-", paste0("Y", numDif,collapse = ""), " ON ", "X1-", paste0("X", numCaus, "*",magDif,";",collapse = ""),"\n")
  ## Assig magnitude of effect here
  sec.six <- paste0("\n f ON ", paste0(" X1-", paste0("X", numCaus,collapse = ""),"*", magEff,";\n"))
  ## Combine all of these
  all.out <- list(outText=paste0(sec.one, sec.two, sec.two.dif,sec.thr, sec.fou, sec.fiv, sec.six), outCau = numCaus, outDif=numDif, outMag = magDif, outInd=numInd, magEff=magEff)
  return(all.out)
}


modTest <- function(modPop = modPopulation()){
  ## This function will reutnr a list of all of the unifrom DIF models to test given a set of input
  ## causal varaibles and indicator variables simulated through a population model
  ## First identify all of the indicator variables to iterate through
  indVals <- modPop$outInd
  causVals <- modPop$outCau
  ## Now prepare the model statements
  u_mod <- function(item_num){paste0(
    "MODEL: f by \n",
    paste0("V",(1:indVals)[-item_num],"*",collapse = "\n"),
    paste0("\nV",item_num," (alpha);"),
    "\n !Factor variances set to 1;
    f@1;
    !Group predicts latent factor;\n",
    paste0("f on ",paste0("V" ,seq(indVals+1, indVals+causVals, 1)), ";"),
    paste0("\nV",item_num, " on ", paste0("V" ,seq(indVals+1, indVals+causVals, 1)), ";")
  )
  }
  ## Now grab all models
  all.mods <- lapply(1:indVals, u_mod)
  ## Now return these
  return(all.mods)
}

## This function will go ahead and take the population model
## and create a model for each varaible from all simulated datasets
## This is going to be a slow process...
## I will want to figure out to make this parallel at some point
modSamp <- function(dataDir="./data/trainedModels/", modPop=modPopulation(), n=500){
  ## First load all of the data
  all.dat <- paste0(dataDir, "MIMICrep_", "n_", n, "difMag_", modPop$outMag,"magEff",modPop$magEff ,"numInd_", modPop$outInd, "numCause_", modPop$outCau, "list.dat")
  all.dat <- read.table(all.dat)
  ## Grab the models
  all.mods <- modTest(modPop)
  ## Now create all of our models
  # create an output for all files
  all.out <- list()
  for(i in 1:length(all.dat$V1)){
    ## First load the data
    in.dat <- read.table(paste(dataDir, "/", all.dat$V1[i], sep=''))
    ## Now loop through all possible models
    for(w in 1:length(all.mods)){
      ## Create the models
      mimic_mod <- mplusObject(
        autov = FALSE,
        TITLE  = "MIMIC DIF sim;",
        #VARIABLE = paste0("V",1:dim(in.dat)[2]),
        MODEL = all.mods[[w]],
        rdata=in.dat,
        usevariables = paste0("V",1:dim(in.dat)[2]),
        imputed = FALSE
      )
      ## Now train the model
      nu_model <- mplusModeler(mimic_mod,dataout="numod.txt",modelout="numod.inp",run=1)
      ## Now export everything we need
      out.vals <- nu_model$results$parameters$unstandardized
      out.vals$modelVal <- w
      all.out <- out.vals
    }  
  }
}

## Create a function which will take the list of models & the location of the MPlus list file and train all
## of the models -- this will be used within mclapply
modMCL <- function(all.mod=modTest(), dataDir="./data/trainedModels/", n=500, modPop=NULL, orgWd = orgwd){
  ## First identify the data files
  all.dat <- paste0(dataDir, "MIMICrep_", "n_", n, "difMag_", modPop$outMag, "numInd_", modPop$outInd, "numCause_", modPop$outCau, "list.dat")
  all.dat <- read.table(all.dat)
  ## Now iterate through each of these and train the models
  out.params <- list()
  iterator <- 1
  for(q in all.dat$V1){
    ## Train the models
    in.dat <- read.table(paste(dataDir, "/", all.dat$V1[iterator], sep=''))
    mimic_mod <- mplusObject(
      autov = FALSE,
      TITLE  = "MIMIC DIF sim;",
      #VARIABLE = paste0("V",1:dim(in.dat)[2]),
      MODEL = all.mod,
      VARIABLE = paste0("categorical = ",paste0("V",1:modPop$outInd,collapse = "\n"),";"),
      rdata=in.dat,
      usevariables = paste0("V",1:dim(in.dat)[2]),
      imputed = FALSE)
    ## Now train the model
    # create random file names
    tmp.val <- sample(1:10000000, size = 1)
    dataOuttmp <- paste("numod", tmp.val, ".txt", sep='')
    modelOuttmp <- paste("numod", tmp.val, ".inp", sep='')
    modelOutrestmp <- paste("numod", tmp.val, ".out", sep='')
    setwd(dataDir)
    nu_model <- mplusModeler(mimic_mod,dataout=dataOuttmp,modelout=modelOuttmp,run=1)
    setwd(orgWd)
    ## Now grab params of int
    out.params [[iterator]] <- nu_model$results
    iterator <- iterator + 1
  }
  return(out.params)
}


## Now do one with multicore
modSamp <- function(dataDir="./data/trainedModels/", modPop=NULL, n=200){
  ## First load all of the data
  all.dat <- paste0(dataDir, "MIMICrep_", "n_", n, "difMag_", modPop$outMag, "numInd_", modPop$outInd, "numCause_", modPop$outCau, "list.dat")
  all.dat <- read.table(all.dat)
  ## Grab the models
  all.mods <- modTest(modPop)
  ## Now create all of our models
  # create an output for all files
  tmp <- parallel::mclapply(X = all.mods, FUN = function(x) modMCL(all.mod = x, dataDir = dataDir, n = n, modPop = modPop))
  return(tmp)
}

## Create the sim for the first example
mod.simp <- modPopulation(numCaus = 1, numInd = 6, magDif = .2, numDif = 2)
mod.simp.2 <- create_nu_mod(modPop = mod.simp, n = 200)
test <- mplusObject(TITLE="base",
                    MONTECARLO = mod.simp.2,
                    MODELPOPULATION = mod.simp$outText)
out.test <- mplusModeler(test, dataout = "./data/trainedModels/", modelout = "./data/trainedModels/testSim.inp", run = 1)
out.test.1 <- modSamp(dataDir = "./data/trainedModels/", modPop = mod.simp, n=200)



## Now set up all of the sim params
orgwd <- getwd()
n <- c(200,500)
magDif <- c(0)
facLoad.e <- c(.4,.8)
facLoad.o <- c(.4,.8)
percItems <- c(2,4,6)
min.dif <- c(-1,0)
nCause <- c(1)
nIndicator <- c(20)
magEff <- c(.2, .4, .6)
modCount <- 100
all.perms <- expand.grid(n, magDif, facLoad.e, facLoad.o,percItems,nCause, nIndicator, magEff, min.dif)
colnames(all.perms) <- c("n", "magDif", "facLoad.e","facLoad.o",
                         "percItems", "nCause", "nIndicator", 
                         "magEff", "minDif")

## Create the simulations datasets here
## Also provide a percent bar for this
pb <- txtProgressBar(min=0, max = dim(all.perms)[1], style = 3)
for(i in 1:nrow(all.perms)){
  ## identify our variables
  n <- all.perms[i,"n"]
  magDif <- all.perms[i,"magDif"]
  facLoad.e <- all.perms[i,"facLoad.e"]
  facLoad.o <- all.perms[i,"facLoad.o"]
  percItems <- all.perms[i,"percItems"]
  nCause <- all.perms[i,"nCause"]
  nIndicator <- all.perms[i,"nIndicator"]
  magEff <- all.perms[i, "magEff"]
  minDif <- all.perms[i,"minDif"]
  facLoad <- rep(c(facLoad.o, facLoad.e), nIndicator/2)
  ## Now create our simulation
  # First create the output directory
  output.dir <- paste("./data/simulatedDataBin/n_", n, "_magDif_", magDif, 
                      "_facLoad_e_", facLoad.e,"_facLoad_o_",facLoad.o,"_percItems_", percItems,
                      "_nIndicator_", nIndicator,"_nCause_",nCause,
                      "_magEff_", magEff,"_minDif_", minDif,
                      "/", sep='')
  if(!dir.exists(output.dir)){dir.create(output.dir)}
  ## Now create the model population values
  mod.simp <- modPopulation(numCaus = nCause, numInd = nIndicator, magDif = magDif, numDif = percItems, magEff = magEff, facLoad = facLoad, minDif = minDif, maxDif = minDif + 2)
  mod.simp.2 <- create_nu_mod(modPop = mod.simp, n = n, dataDir = output.dir, seed = i, nReps = modCount)
  test <- mplusObject(TITLE="base",
                      MONTECARLO = mod.simp.2,
                      MODELPOPULATION = mod.simp$outText)
  ## Now see if the files have already been simulated - if not run the model
  if(!file.exists(paste(output.dir, "*list.dat", sep=''))){out.test <- mplusModeler(test, dataout = output.dir, modelout = paste(output.dir, "simVals.inp", sep='/'), run = 1)}
  ## Now make the directory if it doesn't exist
  setTxtProgressBar(pb, i)
}

