# ---- declare-functions ---------------------------------------------------------
modPopulation <- function(numCaus=1, numInd=6, magDif=.3, numDif=1, facLoad = .6, magEff=.3, sub.val = .15, add.val=.15, seedVal=16){
  ## Set seed
  set.seed(seedVal)
  ## Now check for null effects
  if(magDif==0){
    add.val=0
    sub.val=0
  }
  ## This function will return the necessary text for the MODEL POPULATION statement in mplus
  ## This first section will return all of the scaled indicator variables
  ## Prep all of the indicator variables
  facLoad.sqr <- facLoad^2
  resid.var <- 1 - facLoad.sqr
  ## check to see if an integer or a character vector was provided for the indicator var 
  if(class(numInd)=="numeric"){
    sec.two <- paste0(
      paste("\n", collapse = "\n"),
      paste0("X", 1:numInd, "~~",resid.var,"*X", 1:numInd ,sep="", collapse = "\n")
    )
  }else if(class(numInd)=="character"){
    sec.two <- paste0(
      paste("\n", collapse = "\n"),
      paste0(numInd, "~~",resid.var,"*", numInd ,sep="", collapse = "\n")
    )    
  }
  ## Now create our factor variable -- only 1 for the moment
  ## Now randomly assign the factor loadings
  min.val <- facLoad-sub.val
  max.val <- facLoad +add.val
  fac.load.vec <- round(runif(n = numInd, min = min.val, max = max.val), 3)
  if(class(numInd)=="numeric"){
    sec.thr <- paste0("\nf =~", paste0(fac.load.vec,"*X", 1:numInd, collapse = "+") ,collapse = "")
  }else if(class(numInd)=="character"){
    sec.thr <- paste0("\nf =~", paste0(fac.load.vec,"*", numInd,collapse = "+") ,collapse = "")  
  }
  ## Assign variance of the factor
  sec.fou <- paste0("\nf~~1*f;\n",collapse = "")
  ##Assign the DIF here
  min.val.dif <- magDif -sub.val
  max.val.dif <- magDif +add.val
  magDif.vec <- round(runif(numDif, min = min.val.dif, max = max.val.dif), digits = 3)
  sec.fiv <- paste0(
    paste("\n", collapse = "\n"),
    paste0("Y", 1:numCaus, "~",paste0(magDif.vec,"*X", 1:numDif ,sep="", collapse = "+"), collapse = "\n")
  )
  ## Assign magnitude of ME here
  sec.six <- paste0(
    paste("\n", collapse = "\n"),
    paste0("f", 1:numCaus, "~",magEff ,"*Y", 1:numCaus ,sep="", collapse = "\n")
  )
  ## Constrain indicator variance here
  magEff.sqr <- magEff^2
  resid.var.me <- 1 - magEff.sqr
  sec.sev <- paste0(
    #paste("\n", collapse = "\n"),
    #paste0("[Y", 1:numInd, "@0];",sep="", collapse = "\n"),
    paste("\n", collapse = "\n"),
    paste0("Y", 1:numCaus, "~~",resid.var.me,"*Y", 1:numCaus ,sep="", collapse = "\n")
  )
  ## Combine all of these
  all.out <- list(outText=paste0(sec.thr,sec.fiv, sec.six), outCau = numCaus, outDif=numDif, outMag = magDif, outInd=numInd, facLoad=facLoad, facLoadVec = fac.load.vec, magDifVec=magDif.vec)
  return(all.out)
}

modBase <- function(modPop = modPopulation()){
  ## This function will reutnr a list of all of the uniform DIF models to test given a set of input
  ## causal varaibles and indicator variables simulated through a population model
  ## First identify all of the indicator variables to iterate through
  indVals <- modPop$outInd
  causVals <- modPop$outCau
  ## Now prepare the model statements
  u_mod <- function(item_num){paste0(
    paste0("f =~",paste0("X", 1:indVals,collapse = "+")),
    paste("\n", collapse = "\n"),
    paste0("f ~",paste0("Y", 1:causVals,collapse = "+")),
    paste("\n", collapse = "\n"))
  }
  ## Now grab all models
  all.mods <- u_mod(item_num = 1) # Providing item num although we don't really need it -- just lazy coding here
  ## Now return these
  return(all.mods)
}

## Declare a function which will create the sum score and IRT scores for the data -- given an input values
all.steps.ts.mod <- function(n = 200, magDif = 0, facLoad.e = .4, facLoad.o = .4, percItems = 2, nCause = 1, nIndicator=6, magEff = .2, minDif = -1, dataDir = "./data/simulatedDataBin/", seedVal = 16, pureMim = FALSE, mplusDat = TRUE){
  ## Load library(s)
  library(dplyr)
  library(lavaan)
  ## Define the output
  output.dir <- paste("./data/simulatedDataBin/n_", n, "_magDif_", magDif, 
                      "_facLoad_e_", facLoad.e,"_facLoad_o_",facLoad.o,"_percItems_", percItems,
                      "_nIndicator_", nIndicator,"_nCause_",nCause,
                      "_magEff_", magEff,"_minDif_", minDif,
                      "/", sep='')
  nDIFItems <- percItems
  facLoad <- rep(c(facLoad.o, facLoad.e), nIndicator/2)
  ## Make the directory if it does not exist
  if(!dir.exists(output.dir)){system(paste("mkdir ", output.dir))}
  ## First simulate the data -- not sure how to handle seed, not sure if worth worrying about it atm
  popMod <- modPopulation(magDif = magDif, numCaus = nCause, numInd = nIndicator, 
                          numDif = percItems, facLoad = facLoad, magEff = magEff, seedVal = seedVal)
  if(!mplusDat){
    set.seed(seedVal)
    myData <- try(simstandard::sim_standardized(m=popMod$outText, n=n, errors = FALSE, latent = FALSE), silent=FALSE)
    
    if(is.null(nrow(myData))){
      ## Now run w/o standardized because cov broke the mvrnorm call
      myData <- try(lavaan::simulateData(model = popMod$outText, sample.nobs = n, model.type = "sem", orthogonal = FALSE, 
                                         std.lv = TRUE, standardized = FALSE, seed = seedVal, return.fit = TRUE), silent=FALSE)
      
    }
  }else{
    myData <- read.table(paste0(output.dir,"MIMICrep_", "n_", n, "difMag_", popMod$outMag, "numInd_", popMod$outInd, "numCause_", popMod$outCau,seedVal,".dat",collapse = ""))
  }
  ## Create the IRT model
  ## Remove any items without any varaince
  if(sum(diag(var(myData))==0)!=0){
    ## Remove any indicators with no variance
    index <- which(diag(var(myData))==0)
    myData <- myData[,-index]
    nIndicator <- nIndicator - length(index)
  }
  ## Now recreate our pure mimic model
  popMod <- modPopulation(magDif = magDif, numCaus = nCause, numInd = nIndicator, 
                          numDif = percItems, facLoad = facLoad, magEff = magEff, seedVal = seedVal)
  mimic.model.syntax <- modBase(popMod)
  colnames(myData) <- c(paste("X", 1:nIndicator, sep=''), paste("Y", nCause, sep=''))
  ## Train models here
  mod.mim <- lavaan::sem(model=mimic.model.syntax, data = myData, std.lv=T, ordered = colnames(myData)[1:nIndicator])
  mirt.model <- mirt::mirt(data = myData[,paste("X", 1:nIndicator, sep='')], model = 1)
  myData$irt.fs <- mirt::fscores(mirt.model)
  ## create the ss values
  myData$ss.fs <- scale(rowSums(myData[,paste("X", 1:nIndicator, sep='')]))[,]
  ## Now train the two models
  irt.ts <- lm(Y1 ~ irt.fs, data=myData)
  ss.ts <- lm(Y1 ~ ss.fs, data=myData)
  
  ## Now write these coefficents out
  vals.one <- data.frame(summary(irt.ts)$coefficients)
  vals.one$model <- "IRT"
  vals.two <- data.frame(summary(ss.ts)$coefficients)
  vals.two$model <- "SS"
  vals.thr <- data.frame(lavaan::parameterestimates(mod.mim))
  vals.thr <- vals.thr[which(vals.thr$op=="~"),]
  vals.thr <- vals.thr[,c("est", "se", "z", "pvalue")]
  vals.thr$model <- "mimic"
  colnames(vals.thr) <- colnames(vals.two)
  rownames(vals.thr) <- "mimicFS"
  all.out <- rbind(vals.one, vals.two, vals.thr)
  all.out$seedVal <- seedVal
  ## Also attach all of the sim values
  all.out$sampSize <- n
  all.out$facLoad.e <- facLoad.e
  all.out$facLoad.o <- facLoad.o
  all.out$minDif <- minDif
  all.out$magEff <- magEff
  # and fin
  return(all.out)
}
