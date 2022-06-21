# ---- load-packages -----------------------------------------------------------
library(lavaan)
library(doParallel)
library(parallel)
library(foreach)

# ---- declare-globals ---------------------------------------------------------
orgwd <- getwd()
n <- c(200,500)
magDif <- c(0)
facLoad.e <- c(.4,.8)
facLoad.o <- c(.4,.8)
min.dif <- c(-1,0)
nCause <- c(1)
nIndicator <- c(20)
magEff <- c(.2, .4, .6)
modCount <- 100
all.perms <- expand.grid(n, magDif, facLoad.e, facLoad.o,nCause, nIndicator, magEff, min.dif, 1:modCount)
colnames(all.perms) <- c("n", "magDif", "facLoad.e","facLoad.o",
                         "nCause", "nIndicator", 
                         "magEff", "minDif", "seedVal")
all.perms.list <- split(all.perms, f=1:nrow(all.perms))

# ---- train-models ---------------------------------------------------------
all.dat <- bettermc::mclapply(all.perms.list,mc.cores = 4,
                              FUN = function(x) try(all.steps.ts.mod(n = x$n, magDif = x$magDif, facLoad.e = x$facLoad.e,facLoad.o = x$facLoad.o,
                                                              percItems = x$percItems, nCause = x$nCause, nIndicator = x$nIndicator,
                                                              magEff = x$magEff, dataDir = "./data/simulatedDataBin/", seedVal = x$seedVal, minDif=x$minDif,
                                                              mplusDat = TRUE)),
                              mc.preschedule = FALSE,
                              mc.progress=TRUE
)
saveRDS(all.dat, "./data/tsCompare.RDS")