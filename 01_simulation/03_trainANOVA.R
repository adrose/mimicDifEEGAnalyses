# ---- load-packages -----------------------------------------------------------
library(magrittr) #Pipes
library(lavaan)
library(stringr)
library(dplyr)
library(mirt)

# ---- declare-globals ---------------------------------------------------------
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
all.perms <- expand.grid(n, magDif, facLoad.e, facLoad.o,percItems,nCause, nIndicator, magEff, min.dif, 1:modCount)
colnames(all.perms) <- c("n", "magDif", "facLoad.e","facLoad.o",
                         "percItems", "nCause", "nIndicator", 
                         "magEff", "minDif", "seedVal")
all.perms.list <- split(all.perms, f=1:nrow(all.perms))

all.dat <- readRDS("./data/tsCompare.RDS")
dim.list <- lapply(all.dat, dim)
sum.list <- lapply(match.list, sum)
index <- which(sum.list!=2)
all.dat <- all.dat[-index]
all.dat <- bind_rows(all.dat)

## Now create the 
all.dat$diffVal <- 0
all.dat$Estimate[grep(x = rownames(all.dat), pattern = "mimicFS")] <- abs(all.dat$Estimate[grep(x = rownames(all.dat), pattern = "mimicFS")])
all.dat$diffVal <- all.dat$magEff - all.dat$Estimate
all.dat$diffVal[grep(x = rownames(all.dat), pattern = "mimicFS")]<- all.dat$magEff[grep(x = rownames(all.dat), pattern = "mimicFS")] - abs(all.dat$Estimate[grep(x = rownames(all.dat), pattern = "mimicFS")])
all.dat$diffValsq <- all.dat$diffVal^2
all.dat$diffValrt <- sqrt(all.dat$diffValsq)
all.dat <- all.dat[-grep(x = rownames(all.dat), pattern = "Intercept"),]

# ---- model-bias ---------------------------------------------------------
all.dat <- all.dat %>% 
  mutate(model = factor(model, levels=c("SS", "IRT", "mimic")),
         sampSize = as.factor(sampSize),
         facLoad.e = as.factor(facLoad.e),
         facLoad.o = as.factor(facLoad.o),
         minDif = as.factor(minDif),
         magEff = factor(magEff, levels=c(0.2, 0.4, 0.6)),
         facLoad = factor(paste(facLoad.e, facLoad.o),levels = c("0.4 0.4", "0.4 0.8", "0.8 0.4", "0.8 0.8")))

mod <- lm(diffVal ~ (model + sampSize + facLoad + minDif + magEff)^4, data=all.dat[-which(all.dat$model=="mimic"),])
mod2 <- lm(diffValrt ~ (model + sampSize + facLoad + minDif + magEff)^4, data=all.dat[-which(all.dat$model=="mimic"),])
mod.4.aov <- aov(diffVal ~(model + sampSize + facLoad + minDif + magEff)^4, data=all.dat[-which(all.dat$model=="mimic"),])
eta.vals <- effectsize::eta_squared(mod.4.aov, partial = TRUE)
eta.vals <- data.frame(eta.vals)
cof.vals <- effectsize::cohens_f(mod.4.aov, partial = TRUE)
cof.vals <- data.frame(cof.vals)
val.4 <- car::Anova(mod2)