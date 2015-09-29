####################################
####################################
## Analyse of effect on climatic variables
## Defossez et al. Oikos in revision
Sys.setlocale('LC_ALL','C')

## read data
fun.canop.clim.jags <- function(var.n,
                                clim.data.name,
                                vars.clim = c('SWC_effect', 'VPD_effect',
                                    'Tmin_effect'),
                                output.dir = 'output',
                                jags.model.dir = 'jags.model'){
library(R2jags)

data.clim <- read.csv(clim.data.name)

data.canopy <- fun.data.canopy(data.clim)

var.c <- vars.clim[var.n]


# NUMBER OF CHAINS TO RUN
nchains <-  4
 #### format data per climatic variables
 jags.data <- format.jags.data.clim(data.canopy, var.c)

 ##################
 #### NULL MODEL
 #### Inits
 jags.inits <- lapply(1:nchains, inits.clim, n.param = 1)

 ### SEND to jags
 clim.NULL <-jags(data=jags.data,
                  inits=jags.inits,
                  model.file = file.path(jags.model.dir,
                      "jags.model.clim.NULL"),
                  parameters.to.save = c("param","tauC"),
                  n.chains = nchains,
                  n.iter = 70000,
                  n.burnin=20000,n.thin=50)

 ##################
 #### L MODEL
 #### Inits
 jags.inits <- lapply(1:nchains, inits.clim, n.param = 2)

 ### SEND to jags
 clim.L <-jags(data=jags.data,
                        inits=jags.inits,
                        model.file = file.path(jags.model.dir,
                                               "jags.model.clim.L"),
                        parameters.to.save = c("param","tauC"),
                        n.chains = nchains,
                        n.iter = 70000,
                        n.burnin=20000,n.thin=50)

 # WARNINGS OF BAD CONVERGENCE
if(any(clim.L$BUGSoutput$summary[, 'Rhat']>1.1))
    stop('badconvergence Rhat > 1.1')
 # SAVE OUTPUTS
 obj.save <- list(clim.NULL, clim.L)
 save(obj.save,
      file = file.path(output.dir,
                      paste0("canop.res.", var.c, ".Rdata")))
}


####################
## ADD DIC TABLE
fun.dic.table.canop <- function(vars.clim = c('SWC_effect', 'VPD_effect',
                                    'Tmin_effect'),
                               output.dir = 'output',
                               jags.model.dir = 'jags.model'){

CANOPY.DIC.table <- matrix(NA,nrow=length(vars.clim),ncol=2)
rownames(CANOPY.DIC.table) <- vars.clim
colnames(CANOPY.DIC.table) <- c("NULL", "CANOP")
   for (i in vars.clim)
{
    load(file.path(output.dir,
                   paste0("canop.res.", i, ".Rdata")))
    CANOPY.DIC.table[i,1:2] <- c(obj.save[[1]]$BUGSoutput$DIC,
                                 obj.save[[2]]$BUGSoutput$DIC)
}

write.csv(CANOPY.DIC.table,file.path(output.dir, "CANOPY.DIC.table.csv"))

(CANOPY.DIC.table) -apply(CANOPY.DIC.table, MARGIN = 1, min)
}

fun.canop.cover.jags <- function(data.herb.light.name,
                               output.dir = 'output',
                               jags.model.dir = 'jags.model'){

#####################
##### HERB INTERCEPTED LIGHT
data.herb.light <- read.csv(data.herb.light.name)
jags.data <- format.jags.data.herb(data.herb.light, 'cover')

 ##################
 #### NULL MODEL
 #### Inits
 nchains <-  4
 jags.inits <- lapply(1:nchains, inits.clim, n.param = 1)

 ### SEND to jags
 clim.NULL <-jags(data=jags.data,
                  inits=jags.inits,
                  model.file = file.path(jags.model.dir,
                      "jags.model.clim.NULL"),
                  parameters.to.save = c("param","tauC"),
                  n.chains = nchains,
                  n.iter = 70000,
                  n.burnin=20000,n.thin=50)

 ##################
 #### L MODEL
 #### Inits
 jags.inits <- lapply(1:nchains, inits.clim, n.param = 2)

 ### SEND to jags
 clim.L <-jags(data=jags.data,
                        inits=jags.inits,
                        model.file = file.path(jags.model.dir,
                            "jags.model.clim.L"),
                        parameters.to.save = c("param","tauC"),
                        n.chains = nchains,
                        n.iter = 70000,
                        n.burnin=20000,n.thin=50)


 # WARNINGS OF BAD CONVERGENCE
 if(any(clim.L$BUGSoutput$summary[, 'Rhat']>1.1))
     stop('badconvergence Rhat > 1.1')
 # SAVE OUTPUTS
 obj.save <- list(clim.NULL, clim.L)
 save(obj.save,
      file = file.path(output.dir,
                      paste0("canop.res.", 'cover', ".Rdata")))

### DIC
c(obj.save[[1]]$BUGSoutput$DIC,
  obj.save[[2]]$BUGSoutput$DIC)
}

## PLOT
fun.predict.plot.canop <- function(clim.data.name,data.herb.light.name, vars.clim= c('SWC_effect', 'VPD_effect',
                                    'Tmin_effect'),
                               output.dir = 'output',
                               jags.model.dir = 'jags.model'){
ylab.vec <- c("Site centred SWC (%)",
              "Site centred VPD (Pa)",
              expression(paste0('Site centred Tmin ',degree,'C')))
names(ylab.vec) <- vars.clim
lab.vec <- c('(a)', '(b)', '(c)')
names(lab.vec) <- vars.clim
y.t.lab <- c(13, 260, 1.15)
names(y.t.lab) <- vars.clim

data.clim <- read.csv(clim.data.name)

data.canopy <- fun.data.canopy(data.clim)
data.herb.light <- read.csv(data.herb.light.name)


par(mfrow=c(2,2))
for (var.c in vars.clim){

load(file.path(output.dir,
               paste0("canop.res.", var.c, ".Rdata")))
mat.res <- obj.save[[2]]$BUGSoutput$sims.matrix

data.canopy.t <- data.canopy[!is.na(data.canopy[, var.c]) &
                             !is.na(data.canopy[, 'transmitted_light']), ]
temp.box <- cbind(data.canopy.t[[var.c]],
                  100-data.canopy.t$transmitted_light)
temp.herb.box<-  temp.box[,1]
temp.herb2 <- temp.herb.box[order(temp.box[,2])]
x<-hist(temp.box[,2],plot=F,br=6)
mids<-rep(x$mids,x$counts)
box.res <- boxplot(temp.herb2~mids, plot=FALSE,
                   at=x$mids,outline=F,boxwex=4,xaxt = "n")

light.m <- mean(100-data.canopy.t$transmitted_light, na.rm = TRUE)
light.sd <- sd(100-data.canopy.t$transmitted_light, na.rm = TRUE)

### light prediction vector
light.plot<- 1:100
light.plot.scale <- (light.plot - light.m)/light.sd

## prediction pour chaque mcmc
data.herbt<-  matrix(0,nrow=length(mat.res[,2]),ncol=length(light.plot))

for (i in 1:length(mat.res[,2]))
{
data.herbt[i,] <-   mat.res[i,2] + mat.res[i,3]*light.plot.scale
}


mean.vec.herb <- apply(data.herbt,MARGIN=2,FUN=mean)
mean.quant.herb <- apply(data.herbt,MARGIN=2,FUN=quantile,probs=c(0.025,0.975))

# plot


plot(light.plot,mean.vec.herb,type="l",col="black",lwd=5,
     ylim=range(box.res$stats, na.rm = TRUE),
     xlab="Canopy intercepted light (%)", ylab=ylab.vec[var.c],
     main="")

lines(light.plot,mean.vec.herb,col="black",lwd=5)
lines(light.plot,mean.quant.herb[1,],lty=2,lwd=1)
lines(light.plot,mean.quant.herb[2,],lty=2,lwd=1)

boxplot(temp.herb2~mids,add=T,at=x$mids,outline=F,boxwex=4,xaxt = "n")
text(98, y.t.lab[var.c], lab.vec[var.c])

}

###cover herb
load(file.path(output.dir,
               paste0("canop.res.", 'cover', ".Rdata")))
mat.res <- obj.save[[2]]$BUGSoutput$sims.matrix
data.herb.t <- data.herb.light[!is.na(data.herb.light[, 'cover']) &
                               !is.na(data.herb.light[, 'transmitted_light']), ]
temp.box <- cbind(data.herb.t[['cover']],
                  100-data.herb.t$transmitted_light)
temp.herb.box<-  temp.box[,1]
temp.herb2 <- temp.herb.box[order(temp.box[,2])]
x<-hist(temp.box[,2],plot=F,br=6)
i<-rep(x$mids,x$counts)
box.res <- boxplot(temp.herb2~i, plot = FALSE,at=x$mids,outline=F,boxwex=4,xaxt = "n")
### light prediction vector
light.plot <- 1:100
## prediction pour chaque mcmc
data.herbt<-  matrix(0,nrow=length(mat.res[,2]),ncol=length(light.plot))

for (j in 1:length(mat.res[,2]))
{
data.herbt[j,] <-   mat.res[j,2] + mat.res[j,3]*log(light.plot)
}

mean.vec.herb <- apply(data.herbt,MARGIN=2,FUN=mean)
mean.quant.herb <- apply(data.herbt,MARGIN=2,FUN=quantile,probs=c(0.025,0.975))

# plot
plot(light.plot,rev(mean.vec.herb),type="l",col="black",lwd=5,
     ylim=range(box.res$stats, na.rm = TRUE),
     xlab="Canopy intercepted light (%)",
     ylab=" Ground vegetation cover (%)",
     main="")
lines(light.plot,rev(mean.vec.herb),col="black",lwd=5)
lines(light.plot,rev(mean.quant.herb[1,]),lty=2,lwd=1)
lines(light.plot,rev(mean.quant.herb[2,]),lty=2,lwd=1)
boxplot(temp.herb2~i,add=T,at=x$mids,outline=F,boxwex=4,xaxt = "n")
text(98, 98, '(d)')

}






