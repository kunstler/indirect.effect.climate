####################################
####################################
## Analyse of effect on climatic variables
## Defossez et al. Oikos in revision


library(R2jags)
source(file.path('R', 'fun.fit.R'))

#####
jags.dir <- 'jags.output'
dir.create(jags.dir, showWarnings = FALSE)
output.dir <- 'output'
dir.create(output.dir, showWarnings = FALSE)
figs.dir <- 'figs'
dir.create(figs.dir, showWarnings = FALSE)
data.dir <- 'data'

## read data

data.clim <- read.csv(file.path(data.dir, 'data.climate.all.csv'))

## format data
data.herb <- fun.data.herb(data.clim)

#############################
## SEND JAGS MODEL HERBS
vars.clim <- c('delta.SWC', 'delta.VPD', 'delta.Tmin')

# NUMBER OF CHAINS TO RUN
nchains <-  4
for (var.c in vars.clim){
 #### format data per climatic variables
 jags.data <- format.jags.data.clim(data.herb, var.c)

 ##################
 #### NULL MODEL
 #### Inits
 jags.inits <- lapply(1:nchains, inits.clim, n.param = 1)

 ### SEND to jags
 clim.NULL <-jags(data=jags.data,
                  inits=jags.inits,
                  model.file = file.path("jags.model.clim.NULL"),
                  parameters.to.save = c("param","tauC"),
                  n.chains = nchains,
                  n.iter = 70000,
                  working.directory =jags.dir,
                  n.burnin=20000,n.thin=50)

 ##################
 #### L MODEL
 #### Inits
 jags.inits <- lapply(1:nchains, inits.clim, n.param = 2)

 ### SEND to jags
 clim.L <-jags(data=jags.data,
                        inits=jags.inits,
                        model.file = file.path(jags.dir, "jags.model.clim.L"),
                        parameters.to.save = c("param","tauC"),
                        n.chains = nchains,
                        n.iter = 70000,
                        n.burnin=20000,n.thin=50)

 # SAVE PLOT OF RESULTS
 pdf(file.path(jags.dir,paste0('jags.herb', var.c, '.pdf')))
 plot(clim.L)
 dev.off()

 # WARNINGS OF BAD CONVERGENCE
 if(any(clim.L$BUGSoutput$summary[, 'Rhat']>1.1)) stop('badconvergence Rhat > 1.1')
 # SAVE OUTPUTS
 obj.save <- list(clim.NULL, clim.L)
 save(obj.save,
      file = file.path(output.dir,
                      paste0("herb.res.", var.c, ".Rdata")))
}


####################
## DIC TABLE
HERB.DIC.table <- matrix(NA,nrow=length(vars.clim),ncol=2)
rownames(HERB.DIC.table) <- vars.clim
colnames(HERB.DIC.table) <- c("NULL", "CANOP")
   for (i in vars.clim)
{
    load(file.path(output.dir,
                   paste0("herb.res.", i, ".Rdata")))
    HERB.DIC.table[i,1:2] <- c(obj.save[[1]]$BUGSoutput$DIC,
                               obj.save[[2]]$BUGSoutput$DIC)
}

write.csv(HERB.DIC.table,file.path(output.dir, "HERB.DIC.table.csv"))

(HERB.DIC.table) -apply(HERB.DIC.table, MARGIN = 1, min)


#####################
##### HERB INTERCEPTED LIGHT
data.herb.light <- read.csv(file.path(data.dir, 'data.herb.light.csv'))
data.herb.light$light <-  100 - 100*(data.herb.light$light.H0  / data.herb.light$light.H20)
jags.data <- format.jags.data.herb(data.herb.light, 'light')

 ##################
 #### NULL MODEL
 #### Inits
 jags.inits <- lapply(1:nchains, inits.clim, n.param = 1)

 ### SEND to jags
 clim.NULL <-jags(data=jags.data,
                  inits=jags.inits,
                  model.file = file.path("jags.model.clim.NULL"),
                  parameters.to.save = c("param","tauC"),
                  n.chains = nchains,
                  n.iter = 70000,
                  working.directory =jags.dir,
                  n.burnin=20000,n.thin=50)

 ##################
 #### L MODEL
 #### Inits
 jags.inits <- lapply(1:nchains, inits.clim, n.param = 2)

 ### SEND to jags
 clim.L <-jags(data=jags.data,
                        inits=jags.inits,
                        model.file = file.path(jags.dir, "jags.model.clim.L"),
                        parameters.to.save = c("param","tauC"),
                        n.chains = nchains,
                        n.iter = 70000,
                        n.burnin=20000,n.thin=50)

 # SAVE PLOT OF RESULTS
 pdf(file.path(jags.dir,paste0('jags.herb', 'light', '.pdf')))
 plot(clim.L)
 dev.off()

 # WARNINGS OF BAD CONVERGENCE
 if(any(clim.L$BUGSoutput$summary[, 'Rhat']>1.1)) stop('badconvergence Rhat > 1.1')
 # SAVE OUTPUTS
 obj.save <- list(clim.NULL, clim.L)
 save(obj.save,
      file = file.path(output.dir,
                      paste0("herb.res.", 'light', ".Rdata")))

### DIC
c(obj.save[[1]]$BUGSoutput$DIC,
  obj.save[[2]]$BUGSoutput$DIC)




######################
# Compute PREDICTION
# an plot


ylab.vec <- c(expression(Delta ~ "SWC"["V" ~ -~"NV"] ~ "(%)"),
              expression(Delta ~ "VPD"["V" ~ -~"NV"] ~ "(Pa)"),
              expression(Delta ~ "Tmin"["V" ~ -~"NV"] ~ "(°C)"))
names(ylab.vec) <- vars.clim
lab.vec <- c('(a)', '(b)', '(c)')
names(lab.vec) <- vars.clim
y.t.lab <- c(11.5, 200, 0.4)
names(y.t.lab) <- vars.clim

pdf(file.path(figs.dir, 'fig.effect.clim.herb.pdf'))

par(mfrow=c(2,2))
for (var.c in vars.clim){
load(file.path(output.dir,
               paste0("herb.res.", var.c, ".Rdata")))
mat.res <- obj.save[[2]]$BUGSoutput$sims.matrix
data.herb.t <- data.herb[!is.na(data.herb[, var.c]) & !is.na(data.herb[, 'transmitted_light']), ]
temp.box <- cbind(data.herb.t[[var.c]],
                  100-data.herb.t$transmitted_light)
temp.herb.box<-  temp.box[,1]
temp.herb2 <- temp.herb.box[order(temp.box[,2])]
x<-hist(temp.box[,2],plot=F,br=6)
i<-rep(x$mids,x$counts)
box.res <- boxplot(temp.herb2~i, plot = FALSE,at=x$mids,outline=F,boxwex=4,xaxt = "n")
light.m <- mean(100 - data.herb.t$transmitted_light, na.rm = TRUE)
light.sd <- sd(100 - data.herb.t$transmitted_light, na.rm = TRUE)
### light prediction vector
light.plot <- 1:100
light.plot.scale <- (light.plot - light.m)/light.sd
## prediction pour chaque mcmc
data.herbt<-  matrix(0,nrow=length(mat.res[,2]),ncol=length(light.plot))

for (j in 1:length(mat.res[,2]))
{
data.herbt[j,] <-   mat.res[j,2] + mat.res[j,3]*light.plot.scale
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
boxplot(temp.herb2~i,add=T,at=x$mids,outline=F,boxwex=4,xaxt = "n")
text(98, y.t.lab[var.c], lab.vec[var.c])

}

###light herb
load(file.path(output.dir,
               paste0("herb.res.", 'light', ".Rdata")))
mat.res <- obj.save[[2]]$BUGSoutput$sims.matrix
data.herb.t <- data.herb.light[!is.na(data.herb.light[, 'light']) &
                               !is.na(data.herb.light[, 'transmitted_light']), ]
temp.box <- cbind(data.herb.t[['light']],
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
     ylab=expression(Delta ~ "Intercepted light"["V" ~ -~"NV"] ~ "(%)"),
     main="")
lines(light.plot,rev(mean.vec.herb),col="black",lwd=5)
lines(light.plot,rev(mean.quant.herb[1,]),lty=2,lwd=1)
lines(light.plot,rev(mean.quant.herb[2,]),lty=2,lwd=1)
boxplot(temp.herb2~i,add=T,at=x$mids,outline=F,boxwex=4,xaxt = "n")
text(98, 88, '(d)')

dev.off()

