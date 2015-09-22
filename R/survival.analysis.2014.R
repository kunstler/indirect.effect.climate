####################################
####################################
## Analyse of first years
## seedling survival
## Defossez et al. Oikos in revision

fun.surv.jags <- function(sp.n, data.seedling.survival.name,
                          species = c( "Pinus.uncinata", "Larix.decidua",
                              "Abies.alba", "Fagus.sylvatica",
                              "Quercus.petraea"),
                          jags.model.dir = 'jags.model',
                          output.dir = 'output'){
library(R2jags)

data.survival <- read.csv(file=data.seedling.survival.name,
                          sep = ',')

### remove indiv with missing data
data.survival <- subset(data.survival, !is.na(transmitted_light))# light
data.survival <- subset(data.survival, !is.na(dds))# DDS
data.survival <- subset(data.survival, !is.na(seedlings_survival))# survival

# species selected
sp <- species[sp.n]
#################################
##### MCMC ESTIMATION WITH JAGS

# NUMBER OF CHAINS TO RUN
nchains <-  4
 #### format data per species
 jags.data <- format.data.survival(sp, data.survival)

 ##################
 #### NULL MODEL
 #### Inits
 jags.inits <- lapply(1:nchains, inits.survival.null, Nbloc = jags.data$Nbloc)

 ### SEND to jags
 Survival.null <-jags(data=jags.data,
                        inits=jags.inits,
                        model.file = file.path(jags.model.dir, "SURVIVAL.MODEL.NULL"),
                        parameters.to.save = c("BLOC","tauBLOC"),
                        n.chains = nchains,
                        n.iter = 70000,
                        n.burnin=20000,n.thin=50)

 ##################
 #### T MODEL
 #### Inits
 jags.inits <- lapply(1:nchains, inits.survival.T, Nbloc = jags.data$Nbloc)

 ### SEND to jags
 Survival.T <-jags(data=jags.data,
                        inits=jags.inits,
                        model.file = file.path(jags.model.dir, "SURVIVAL.MODEL.T"),
                        parameters.to.save = c("BLOC","tauBLOC", "param"),
                        n.chains = nchains,
                        n.iter = 70000,
                        n.burnin=20000,n.thin=50)

 ##################
 #### INTER MODEL
 #### Inits
 jags.inits <- lapply(1:nchains, inits.survival.inter, Nbloc = jags.data$Nbloc)

 ### SEND to jags
 Survival.inter <-jags(data=jags.data,
                        inits=jags.inits,
                        model.file = file.path(jags.model.dir, "SURVIVAL.MODEL.INTER"),
                        parameters.to.save = c("param","BLOC","tauBLOC"),
                        n.chains = nchains,
                        n.iter = 70000,
                        n.burnin=20000,n.thin=50)


 # WARNINGS OF BAD CONVERGENCE
if(any(Survival.inter$BUGSoutput$summary[, 'Rhat']>1.1)) stop('badconvergence Rhat > 1.1')
 # SAVE OUTPUTS
 obj.save <- list(Survival.null, Survival.T, Survival.inter)
 save(obj.save,
      file = file.path(output.dir,
                      paste0("survival.res.", sp, ".Rdata")))
}


fun.dic.table.surv <- function(species = c( "Pinus.uncinata", "Larix.decidua",
                                "Abies.alba","Fagus.sylvatica",
                                "Quercus.petraea"),
                               output.dir = 'output',
                               jags.model.dir = 'jags.model'){
## BUILD DIC TABLE
SURVIVAL.DIC.table <- matrix(NA,nrow=length(species),ncol=3)
rownames(SURVIVAL.DIC.table) <- species
colnames(SURVIVAL.DIC.table) <- c("NULL", "T", "ALL_INTER")
   for (i in species)
{
    load(file.path(output.dir,
                   paste0("survival.res.", i, ".Rdata")))
    SURVIVAL.DIC.table[i,1:3] <- c(obj.save[[1]]$BUGSoutput$DIC,
                                 obj.save[[2]]$BUGSoutput$DIC,
                                 obj.save[[3]]$BUGSoutput$DIC  )
}

write.csv(SURVIVAL.DIC.table,file.path(output.dir, "SURVIVAL.DIC.table.csv"))

(SURVIVAL.DIC.table) -apply(SURVIVAL.DIC.table, MARGIN = 1, min)
}




fun.R2.table.surv <- function(data.seedling.survival.name,
                              species = c("Pinus.uncinata", "Larix.decidua",
                                          "Abies.alba", "Fagus.sylvatica",
                                          "Quercus.petraea"),
                               output.dir = 'output',
                               jags.model.dir = 'jags.model'){
#########################
### COMPUTE R2m and R2c
data.survival <- read.csv(file=data.seedling.survival.name,
                          sep = ',')

### remove indiv with missing data
data.survival <- subset(data.survival, !is.na(transmitted_light))# light
data.survival <- subset(data.survival, !is.na(dds))# DDS
data.survival <- subset(data.survival, !is.na(seedlings_survival))# survival

SURVIVAL.R2.table <- matrix(NA,nrow=length(species),ncol=2)
rownames(SURVIVAL.R2.table) <- species
colnames(SURVIVAL.R2.table) <- c("R2m", "R2c")
   for (i in species)
{
    jags.data <- format.data.survival(i, data.survival)

    load(file.path(output.dir,
                   paste0("survival.res.", i, ".Rdata")))
    R2mc.list <- fun.R2mc.jags.full.model(obj.save[[3]],
                                          jags.data$Tair,
                                          jags.data$Lum,
                                          jags.data$Herb)
    SURVIVAL.R2.table[i,1:2] <- c(mean(R2mc.list[[1]]),
                                  mean(R2mc.list[[2]]))
}

write.csv(SURVIVAL.R2.table,file.path(output.dir, "SURVIVAL.R2.table.csv"))
}


#######################################
#### Predict interaction coeffcient

fun.surv.plot <- function( data.seedling.survival.name,
                          species = c( "Pinus.uncinata", "Larix.decidua",
                              "Abies.alba","Fagus.sylvatica",
                              "Quercus.petraea"),
                          jags.model.dir = 'jags.model',
                          output.dir = 'output'){

data.survival <- read.csv(file=data.seedling.survival.name,
                          sep = ',')

### remove indiv with missing data
data.survival <- subset(data.survival, !is.na(transmitted_light))# light
data.survival <- subset(data.survival, !is.na(dds))# DDS
data.survival <- subset(data.survival, !is.na(seedlings_survival))# survival




for (j in species)
{
    data.temp <- data.survival[data.survival$species == j,]
    Tair <- data.temp$dds

    load(file.path(output.dir,
                   paste0("survival.res.", j, ".Rdata")))

    res <- obj.save[[3]]$BUGSoutput$sims.matrix
    temp.plot<- 600:2500
    temp.plot.CR <- (temp.plot -mean(Tair))/sd(Tair)

    ## herbaceous effect
    list.veg <- inter.param.mean.ci.vs.T('param[4]', 'param[6]',
                                          res, temp.plot.CR, temp.plot)
    ## canopy effect
    list.can <- inter.param.mean.ci.vs.T('param[3]', 'param[7]',
                                         res, temp.plot.CR, temp.plot)
    ## indirect effect
    list.ind <- inter.param.mean.ci.vs.T('param[5]', 'param[8]',
                                         res, temp.plot.CR, temp.plot)
    obj.plot <-list(list.veg,
                    list.can,
                    list.ind)
    save(obj.plot,file= file.path(output.dir, paste0("obj.plot.survival",j)))
}

##
species.lat <-  c("Pinus uncinata", "Larix decidua",
                  "Abies alba", "Fagus sylvatica")

df.tot <- do.call('rbind',
                  lapply(seq_len(length(species[-5])),
                         format.df,
                         species[-5],
                         species.lat,
                         "obj.plot.survival",
                         output.dir)
                  )
# get good order of param for plot
neworder <- c("Ground vegetation direct", "Canopy direct", "Canopy indirect")
df.tot$param <- factor(df.tot$param,levels=neworder)
levels(df.tot$param) <- c("Ground vegetation direct effect",
                          "Canopy direct effect", "Canopy indirect effect")


## Figure 3: parameters of survival model f(SGDD)
library(ggplot2)

dat <-  expand.grid(temp = 600,
                    mean = 3,
                    species = unique(df.tot$species),
                    param = unique(df.tot$param))
dat$labs <- dat$species
dat$labs <- as.character(dat$labs)
dat$labs[dat$param != 'Ground vegetation direct effect'] <- ''


ggplot(df.tot, aes(x = temp, y = mean)) +
  geom_line(aes(x = temp, y = mean)) +
  geom_ribbon(aes(ymin=quant.l, ymax=quant.h, show_guide=FALSE, linetype = NA),
              alpha=0.2) +
  geom_line(aes(x = temp, y = quant.l), linetype = 3) +
  geom_line(aes(x = temp, y = quant.h), linetype = 3) +
  facet_grid(species ~ param) +  theme_simple() +
  geom_text(aes(x=temp, y=mean, label=labs, group=NULL),
            data=dat, fontface = 'bold.italic', size = 4, hjust = 0) +
  geom_hline(yintercept = 0, col = 'black',linetype=2) +
  coord_cartesian(ylim = c(-3.8, 3.8)) +
  xlab('Degree Day Sum (>5.5°C)') +
  ylab('Interaction coefficient estimates from seedling survival')+
  theme(strip.text.y = element_text(size = 0), strip.text.x = element_text(size = 9, face = 'bold'))
}

