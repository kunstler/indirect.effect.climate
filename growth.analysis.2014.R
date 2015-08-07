####################################
####################################
## ANALYSIS of first years
## seedling growth
## Defossez et al. Oikos

jags.dir <- 'jags.output'
dir.create(jags.dir, showWarnings = FALSE)
output.dir <- 'output'
dir.create(output.dir, showWarnings = FALSE)
figs.dir <- 'figs'
dir.create(figs.dir, showWarnings = FALSE)
data.dir <- 'data'

source(file.path('R', 'fun.fit.R'))
# Read Data

data.growth <- read.csv(file.path(data.dir, "data.seedling.growth.csv"))


#################################################
## JAGS CODE to run model for each species
library(R2jags)

# species vector
species <- c( "Pinus.uncinata", "Larix.decidua", "Abies.alba",
             "Fagus.sylvatica", "Quercus.petraea")


# number of chains
nchains <- 4
## loop over all species to fit model
for (i in species)
  {
  # data
  jags.data <- format.data.growth(i, data.growth)

  ##################
  ### Fit Null Model
  #### Inits
  jags.inits <- lapply(1:nchains, inits.growth.null, Nbloc = jags.data$Nbloc)

  ### SEND to jags
  GROWTH.NULL<-jags(data=jags.data,
                    inits=jags.inits,
                    model.file =  file.path(jags.dir, "GROWTH.MODEL.NULL"),
                    parameters.to.save = c("tauBLOC","tauPROC"),
                    n.chains = 4,
                    n.iter = 70000,
                    n.burnin=20000,
                    n.thin=50)

  ##################
  ### Fit T Model
  #### Inits
  jags.inits <- lapply(1:nchains, inits.growth.T, Nbloc = jags.data$Nbloc)

  ### SEND to jags
  GROWTH.T<-jags(data=jags.data,
                    inits=jags.inits,
                    model.file =  file.path(jags.dir, "GROWTH.MODEL.T"),
                    parameters.to.save = c("tauBLOC","tauPROC", "param"),
                    n.chains = 4,
                    n.iter = 70000,
                    n.burnin=20000,
                    n.thin=50)

  ###################
  ### Fit Inter Model
  #### Inits
  jags.inits <- lapply(1:nchains, inits.growth.inter, Nbloc = jags.data$Nbloc)

  ### SEND to jags
  GROWTH.INTER <-jags(data = jags.data,
                        inits = jags.inits,
                        model.file =  file.path(jags.dir, "GROWTH.MODEL.INTER"),
                        parameters.to.save = c("param","BLOC",
                                               "tauBLOC","tauPROC"),
                        n.chains = nchains,
                        n.iter =70000,
                        n.burnin=20000,
                        n.thin=50)
 # TODO ADD A WARNINGS OF BAD CONVERGENCE
if(any(GROWTH.INTER$BUGSoutput$summary[, 'Rhat']>1.1)) stop('bad convergence Rhat > 1.1')

    obj.save <- list(GROWTH.NULL, GROWTH.T, GROWTH.INTER)
    save(obj.save,
         file = file.path(output.dir,
                          paste0("growth.res.", i, ".Rdata")))
}

##
GROWTH.DIC.table <- matrix(NA,nrow=length(species),ncol=3)
rownames(GROWTH.DIC.table) <- species
colnames(GROWTH.DIC.table) <- c("NULL", "T", "ALL_INTER")
   for (i in species)
{
    load(file.path(output.dir,
                   paste0("growth.res.", i, ".Rdata")))
    GROWTH.DIC.table[i,1:3] <- c(obj.save[[1]]$BUGSoutput$DIC,
                                 obj.save[[2]]$BUGSoutput$DIC,
                                 obj.save[[3]]$BUGSoutput$DIC)
}

write.csv(GROWTH.DIC.table,file.path(output.dir, "GROWTH.DIC.table.csv"))


#######################################
#### Predict interaction coeffcient
 
 for (j in species)
{
    data.temp <- data.growth[data.growth$species==j,]
    Tair <- data.temp$dds

    load(file.path(output.dir,
                   paste0("growth.res.", j, ".Rdata")))

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
    save(obj.plot,file= file.path(output.dir, paste0("obj.plot.growth",j)))
}

## PREDICT FOR PLOT
species.lat <-  c("Pinus uncinata", "Larix decidua",
                  "Abies alba", "Fagus sylvatica",  "Quercus petraea")


df.tot <- do.call('rbind',
                  lapply(seq_len(length(species)),
                         format.df,
                         species,
                         species.lat,
                         "obj.plot.growth")
                  )
# get good order of param for plot
neworder <- c("Ground vegetation direct", "Canopy direct", "Canopy indirect")
df.tot$param <- factor(df.tot$param,levels=neworder)
levels(df.tot$param) <- c("Ground vegetation direct effect", "Canopy direct effect", "Canopy indirect effect")


## PLOT
library(ggplot2)


dat <-  expand.grid(temp = 600,
                    mean = 2.3,
                    species = unique(df.tot$species),
                    param = unique(df.tot$param))
dat$labs <- dat$species
dat$labs <- as.character(dat$labs)
dat$labs[dat$param != 'Ground vegetation direct effect'] <- ''


pdf(file.path(figs.dir, 'fig.effect.growth.ggplot.pdf'))
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
  coord_cartesian(ylim = c(-2.8, 2.8)) +
  xlab('Degree Day Sum (>5.5°C)') +
  ylab('Interaction coefficient estimates from seedling growth')+
  theme(strip.text.y = element_text(size = 0), strip.text.x = element_text(size = 9, face = 'bold'))
dev.off()











