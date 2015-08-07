##########################
##########################
## FUNCTIONS FOR JAGS FIT

jags.dir <- 'jags.output'
dir.create(jags.dir, showWarnings = FALSE)


format.data.growth <- function(sp, DF){
  ## Select data species
  data.temp <- DF[DF$species == sp,]

  Lum <- drop(scale(data.temp$transmitted_light))
  Growth <- drop(scale(data.temp$height_growth))
  N.indiv<- length(Growth)
  Tair <- drop(scale(data.temp$dds))
  Herb <- rep(0, nrow(data.temp))
  Herb[data.temp[,'ground_vegetation_treatment']=="V"] <- 1
  num.bloc <- as.vector(unclass(factor(paste(data.temp[,'site'],
                                             data.temp[,'block']))))
  Nbloc <- length(table(num.bloc))
  jags.data<-list("Growth" = Growth,
                  "N.indiv" = N.indiv,
                  "Nbloc" = Nbloc,
                  "num.bloc" = num.bloc,
                  "Tair" = Tair,
                  "Herb" = Herb,
                  "Lum" = Lum)
return(jags.data)
}

format.data.survival <- function(sp, DF){
  ## Select data species
  data.temp <- DF[DF$species == sp,]
  Lum <- as.vector(scale(data.temp$transmitted_light))
  Survie <- data.temp$seedlings_survival
  ## create a vector of 0/1 for the survival
  N.indiv<- length(Survie)
  Tair <- as.vector(scale(data.temp$dds))
  Herb <- rep(0, nrow(data.temp))
  Herb[data.temp[,'ground_vegetation_treatment']=="V"] <- 1
  num.bloc <- as.vector(unclass(factor(paste(data.temp[,'site'],
                                             data.temp[,'block']))))
  Nbloc <- length(unique(num.bloc))

  jags.data<-list("Survie" = Survie,
                  "N.indiv" = N.indiv,
                  "Nbloc" = Nbloc,
                  "num.bloc" = num.bloc,
                  "Tair" = Tair,
                  "Herb" = Herb,
                  "Lum" = Lum)
return(jags.data)
}

###############################
## WRITE JAGS MODELS SURVIVAL
MODEL.INTER<-
 "#############################################################################
 ######################## SURVIVAL model with BUGGS ###########################
 model {
 ############ Likelihood ###################
     for (i in 1:N.indiv) {
 ## logistic regression
 Survie[i] ~ dbern(ppp[i])
 ##enviro effect on the bernouilli
 logit(ppp[i]) <- param[1] +
                param[2]*Tair[i] +
                param[3]*Lum[i] +
                param[4]*Herb[i] +
                param[5]*Lum[i]*Herb[i] +
                param[6]*Tair[i]*Herb[i] +
                param[7]*Tair[i]*Lum[i] +
                param[8]*Tair[i]*Lum[i]*Herb[i] +
                BLOC[num.bloc[i]]
}
 ####  priors
 ################################################
 ########### Hierarchical parameters ########
for (n in 1:Nbloc)
{
    BLOC[n]~dnorm(0,tauBLOC)
}
 ###############################################
 ########### Non-hierarchical parameters ########
for (j in 1:8)
{
param[j] ~dnorm(0,1.0E-6)
}
tauBLOC ~ dgamma(0.001,0.001)
 } # End of the jags model
 "

cat(MODEL.INTER,
    file = file.path(jags.dir, "SURVIVAL.MODEL.INTER"),
    sep=" ", fill = FALSE,
    labels = NULL, append = FALSE)

MODEL.NULL<-
 "
 ####### SURVIVAL model with BUGGS ######
 model {
 ####### Likelihood ######
     for (i in 1:N.indiv) {
 ## logistic regression
 Survie[i] ~ dbern(ppp[i])
 ##enviro effect on the bernouilli
 logit(ppp[i]) <- param[1] + BLOC[num.bloc[i]]
}
 #######################################
 ####  priors

 ### Hierarchical parameters ###
for (n in 1:Nbloc)
{
    BLOC[n]~dnorm(0,tauBLOC)
}
 ### Non-hierarchical parameters ###
for (j in 1:1)
{
param[j] ~dnorm(0,1.0E-6)
}
tauBLOC ~ dgamma(0.001,0.001)
 } # End of the jags model
 "

cat(MODEL.NULL,
    file = file.path(jags.dir, "SURVIVAL.MODEL.NULL"),
    sep=" ", fill = FALSE,
    labels = NULL, append = FALSE)

MODEL.T<-
 "
 #### SURVIVAL model with BUGGS ###
 model {
 ### Likelihood ###
     for (i in 1:N.indiv) {
 ## logistic regression
 Survie[i] ~ dbern(ppp[i])
 ##enviro effect on the bernouilli
 logit(ppp[i]) <- param[1] +
                param[2]*Tair[i] + BLOC[num.bloc[i]]
}
 ##################
 ####  priors

 #### Hierarchical parameters ###
for (n in 1:Nbloc)
{
    BLOC[n]~dnorm(0,tauBLOC)
}
 ### Non-hierarchical parameters ###
tauBLOC ~ dgamma(0.001,0.001)
for (j in 1:2)
{
param[j] ~dnorm(0,1.0E-6)
}

 } # End of the jags model
 "

cat(MODEL.T,
    file = file.path(jags.dir, "SURVIVAL.MODEL.T"),
    sep=" ", fill = FALSE,
    labels = NULL, append = FALSE)



##################################
##################################
### WRITE JAGS MODELS GROWTH

## NULL MODEL
NULL.MODEL <-
 "model {
    ############ Likelihood ###################
         for (i in 1:N.indiv) {
    ## linear regression
    Growth[i] ~ dnorm(theo.growth[i],tauPROC)
    ## enviro effect
        theo.growth[i] <- param[1]+BLOC[num.bloc[i]]
    }
    ########### Hierarchical parameters ########
    for (n in 1:Nbloc)
    {
        BLOC[n]~dnorm(0,tauBLOC)
    }
    ########## Non-hierarchical parameters ########
    for (j in 1:1)
    {
        param[j] ~dnorm(0,1.0E-6)
    }
    tauBLOC ~ dgamma(0.001,0.001)
    tauPROC ~ dgamma(0.001,0.001)
 } # End of the jags model
 "

cat(NULL.MODEL , file = file.path(jags.dir, "GROWTH.MODEL.NULL"), sep=" ",
    fill = FALSE, labels = NULL, append = FALSE)

### Full models
MODEL.T<-
 "model {
    ############ Likelihood ###################
    for (i in 1:N.indiv) {
        ## linear regression
        Growth[i] ~ dnorm(theo.growth[i],tauPROC)
        theo.growth[i] <-param[1]
                       + param[2]*Tair[i]
                       + BLOC[num.bloc[i]]
    }
    ########### Hierarchical parameters ########
    for (n in 1:Nbloc)
    {
        BLOC[n]~dnorm(0,tauBLOC)
    }
    ########### Non-hierarchical parameters #######
    for (j in 1:2)
    {
        param[j] ~dnorm(0,1.0E-6)
    }
    tauBLOC ~ dgamma(0.001,0.001)
    tauPROC ~ dgamma(0.001,0.001)
} # End of the jags model
 "

cat(MODEL.T, file = file.path(jags.dir,"GROWTH.MODEL.T"), sep=" ",
    fill = FALSE, labels = NULL, append = FALSE)

### Full models
MODEL.INTER<-
 "model {
    ############ Likelihood ###################
    for (i in 1:N.indiv) {
        ## linear regression
        Growth[i] ~ dnorm(theo.growth[i],tauPROC)
        theo.growth[i] <-param[1]
                       + param[2]*Tair[i]
                       + param[3]*Lum[i]
                       + param[4]*Herb[i]
                       + param[5]*Lum[i]*Herb[i]
                       + param[6]*Tair[i]*Herb[i]
                       + param[7]*Tair[i]*Lum[i]
                       + param[8]*Tair[i]*Lum[i]*Herb[i]
                       + BLOC[num.bloc[i]]
    }
    ############################################
    ########### Hierarchical parameters ########
    for (n in 1:Nbloc)
    {
        BLOC[n]~dnorm(0,tauBLOC)
    }
    ###############################################
    ########### Non-hierarchical parameters #######
    for (j in 1:8)
    {
        param[j] ~dnorm(0,1.0E-6)
    }
    tauBLOC ~ dgamma(0.001,0.001)
    tauPROC ~ dgamma(0.001,0.001)
} # End of the jags model
 "

cat(MODEL.INTER, file = file.path(jags.dir,"GROWTH.MODEL.INTER"), sep=" ",
    fill = FALSE, labels = NULL, append = FALSE)


# generate inits function
inits.growth.null <- function(i, Nbloc){
  tauBLOC <- 200 +runif(1, min = -199, max = 199)
  inits<- list(param = runif(1,-1, 1),
               tauBLOC = tauBLOC,
               tauPROC = 1 + runif(1, min = -0.5, max = 0.5),
               BLOC = rnorm(Nbloc,mean=0,sd=sqrt(1/tauBLOC)))
return(inits)
  }

inits.growth.T <- function(i, Nbloc){
  tauBLOC <- 200 +runif(1, min = -199, max = 199)
  inits<- list(param = runif(2,-1, 1),
               tauBLOC = tauBLOC,
               tauPROC = 1 + runif(1, min = -0.5, max = 0.5),
               BLOC = rnorm(Nbloc,mean=0,sd=sqrt(1/tauBLOC)))
return(inits)
  }

inits.growth.inter <- function(i, Nbloc){
  tauBLOC <- 200 +runif(1, min = -199, max = 199)
  inits<- list(param = runif(8,-1, 1),
               tauBLOC = tauBLOC,
               tauPROC = 1 + runif(1, min = -0.5, max = 0.5),
               BLOC = rnorm(Nbloc,mean=0,sd=sqrt(1/tauBLOC)))
return(inits)
  }


inits.survival.null <- function(i, Nbloc){
  tauBLOC <- 1 +runif(1, min = -0.5, max = 0.5)
  inits<- list(param = runif(1,-1, 1),
               tauBLOC = tauBLOC,
               BLOC = rnorm(Nbloc,mean=0,sd=sqrt(1/tauBLOC)))
return(inits)
  }

inits.survival.inter <- function(i, Nbloc){
  tauBLOC <- 1 +runif(1, min = -0.5, max = 0.5)
  inits<- list(param = runif(8,-1, 1),
               tauBLOC = tauBLOC,
               BLOC = rnorm(Nbloc,mean=0,sd=sqrt(1/tauBLOC)))
return(inits)
  }

inits.survival.T <- function(i, Nbloc){
  tauBLOC <- 1 +runif(1, min = -0.5, max = 0.5)
  inits<- list(param = runif(2,-1, 1),
               tauBLOC = tauBLOC,
               BLOC = rnorm(Nbloc,mean=0,sd=sqrt(1/tauBLOC)))
return(inits)
  }



inter.param.mean.ci.vs.T <- function(p.name.1, p.name.2,
                                     res, temp.plot.CR, temp.plot,
                                     probs.CI = c(0.025,0.975)){
    matrix.res <-  matrix(0,nrow=length(res[,1]),ncol=length(temp.plot.CR))

    for (i in 1:length(res[,1]))
    {
        matrix.res[i, ] <-   res[i, p.name.1] + res[i, p.name.2]*temp.plot.CR
    }

    mean.vec <- apply(matrix.res,
                      MARGIN=2,
                      FUN=mean)
    mean.quant <- apply(matrix.res,
                        MARGIN=2,
                        FUN=quantile,
                        probs= probs.CI)
return(data.frame(temp= temp.plot,
                  mean = mean.vec,
                  quant.l = t(mean.quant)[, '2.5%'],
                  quant.h = t(mean.quant)[, '97.5%']))
}




format.df <-  function(i, sp.run2, species.lat, name.obj){
    load(file.path(output.dir, paste0(name.obj,sp.run2[i])))
    df <- data.frame(do.call('rbind', obj.plot),
                     param = c(rep('Ground vegetation direct', nrow(obj.plot[[1]])),
                               rep('Canopy direct', nrow(obj.plot[[1]])),
                               rep('Canopy indirect', nrow(obj.plot[[1]]))),
                     species = rep(species.lat[i], 3*nrow(obj.plot[[1]])))
    df[df$param != 'Ground vegetation direct', c('mean', 'quant.l', 'quant.h')] <-
        -df[df$param != 'Ground vegetation direct', c('mean', 'quant.l', 'quant.h')]
return(df)
}

theme_simple <-  function(){
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    strip.background = element_blank(),
    legend.key = element_blank(),
    strip.text.y = element_text(face = 'italic')
  ) +
  #draws x and y axis line
  theme(axis.line = element_line(color = 'black'))
}


#################################
##### JAGS for micro environment

jags.model.clim.L<-
 "
###################### Light model with  jags ######################
 model {
############ Likelihood ###################
  for (i in 1:N.indiv) {
    varmicroclim[i] ~ dnorm(mu[i],tauC)
    mu[i] <-  param[1] + param[2]*Lum[i] #
  }

## priors
  for (j in 1:2){
    param[j] ~dnorm(0,1.0E-6)#I(-5,5)
  }
  tauC~dgamma(0.001,0.001) # Non-informative prior
} # End of the jags model
 "
cat(jags.model.clim.L, file = file.path(jags.dir,"jags.model.clim.L"), sep=" ",
    fill = FALSE, labels = NULL, append = FALSE)


jags.model.clim.NULL<-
 "
###################### Light model with  jags ######################
 model {
############ Likelihood ###################
  for (i in 1:N.indiv) {
    varmicroclim[i] ~ dnorm(mu[i],tauC)
    mu[i] <-  param[1]
  }

## priors
  for (j in 1:1){
    param[j] ~dnorm(0,1.0E-6)#I(-5,5)
  }
  tauC~dgamma(0.001,0.001) # Non-informative prior
} # End of the jags model
 "
cat(jags.model.clim.NULL, file = file.path(jags.dir,"jags.model.clim.NULL"), sep=" ",
    fill = FALSE, labels = NULL, append = FALSE)

inits.clim <- function(i, n.param){
  tauC <- 1 +runif(1, min = -0.5, max = 0.5)
  param <- rep(0, n.param) +runif(n.param, min = -0.3, max = 0.3)
  inits<- list(param = param,
               tauC = tauC)
  return(inits)
  }



#######################
#######################
## FORMAT DATA CLIM


fun.data.season <-  function(df, type = 'quantile'){
require(dplyr)
data.clim <-  mutate(df, ID = paste(block, year))

if(type == 'quantile'){
data.season <-  group_by(data.clim, ID) %>%
        summarise(site= site[1],
                  block = block[1],
                  plot = plot[1],
                  year = year[1],
                  ground_vegetation_treatment = ground_vegetation_treatment[1],
                  transmitted_light= transmitted_light[1],
                  dds = dds[1],
                  VPD.season = quantile(VPD.season, probs = 0.8, na.rm = TRUE),
                  SWC.season = mean(SWC.season, na.rm = TRUE),
                  min.t.air = quantile(min.t.air, probs = 0.2, na.rm = TRUE))
}
return(data.season)
}

fun.data.canopy <-  function(df){
data.season <- fun.data.season(df)

data.season.NV <- filter(data.season, ground_vegetation_treatment == 'NV')
data.season.V <- filter(data.season, ground_vegetation_treatment == 'V')
data.season.NV <- mutate(data.season.NV, ID2 = paste(site, plot, year))
data.season.V <- mutate(data.season.V, ID2 = paste(site, plot, year))
## remove negative SWC
data.season.V$SWC.season[data.season.V$SWC.season<0] <-  NA
data.season.NV$SWC.season[data.season.NV$SWC.season<0] <-  NA
# data for canop effect
data.season.NV.effect <- group_by(data.season.NV, paste(site, year)) %>%
                  mutate(VPD.season.m = mean(VPD.season, na.rm =TRUE),
                         SWC.season.m = mean(SWC.season, na.rm = TRUE),
                         min.t.air.m = mean(min.t.air, na.rm = TRUE),
                         VPD_effect = VPD.season - VPD.season.m,
                         SWC_effect = SWC.season - SWC.season.m,
                         Tmin_effect = min.t.air - min.t.air.m
                         )

# compute delata V NV
data.canopy <-  as.data.frame(data.season.NV.effect)
return(data.canopy)
}

fun.data.herb <-  function(df){
data.season <- fun.data.season(df)
data.season.NV <- filter(data.season, ground_vegetation_treatment == 'NV')
data.season.V <- filter(data.season, ground_vegetation_treatment == 'V')
data.season.NV <- mutate(data.season.NV, ID2 = paste(site, plot, year))
data.season.V <- mutate(data.season.V, ID2 = paste(site, plot, year))
## remove negative SWC
data.season.V$SWC.season[data.season.V$SWC.season<0] <-  NA
data.season.NV$SWC.season[data.season.NV$SWC.season<0] <-  NA

data.season.V$delta.VPD <- data.season.V$VPD.season - data.season.NV$VPD.season
data.season.V$delta.SWC <- data.season.V$SWC.season - data.season.NV$SWC.season
data.season.V$delta.Tmin <- data.season.V$min.t.air - data.season.NV$min.t.air

data.herb <- as.data.frame(data.season.V)
return(data.herb)
}


format.jags.data.clim <- function(df, var.c){

data.herb.t <- df[!is.na(df[, var.c]) & !is.na(df[, 'transmitted_light']), ]
varmicroclim <-  data.herb.t[, var.c]
N.indiv<- length(varmicroclim) #
Lum <- 100- data.herb.t$transmitted_light
DDS <- data.herb.t$dds


jags.data<-list('varmicroclim' = varmicroclim,
                'N.indiv' = N.indiv,
                'Lum' = as.numeric(scale(Lum)),
                'DDS' = as.numeric(scale(DDS)))
return(jags.data)
}



format.jags.data.herb <- function(df, var.c){

data.herb.t <- df[!is.na(df[, var.c]) & !is.na(df[, 'transmitted_light']), ]
varmicroclim <-  data.herb.t[, var.c]
N.indiv<- length(varmicroclim) #
Lum <- data.herb.t$transmitted_light
DDS <- data.herb.t$dds


jags.data<-list('varmicroclim' = varmicroclim,
                'N.indiv' = N.indiv,
                'Lum' = log(Lum),
                'DDS' = DDS)
return(jags.data)
}
