model {
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
 