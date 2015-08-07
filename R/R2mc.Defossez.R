
##' .. compute R2m and R2c based on Nakagawa and Schielzeth (2013) MEE ..
##'
##' .. function work  with the R2jags object resulting from the most complex model of Defossez (Survival_T_L_H_LH_TH_TL_LTH_chene_ete_guiz for instance) and the three explicative variables Tair Lum Herb used in fitting the model..
##' @title 
##' @param model.res jags output from R2jags 
##' @param Tair explicative var used in model fit
##' @param Lum explicative var used in model fit
##' @param Herb explicative var used in model fit
##' @return a list $R2m is marginal R2 for each mcmc sample $R2c is conditional R2 for each mcmc sample
##' @author Kunstler
fun.R2mc.jags.full.model <- function(model.res,Tair,Lum,Herb){
    
param.name <- paste("param[",1:8,"]",sep='')
if(!all(length(Tair)==length(Herb),
    length(Tair)==length(Lum),
    length(Lum)==length(Herb)))
    stop('error Tair Herb and Lum have not the same length')
if(sum(! param.name %in%  colnames(model.res$BUGSoutput$sims.matrix))>0)
    stop('model.res do not have the good param param[1] to param[8]')
# create design matrix
Xmat <- cbind(rep(1, length(Tair)), Tair,Lum,Herb,Lum*Herb,Tair*Herb,Tair*Lum,Tair*Lum*Herb)
## create param matrix
beta.mat <- model.res$BUGSoutput$sims.matrix[,param.name]
# compute fixed effect var
fixed.expect <- array(dim = c(nrow(beta.mat), length(Tair)))
varF <- rep(NA, nrow(beta.mat))
for (i in 1:nrow(beta.mat)) {
    fixed.expect[i, ] <- beta.mat[i, ] %*% t(Xmat)
    varF[i] <- var(fixed.expect[i, ])
}
# var bloc effect and var for binomila model with logit
varBLOC <- 1/model.res$BUGSoutput$sims.matrix[,'tauBLOC']
varDist <- (pi^2)/3
# Calculate marginal R-squared
    Rm <- varF/(varF+varBLOC+varDist)
    # Calculate conditional R-squared (fixed effects+random effects/total variance)
    Rc <- (varF+varBLOC)/(varF+varBLOC+varDist)

return(list(R2m=Rm,R2c=Rc))
}


