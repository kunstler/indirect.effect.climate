####################################
####################################
## Analyse of climate variables
## and PCA
## Defossez et al. Oikos in revision

pca.clim <- function(clim_mean_name){
clim_mean <- read.csv(clim_mean_name)
## PCA of the climatic variables
rownames(clim_mean) <- clim_mean[,1]
pca.res <- princomp(clim_mean[,-1])
## Table for appendix S2
mat_pca <- matrix(NA, nrow = 5, ncol = 4)
mat_pca[2:5, ] <- unclass(loadings(pca.res))
mat_pca[1, ] <- pca.res$sdev^2/sum(pca.res$sdev^2)
rownames(mat_pca) <- c('Perc variance',  rownames(unclass(loadings(pca.res))))
colnames(mat_pca) <- colnames(unclass(loadings(pca.res)))
write.csv(mat_pca, 'output/mat_pca.csv', row.names = FALSE)
mat_pca
}

fig.clim.cor <- function(clim_mean_name, mat_pca_name){
clim_mean <- read.csv(clim_mean_name)
mat_pca <- as.matrix(read.csv(mat_pca_name))
# Figure S1 in Appendix 2: plot pairs correlations and PCA variance
par(mfrow=c(2,2), mgp = c(2, 0.5, 0), mar = c(3.5,3,1,1))
plot(clim_mean$DDS,clim_mean$SWC,xlab="Degree-day sum (>5.54°C)",ylab="SWC (%)")
text(2000,25,paste("rho=",round(cor.test(clim_mean$DDS,clim_mean$SWC,method="spearman")$estimate,2)," *"))
plot(clim_mean$DDS,clim_mean$VPD,xlab="Degree-day sum (>5.54°C)",ylab="VPD (Pa)")
text(1000,800,paste("rho=",round(cor.test(clim_mean$DDS,clim_mean$VPD,method="spearman")$estimate,2)," *"))
plot(clim_mean$DDS,clim_mean$Tmin,xlab="Degree-day sum (>5.54°C)",ylab="Tmin (°C)")
text(1000,8,paste("rho=",round(cor.test(clim_mean$DDS,clim_mean$Tmin,method="spearman")$estimate,2)," *"))
 par(mar = c(3.5,4,1,1))                 
barplot(mat_pca[1, ], ylab = "Percentage of variance \n per components", las = 3 )
}
