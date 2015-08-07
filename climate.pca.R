####################################
####################################
## Analyse of climate variables
## and PCA
## Defossez et al. Oikos in revision

#####
output.dir <- 'output'
dir.create(output.dir, showWarnings = FALSE)
figs.dir <- 'figs'
dir.create(figs.dir, showWarnings = FALSE)
data.dir <- 'data'


## read mean annual climate data
clim.mean <- read.csv(file.path(data.dir, "data.clim.mean.csv"),sep=",")



## PCA of the climatic variables
rownames(clim.mean) <- clim.mean[,1]
pca.res <- princomp(clim.mean[,-1])

## Table for appendix S2
mat.pca <- matrix(NA, nrow = 5, ncol = 4)
mat.pca[2:5, ] <- unclass(loadings(pca.res))
mat.pca[1, ] <- pca.res$sdev^2/sum(pca.res$sdev^2)
rownames(mat.pca) <- c('Perc variance',  rownames(unclass(loadings(pca.res))))
colnames(mat.pca) <- colnames(unclass(loadings(pca.res)))

write.csv(mat.pca, file = file.path(output.dir, "res.pca.climate.csv"))

# Figure S1 in Appendix 2: plot pairs correlations and PCA variance
pdf(file.path(figs.dir, "clim.cor.pdf"))
par(mfrow=c(2,2), mgp = c(2, 0.5, 0))
plot(clim.mean$DDS,clim.mean$SWC,xlab="Degree-day sum (>5.54°C)",ylab="SWC (%)")
text(2000,25,paste("rho=",round(cor.test(clim.mean$DDS,clim.mean$SWC,method="spearman")$estimate,2)," *"))
plot(clim.mean$DDS,clim.mean$VPD,xlab="Degree-day sum (>5.54°C)",ylab="VPD (Pa)")
text(1000,800,paste("rho=",round(cor.test(clim.mean$DDS,clim.mean$VPD,method="spearman")$estimate,2)," *"))
plot(clim.mean$DDS,clim.mean$Tmin,xlab="Degree-day sum (>5.54°C)",ylab="Tmin (°C)")
text(1000,8,paste("rho=",round(cor.test(clim.mean$DDS,clim.mean$Tmin,method="spearman")$estimate,2)," *"))
barplot(mat.pca[1, ], ylab = "Percentage of variance \n per components", las = 3 )
dev.off()
