# Simulate data for BFSO_D-15-00532
library(vegan)
library(ggplot2)

# nifH 
nifh_nh <- as.matrix(rnorm(3, mean=2.32e7, sd=0.02e7 * sqrt(3)))
nifh_yl <- as.matrix(rnorm(3, mean=4.36e7, sd=0.10e7 * sqrt(3)))
nifh_zl <- as.matrix(rnorm(3, mean=1.37e8, sd=0.06e8 * sqrt(3)))
nifh_cl <- as.matrix(rnorm(3, mean=1.89e8, sd=0.05e8 * sqrt(3)))
# AOB
aob_nh <- as.matrix(rnorm(3, mean=1.90e6, sd=0.05e6 * sqrt(3)))
aob_yl <- as.matrix(rnorm(3, mean=4.58e6, sd=0.03e6 * sqrt(3)))
aob_zl <- as.matrix(rnorm(3, mean=8.24e6, sd=0.04e6 * sqrt(3)))
aob_cl <- as.matrix(rnorm(3, mean=6.13e6, sd=0.09e6 * sqrt(3)))
# AOA
aoa_nh <- as.matrix(rnorm(3, mean=1.40e3, sd=0.01e3 * sqrt(3)))
aoa_yl <- as.matrix(rnorm(3, mean=2.16e3, sd=0.02e3 * sqrt(3)))
aoa_zl <- as.matrix(rnorm(3, mean=7.15e3, sd=0.01e3 * sqrt(3)))
aoa_cl <- as.matrix(rnorm(3, mean=4.12e3, sd=0.03e3 * sqrt(3)))
# nirK
nirk_nh <- as.matrix(rnorm(3, mean=3.56e5, sd=0.10e5 * sqrt(3)))
nirk_yl <- as.matrix(rnorm(3, mean=7.45e5, sd=0.02e5 * sqrt(3)))
nirk_zl <- as.matrix(rnorm(3, mean=1.39e6, sd=0.04e6 * sqrt(3)))
nirk_cl <- as.matrix(rnorm(3, mean=2.48e6, sd=0.05e6 * sqrt(3)))
# narG
narg_nh <- as.matrix(rnorm(3, mean=1.59e6, sd=0.02e6 * sqrt(3)))
narg_yl <- as.matrix(rnorm(3, mean=5.37e6, sd=0.06e6 * sqrt(3)))
narg_zl <- as.matrix(rnorm(3, mean=7.17e6, sd=0.04e6 * sqrt(3)))
narg_cl <- as.matrix(rnorm(3, mean=1.60e7, sd=0.04e7 * sqrt(3)))

#Combine data to dataframe
nifh <- rbind(nifh_nh, nifh_yl, nifh_zl, nifh_cl)
aoa <- rbind(aoa_nh, aoa_yl, aoa_zl, aoa_cl)
aob <- rbind(aob_nh, aob_yl, aob_zl, aob_cl)
nirk <- rbind(nirk_nh, nirk_yl, nirk_zl, nirk_cl)
narg <- rbind(narg_nh, narg_yl, narg_zl, narg_cl)

ngenes <- data.frame(nifh, aoa, aob, nirk, narg)
ngenes["site"] <- rep(c("nh", "yl", "zl", "cl"), each=3)

#run PCA
ngenes_pca <- princomp(ngenes[-6], cor=TRUE)

#format scores for PCA plot
var_x <- ngenes_pca$sdev[1]^2/sum(ngenes_pca$sdev[1]^2, ngenes_pca$sdev[2]^2, ngenes_pca$sdev[3]^2, ngenes_pca$sdev[4]^2, ngenes_pca$sdev[5]^2)
var_y <- ngenes_pca$sdev[2]^2/sum(ngenes_pca$sdev[1]^2, ngenes_pca$sdev[2]^2, ngenes_pca$sdev[3]^2, ngenes_pca$sdev[4]^2, ngenes_pca$sdev[5]^2)
                                                                    
ngenes_x <- as.data.frame(ngenes_pca$scores[,1])
ngenes_y <- as.data.frame(ngenes_pca$scores[,2])
site <- ngenes["site"]
ngenes_pca_scores <- data.frame(ngenes_x, ngenes_y, site)
colnames(ngenes_pca_scores) <- c("x", "y", "site")

#plot PCA
ggplot(ngenes_pca_scores, aes(x = x, y = y, color = site)) + 
  geom_point(size = 3) + 
  xlab(paste("PCA1 (", round(var_x, 3), "%)")) +
  ylab(paste("PCA2 (", round(var_y, 3), "%)"))

#ANOVA
#Run 1-way ANOVA
nifh_anova = aov(nifh ~ site, data = ngenes)
aoa_anova = aov(aoa ~ site, data = ngenes)
aob_anova = aov(aob ~ site, data = ngenes)
nirk_anova = aov(nirk ~ site, data = ngenes)
narg_anova = aov(narg ~ site, data = ngenes)

#Post-hoc Test
summary(nifh_anova); TukeyHSD(nifh_anova)
summary(aoa_anova); TukeyHSD(aoa_anova)
summary(aob_anova); TukeyHSD(aob_anova)
summary(nirk_anova); TukeyHSD(nirk_anova)
summary(narg_anova); TukeyHSD(narg_anova)

par(mfrow = c(1, 5))
boxplot(nifh ~ site, data = ngenes, main = "nifH")
boxplot(aoa ~ site, data = ngenes, main = "AOA")
boxplot(aob ~ site, data = ngenes, main = "AOB")
boxplot(nirk ~ site, data = ngenes, main = "nirK")
boxplot(narg ~ site, data = ngenes, main = "narG")
