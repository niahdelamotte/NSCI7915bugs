#loadpackages
library(car)
library(lme4)
library(lmerTest)
library(MuMIn)
library(emmeans)
library(ggplot2)
library(ggfortify)
library(ggbiplot)
library(plyr)
library(multcomp)
library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(tidyverse)
library(cluster)
library(dendextend)

#removing NA values from dataset and naming
bugs <- na.omit(read.csv('NSCI7915bugsmeeting.csv')[,1:15])

#looking at pca assumptions
plot(lm(bugs[,5:9]))
#generating multi-panel plot of the above assumption plots
par(mfrow = c(1, 1))
par(mfrow = c(2, 2), mar = c(4,4,4,4))
g1 <- plot(lm(bugs[,5:9]), which = c(1))
g2 <- plot(lm(bugs[,5:9]), which = c(2))
g3 <- plot(lm(bugs[,5:9]), which = c(3))
g4 <- plot(lm(bugs[,5:9]), which = c(5))

#checking for correlation between developmental stages
growmaturity <- bugs[,5:9]
head(growmaturity)
ngrowmat <- scale(growmaturity)
head(ngrowmat)
corr_matrix <- cor(ngrowmat)
#shows some low positive correlations, depending on which developmental stage
ggcorrplot(corr_matrix, type = "lower", lab = TRUE)

#growth rates columns selected and PCA findings where variance occurs most
princomp(bugs[,5:9],scores=TRUE)
summary(princomp(bugs[,5:9]))
#naming growthrate principal  components  for lmer test
growthrate <- princomp(bugs[,5:9])$scores[,1]
#lmer test, fixed effects are predictors temp and diet, sex to account for 
#effect of sex at adulthood following developmental tracking
growmodel <- lmer(growthrate ~ Temp + Diet + (1|Sex),data=bugs)
summary(growthrate ~ Temp + Diet + (1|Sex),data=bugs)
grow_intervals <- confint(growmodel, level = 0.95)
print(grow_intervals)
#strength of random effect:
rand((lmer(growthrate ~ Temp + Diet + (1|Sex),data=bugs)))

#additional pca to get better visualisation of the instar phases
buggrow <- princomp(bugs[,5:9])$scores[,1]
grow <- bugs[,5:9]
head(grow)
ngrow <- scale(grow)
head(ngrow)
corr_matrix <- cor(ngrow)
#generates a correlation matrix of the instar phases
ggcorrplot(corr_matrix)
grow.pca <- princomp(corr_matrix)
#can see the proportion of variance for each principal component
summary(grow.pca)
#how each growth stage  contribute to principal components
grow.pca$loadings[,1:3]

#princomp plots for development
#plotting using scores 
growthrate<- bugs[5:9]
pca_grow <- prcomp(growthrate, scale. = TRUE)
#plot of developmental stages under temperature and diet effects
autoplot(pca_grow, data=bugs, main='Principal Components for Developmental Stages:
         Temperature and Diet Effects', 
         colour='Temp', shape='Diet', size=2,loadings=TRUE, loadings.colour='black', 
         loadings.label = TRUE, loadings.label.size=4, loadings.label.colour = 'black',loadings.label.vjust = 
           1.8,scale=0)
#sex effect on the developmental stages' PCA
autoplot(pca_grow, data=bugs, colour='Sex',size=2,main='Principal Components
         Sex Effect')

#allometric relationships
hist(bugs$PL)
hist(bugs$BL)
hist(bugs$PW)
#will log body measurements to make distribution more normal

#pca assumptions for body measurements
plot(lm(bugs[,13:15]))
#generating a multi-panel plot of the assumption plots
par(mfrow = c(1, 1))
par(mfrow = c(2, 2), mar = c(4,4,4,4))
g1 <- plot(lm(bugs[,13:15]), which = c(1))
g2 <- plot(lm(bugs[,13:15]), which = c(2))
g3 <- plot(lm(bugs[,13:15]), which = c(3))
g4 <- plot(lm(bugs[,13:15]), which = c(5))

#pca of bug body measurements
summary(princomp(bugs[,13:15]))

#linear mixed effects model - allometric relationship of proboscis length and 
#body length accounting for temp, diet and sex
PLBLmodel <- lmer(log(PL) ~ Temp + Diet + log(BL) + (1|Sex),data=bugs)
summary((lmer(log(PL) ~ Temp + Diet + log(BL) + (1|Sex),data=bugs)))
#looking at random effect now for BL-PL relationship
rand((lmer(log(PL) ~ Temp + Diet + log(BL) + (1|Sex),data=bugs)))
#coef var for PL
sd(log(bugs$PL))
#95% CI
PLBL_intervals <- confint(PLBLmodel, level = 0.95)
print(PLBL_intervals)
#naming for biplot - visualise relationships and see influence on PCs
pca_allo <- prcomp(bugs[,13:15], scale. = TRUE)
autoplot(pca_allo, data=bugs, main='Principal Components for Body Measurements', 
         colour='Temp', shape='Diet', size=2, loadings=TRUE, loadings.colour='black', 
         loadings.label = TRUE, loadings.label.size=4,loadings.label.vjust = 
           1.8,scale=0)

#PW and BL allometry
PWBLmodel <-lmer(log(PW) ~ Temp + Diet + log(BL) + (1|Sex),data=bugs)
summary((lmer(log(PW) ~ Temp + Diet + log(BL) + (1|Sex),data=bugs)))
#random effect of sex
rand((lmer(log(PW) ~ Temp + Diet + log(BL) + (1|Sex),data=bugs)))
#95% CI
PWBL_intervals <- confint(PWBLmodel, level = 0.95)
print(PWBL_intervals)

#PW and PL
PWPLmodel <- lmer(log(PW) ~ Temp + Diet + log(PL) + (1|Sex),data=bugs)
summary((lmer(log(PW) ~ Temp + Diet + log(PL) + (1|Sex),data=bugs)))
#random effects
rand((lmer(log(PW) ~ Temp + Diet + log(PL) + (1|Sex),data=bugs)))
#95% CI
PWPL_intervals <- confint(PWPLmodel, level = 0.95)
print(PWPL_intervals)

#pca matrix and variable breakdown
bugbody <- princomp(bugs[,13:15])$scores[,1]
body <- bugs[,13:15]
head(body)
nbody <- scale(body)
head(nbody)
corr_matrix <- cor(nbody)
#correlation matrix generated to look at body measurements
ggcorrplot(corr_matrix,type = "lower", lab = TRUE)
body.pca <- princomp(corr_matrix)
#look at variance across components
summary(body.pca)
#influence of the body measurements on each component
#moderate negative value
body.pca$loadings[,1:3]

#PL
#hierarchial cluster analysis using hclust - want euclidean method for dist, 
#ward.D2 method for hclust
PLhcl <- dist(log(bugs$PL))
PLdf <- PLhcl
PLdf <- scale(PLdf)
head(PLdf)
PLd <- dist(PLdf, method = "euclidean")
PLhc1 <- hclust(PLd, method = "ward.D2" )
#cutting tree into three clusters
fit <- cutree(PLhc1, k = 3 )
par(mfrow = c(1, 1))
#checking difference between start and cut tree
plot(PLhc1) 
hPL <- plot(PLhc1, cex = 0.6, hang = -1, labels=FALSE,main="Dendrogram of Proboscis Length", xlab="Proboscis Length (Ward's Method)")
rect.hclust(PLhc1, k = 3, border = 2:5)
#scree plot to check number of clusters is appropriate - seems to be 3
fviz_nbclust(PLdf, FUN = hcut, method = "wss")
fviz_nbclust(PLdf, FUN = hcut, method = "silhouette")
PLgap_stat <- clusGap(PLdf, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(PLgap_stat)

#BL
#logging body length
BLhcl <- dist(log(bugs$BL))
BLdf <- BLhcl
#scaling before cluster analysis
BLdf <- scale(BLdf)
#checking data
head(BLdf)
#specifying euclidean distancing method
BLd <- dist(BLdf, method = "euclidean")
#ward method for hierarchical clustering 
BLhc1 <- hclust(BLd, method = "ward.D2" )
#cutting tree into three clusters
fit <- cutree(BLhc1, k = 3 )
#checking difference between start and cut tree
plot(BLhc1) 
hBL <- plot(BLhc1, cex = 0.6, hang = -1, labels=FALSE,main="Dendrogram of Body Length",xlab="Body Length (Ward's Method)")
#hierarchical plot with clusters in boxes
rect.hclust(BLhc1, k = 3, border = 2:5)
#scree plot to check number of clusters is appropriate - seems to be 3
fviz_nbclust(BLdf, FUN = hcut, method = "wss")
fviz_nbclust(BLdf, FUN = hcut, method = "silhouette")
BLgap_stat <- clusGap(BLdf, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(BLgap_stat)

#PW
PWhcl <- dist(log(bugs$PW))
PWdf <- PWhcl
PWdf <- scale(PWdf)
head(PWdf)
PWd <- dist(PWdf, method = "euclidean")
PWhc1 <- hclust(PWd, method = "ward.D2" )
#cutting tree into three clusters
fit <- cutree(PWhc1, k = 3)
#checking difference between start and cut tree
plot(PWhc1) 
hPW <- plot(PWhc1, cex = 0.6, hang = -1, labels=FALSE, main="Dendrogram of Pronotum Width", xlab="Pronotum Width (Ward's Method)")
#h cluster plot with boxes around clusters
rect.hclust(PWhc1, k = 6, border = 2:5)
#scree plot to check number of clusters is appropriate - seems to be 6
fviz_nbclust(PWdf, FUN = hcut, method = "wss")
fviz_nbclust(PWdf, FUN = hcut, method = "silhouette")
PWgap_stat <- clusGap(PWdf, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(PWgap_stat)

#dendrograms one figure
par(mfrow = c(1, 1))
par(mfrow = c(3,1), mar = c(mar = c(0.5, 4, 4, 0.5)))
hPW <- plot(PWhc1, cex = 0.6, hang = -1, labels=FALSE, main="Dendrogram of Pronotum Width", xlab="Pronotum Width (Ward's Method)")
rect.hclust(PWhc1, k = 6, border = 2:5)
hBL <- plot(BLhc1, cex = 0.6, hang = -1, labels=FALSE,main="Dendrogram of Body Length",xlab="Body Length (Ward's Method)")
rect.hclust(BLhc1, k = 3, border = 2:5)
hPL <- plot(PLhc1, cex = 0.6, hang = -1, labels=FALSE,main="Dendrogram of Proboscis Length", xlab="Proboscis Length (Ward's Method)")
rect.hclust(PLhc1, k = 3, border = 2:5)
