#1D Penalised Regression ---------------------------------------------------
# Libraries
library(tidyverse) # ggplot(.)
install.packages("glmnet")
library(glmnet) # for regularised regression

# Read dataset
horns <- read.csv("HornsRev.csv")

# Set Impact as factor
horns$Impact <- as.factor(df$Impact)

# Same as glmFitOD3 but covariates are scaled
glmFitOD3Scale <- glm(Nhat ~ Impact*scale(XPos) + Impact*scale(YPos) + scale(Depth), 
                      offset=log(Area), family=quasipoisson, data=df) 

# Create design/model matrix
# Note: ignore the intercept term as glmnet already introduces one by default.
xmatrix <- model.matrix(glmFitOD3Scale)
xmatrix <- xmatrix[, 2:ncol(xmatrix)]
head(xmatrix)

# RIDGE REGRESSION
# Set set
set.seed(123)

# Fit ridge regression model and use cross-validation to select tuning parameter lambda
ridge <- glmnet(xmatrix, df$Nhat, family="poisson", 
                offset=log(df$Area), alpha=0)
cvridge <- cv.glmnet(xmatrix, df$Nhat, family="poisson",
                     offset=log(df$Area), alpha=0, nfolds=10)
# Visualise results
par(mfrow=c(1, 2))
plot(ridge, xvar="lambda") 
abline(v=log(cvridge$lambda.min)) 
plot(cvridge)
abline(v=log(cvridge$lambda.min))
abline(v=log(cvridge$lambda.1se), lty=2) 

# trialled lambdas and chosen minimum lambda value
log(cvridge$lambda)
log(cvridge$lambda.min)

# quasi-Poisson fit
coefGLM <- as.data.frame(coef(glmFitOD3Scale))
colnames(coefGLM) <- "GLM"
coefGLM$Covariate <- row.names(coefGLM)

# Add 95% confidence intervals
confInt <- as.data.frame(confint(glmFitOD3Scale, level=0.95))
colnames(confInt) <- c("CI_Lower", "CI_Upper")
confInt$Covariate <- row.names(confInt)

# Merge together
coefGLM <- dplyr::inner_join(coefGLM, confInt, by="Covariate") 

# Ridge regression
coefRidge <- as.data.frame(as.matrix(coef(cvridge, s="lambda.min")))
colnames(coefRidge) <- "Ridge"
coefRidge$Covariate <- row.names(coefRidge)

# Merge data frames
mdlCoefs <- dplyr::inner_join(coefGLM, coefRidge, by="Covariate")
mdlCoefs=mdlCoefs

# Display differences
print(mdlCoefs)

# Test whether RR coefficients are within 95% CI for the glmFitOD3Scale model
ifelse(mdlCoefs$Ridge > mdlCoefs$CI_Lower & 
         mdlCoefs$Ridge < mdlCoefs$CI_Upper, TRUE, FALSE)

# LASSO

# Set set
set.seed(123)

# Fit LASSO model and use cross-validation to select tuning parameter lambda
lasso <- glmnet(xmatrix, df$Nhat, family="poisson", 
                offset=log(df$Area), alpha=1)
cvlasso <- cv.glmnet(xmatrix, df$Nhat, family="poisson",
                     offset=log(df$Area), alpha=1, nfolds=10)
# Visualise results
par(mfrow=c(1, 2))
plot(lasso, xvar="lambda") 
abline(v=log(cvlasso$lambda.min)) 
plot(cvlasso)
abline(v=log(cvlasso$lambda.min))
abline(v=log(cvlasso$lambda.1se), lty=2) 

# trialled lambdas and chosen minimum lambda value
log(cvlasso$lambda)
log(cvlasso$lambda.min)

# quasi-Poisson fit
coefGLM <- as.data.frame(coef(glmFitOD3Scale))
colnames(coefGLM) <- "GLM"
coefGLM$Covariate <- row.names(coefGLM)

# Add 95% confidence intervals
confInt <- as.data.frame(confint(glmFitOD3Scale, level=0.95))
colnames(confInt) <- c("CI_Lower", "CI_Upper")
confInt$Covariate <- row.names(confInt)

# Merge together
coefGLM <- dplyr::inner_join(coefGLM, confInt, by="Covariate") 

# LASSO
coeflasso <- as.data.frame(as.matrix(coef(cvlasso, s="lambda.min")))
colnames(coeflasso) <- "LASSO"
coeflasso$Covariate <- row.names(coeflasso)

# Merge data frames
mdlCoefs2 <- dplyr::inner_join(mdlCoefs, coeflasso, by="Covariate")
mdlCoefs2=mdlCoefs2

# Display differences
print(mdlCoefs)

# Test whether RR coefficients are within 95% CI for the glmFitOD3Scale model
ifelse(mdlCoefs$LASSO > mdlCoefs$CI_Lower & 
         mdlCoefs$LASSO < mdlCoefs$CI_Upper, TRUE, FALSE)

# ELASTIC NET 

# Set set
set.seed(123)

# Fit Elastic Net model and use cross-validation to select tuning parameter lambda
enet <- glmnet(xmatrix, df$Nhat, family="poisson", 
               offset=log(df$Area), alpha=0.6)
cvenet <- cv.glmnet(xmatrix, df$Nhat, family="poisson",
                    offset=log(df$Area), alpha=0.6, nfolds=10)
# Visualise results
par(mfrow=c(1, 2))
plot(enet, xvar="lambda") 
abline(v=log(cvenet$lambda.min)) 
plot(cvenet)
abline(v=log(cvenet$lambda.min))
abline(v=log(cvenet$lambda.1se), lty=2) 

# trialled lambdas and chosen minimum lambda value
log(cvenet$lambda)
log(cvenet$lambda.min)

# quasi-Poisson fit
coefGLM <- as.data.frame(coef(glmFitOD3Scale))
colnames(coefGLM) <- "GLM"
coefGLM$Covariate <- row.names(coefGLM)

# Add 95% confidence intervals
confInt <- as.data.frame(confint(glmFitOD3Scale, level=0.95))
colnames(confInt) <- c("CI_Lower", "CI_Upper")
confInt$Covariate <- row.names(confInt)

# Merge together
coefGLM <- dplyr::inner_join(coefGLM, confInt, by="Covariate") 

# Elastic Net
coefenet <- as.data.frame(as.matrix(coef(cvenet, s="lambda.min")))
colnames(coefenet) <- "ENet"
coefenet$Covariate <- row.names(coefenet)

# Merge data frames
mdlCoefs3 <- dplyr::inner_join(mdlCoefs2, coefenet, by="Covariate")
mdlCoefs3=mdlCoefs3

# Display differences
print(mdlCoefs3)

# Test whether Elastic Net coefficients are within 95% CI for the glmFitOD3Scale model
ifelse(mdlCoefs3$ENet > mdlCoefs3$CI_Lower & 
         mdlCoefs3$ENet < mdlCoefs3$CI_Upper, TRUE, FALSE)


# Compare three models' regression coefficient
# blue = GLM (with CIs) ------------------------------------------------------
# red = Ridge
# green = LASSO
# purple = Elastic net
ggplot(mdlCoefs3) + 
  geom_point(aes(x=Covariate, y=GLM), col="#377eb8") +
  geom_linerange(aes(x=Covariate, ymin=CI_Lower, ymax=CI_Upper), 
                 col="#377eb8") +
  geom_point(aes(x=Covariate, y=Ridge), col="#e41a1c") +
  geom_point(aes(x=Covariate, y=LASSO), col="#4daf4a") +
  geom_point(aes(x=Covariate, y=ENet), col="#984ea3") +
  ylab("Comparing regression coefficients") + 
  theme(axis.text.x=element_text(angle=90),
        legend.position="none")


# Write two functions that compute performance measures --------------------
RSS <- function(yObs, yFit)
{
  # RSS: Compute the residual sums of squares
  # yObs: the observed response y
  # yFit: the fitted response y
  return(sum((yObs - yFit)^2))
}
ASE <- function(yTruth, yFit)
{
  # ASE: Compute the mean average squared error 
  # yObs: the observed response y
  # yFit: the fitted response y
  return(mean((yTruth - yFit)^2))
}

# Libraries
library(mgcv) # penalised regression splines
install.packages('MRSea')
install.packages('devtools')
# install.packages("devtools")
devtools::install_github("lindesaysh/MRSea", ref="stable")
library(MRSea) # SALSA
install.packages('bias')
library(bias)

# Read simulated dataset from Moodle
df <- read.csv("SimulatedData.csv")
head(df)

#POLYNOMIAL
# Set up blank vectors to store RSS/ASE results
NSets <- max(df$ID)
RSS_Poly <- rep(NA, NSets)
ASE_Poly <- rep(NA, NSets)

# Set up plot with observed data points (in grey) and true function (heavy black line)
plot(df$x, df$response, pch=16, col="grey", cex=0.2,
     xlab="x", ylab="response")
lines(df$x, df$mu, lwd=2)

# Loop across all simulated datasets, fit and plot
for (i in seq(NSets))
{
  # Subset dataset
  dfSub <- subset(df, ID %in% i)
  
  # Fit model
  polyFit <- lm(response ~ x + poly(x, 6), data=dfSub)
  
  # Plot fitted response: lines() will add a line to an existing plot
  lines(dfSub$x, fitted(polyFit), col="#e41a1c", lwd=0.2)
  
  # Store RSS/ASE results
  RSS_Poly[i] <- RSS(dfSub$response, fitted(polyFit))
  ASE_Poly[i] <- ASE(dfSub$mu, fitted(polyFit))
}

#PRS
# Set up plot with observed data points (in grey) and true function (heavy black line)
plot(df$x, df$response, pch=16, col="grey", cex=0.2,
     xlab="x", ylab="response")
lines(df$x, df$mu, lwd=2)

# Loop across all simulated datasets, fit and plot
for (i in seq(NSets))
{
  # Subset dataset
  dfSub <- subset(df, ID %in% i)
  
  # Fit model
  prsFit <- lm(response ~ x + poly(x, 6), data=dfSub)
  
  # Plot fitted response: lines() will add a line to an existing plot
  lines(dfSub$x, fitted(polyFit), col="#e41a1c", lwd=0.2)
  
  # Store RSS/ASE results
  RSS_Poly[i] <- RSS(dfSub$response, fitted(polyFit))
  ASE_Poly[i] <- ASE(dfSub$mu, fitted(polyFit))
}

#SALSA
# Set up blank vectors to store results
NSets <- max(df$ID)
RSS_salsa <- rep(NA, NSets)
ASE_salsa <- rep(NA, NSets)
salsaknot <- rep(NA, NSets)
bias_salsa <- rep(NA, NSets)
var_salsa <- rep(NA, NSets)
# Loop across all sets of data
for (i in seq(100))
{
  # Subset dataset
  dfSub <- subset(df, ID %in% i)
  
  # Fit initial NULL model
  initialModel <- glm(response ~ 1, data=dfSub)
  
  # Set SALSA arguments
  varList <- c("x")
  salsa1DList <- list(fitnessMeasure="BIC", 
                      minKnots_1d=2, maxKnots_1d=40, 
                      startKnots_1d=10, degree=2, 
                      maxIterations=10, gaps=0)
  
  # Run SALSA
  salsa <- MRSea::runSALSA1D(initialModel=initialModel, 
                             salsa1dlist=salsa1DList, 
                             varlist=varList, 
                             factorlist=NULL,
                             datain=dfSub,
                             splineParams=NULL,
                             suppress.printout=TRUE)
  
  # Plot fitted response: lines() will add a line to an existing plot
  lines(dfSub$x, fitted(salsa$bestModel), col="blue", lwd=0.2)
  
  salsamodel <- salsa$bestModel
  
  # Store RSS/ASE results
  RSS_salsa[i] <- RSS(dfSub$response, fitted(salsamodel))
  ASE_salsa[i] <- ASE(dfSub$mu, fitted(salsamodel))
  salsaknot[i] <- length(salsamodel$splineParams[[2]]$knots)
  var_salsa[i] <- var(fitted(salsamodel))
}

# Penalised and Unpenalised regression splines ----------------------------

# Libraries
library(tidyverse) # ggplot(.)
library(mgcv) # gam(.)
library(splines) # bs(.) (B-splines)
library(MuMIn) # dredge(.)
library(fields) # quilt.plot(.)
library(lawstat) # runs.test(.)

# Read dataset
df <- read.csv("HornsRev.csv")

# Set Impact as factor
df$Impact <- as.factor(df$Impact)

# 1. Fit a quasi-Poisson penalised spline based GAM
PRS <- mgcv::gam(Nhat ~ s(XPos) + s(YPos) + s(Depth) + Impact, 
                 data=df, family=quasipoisson, offset=log(Area))
summary(PRS)

# 2. View the partial plots
par(mfrow=c(2,2))
# a
plot(PRS, shade=T, residuals=T, ylim=c(-10,10))
# b
plot(PRS, shade=T, residuals=T)

# 3. Visually compare the Depth relationship obtained using the PRS model to a 
# model that uses a B-spline basis function instead with just one knot at 20m
par(mfrow=c(1,2))

# Replace s(Depth) with bs(Depth, knots=20) in the model
PRS_B <- stats::update(PRS, .~. -s(Depth) + splines::bs(Depth, knots=20))
summary(PRS_B)
stats::termplot(PRS_B, se=T)

# 4. Carry out model selection on the terms in the working GAM using the dredge 
# function. Use QAIC to rank the different models.
library(MuMIn)
options(na.action="na.fail") # fail-safe
head(dredge(PRS, rank="QAIC", chat=summary(PRS)$scale))
summary(PRS)

# 5. Use the “best” model identified to predict over this fine grid and produce 
# a spatial plot of all the results (pre/post-impact both on the link and 
# response scale)
# Read prediction grid
predData <- read.csv("HornsRevPredictionData.csv")
# Predict on the link scale
NhatPredLink <- predict(PRS, newdata=predData, 
                        se=T, type="link")

# Predict on the response scale
NhatPredRes <- predict(PRS, newdata=predData, 
                       se=T, type="response")

# Spatial plots on link scale
library(fields) # quilt.plot(.)
par(mfrow=c(2,2))

for (impact in c("0","1")){
  bWant <- predData$Impact %in% impact # a boolean specifying phase
  fields::quilt.plot(predData$XPos[bWant], 
                     predData$YPos[bWant], 
                     NhatPredLink$fit[bWant], 
                     nrow=25, ncol=60,
                     zlim=range(fitted(PRS)),
                     main=paste0("Impact", impact))
}

# Spatial plots on response scale
library(fields) # quilt.plot(.)

for (impact in c(0,1)){
  bWant <- predData$Impact %in% impact # a boolean specifying phase
  fields::quilt.plot(predData$XPos[bWant], 
                     predData$YPos[bWant], 
                     NhatPredRes$fit[bWant], 
                     nrow=25, ncol=60,
                     zlim=range(fitted(PRS)),
                     main=paste0("Impact", impact))
}

# 6. Fit a GAM with an interaction between Impact and XPos.
PRS_Int <- mgcv::gam(Nhat ~ s(XPos, by=Impact) + s(YPos) + 
                       s(Depth) + Impact, 
                     data=df, family=quasipoisson, offset=log(Area))
summary(PRS_Int)

# 7. Use the model with the interaction term to predict over the 
# fine grid and produce a spatial plot of all the results 
# (pre/post-impact both on the link and response scale).
# Predict on the link scale
NhatPredLink2 <- predict(PRS_Int, newdata=predData, 
                         se=T, type="link")

# Predict on the response scale
NhatPredRes2 <- predict(PRS_Int, newdata=predData, 
                        se=T, type="response")

# Spatial plots on link scale
library(fields) # quilt.plot(.)
par(mfrow=c(2,2))

for (impact in c("0","1")){
  bWant <- predData$Impact %in% impact # a boolean specifying phase
  fields::quilt.plot(predData$XPos[bWant], 
                     predData$YPos[bWant], 
                     NhatPredLink2$fit[bWant], 
                     nrow=25, ncol=60,
                     zlim=range(fitted(PRS_Int)),
                     main=paste0("Impact", impact))
}

# Spatial plots on response scale
library(fields) # quilt.plot(.)

for (impact in c(0,1)){
  bWant <- predData$Impact %in% impact # a boolean specifying phase
  fields::quilt.plot(predData$XPos[bWant], 
                     predData$YPos[bWant], 
                     NhatPredRes2$fit[bWant], 
                     nrow=25, ncol=60,
                     zlim=range(fitted(PRS_Int)),
                     main=paste0("Impact", impact))
}

# 8. 
library(lawstat)
runs.test(residuals(PRS_Int, type="pearson"))


#2D Penalised Regression Splines -------------------------------------------
# Libraries
library(mgcv) # gam(.)
library(fields) # quilt.plot(.)

# Read dataset
df <- read.csv("HornsRev.csv")

# Set Impact as factor
df$Impact <- as.factor(df$Impact)

# 1. Produce pre/post-impact spatial plots
par(mfrow=c(1,2))
for (impact in c("0","1")){
  bWant <- df$Impact %in% impact # a boolean specifying impact
  fields::quilt.plot(df$XPos[bWant], 
                     df$YPos[bWant], 
                     df$Nhat[bWant], 
                     nrow=25, ncol=60,
                     zlim=range(df$Nhat)/10,
                     main=paste0("Impact", impact))
}

# 2. Fit a quasi-Poisson penalised spline based GAM with a two 
# dimensional smoother for the spatial coordinates and a one 
# dimensional smoother for Depth, whilst including Impact
PRS_2D <- mgcv::gam(Nhat ~ s(XPos, YPos) + s(Depth) + Impact, 
                    data=df, family=quasipoisson, offset=log(Area))
summary(PRS_2D)

# 3. View partial plots for PRS_2D
plot(PRS_2D, shade=T, residuals=F)

# 4. Make pre/post-impact predictions on this fine grid using PRS_2D. 
dfpred <- read.csv("HornsRevPredictionData.csv")

# Link scale predictions
PRSPredl<- predict(PRS_2D, newdata=dfpred, type="link")

# Response scale predictions
PRSPredr<- predict(PRS_2D, newdata=dfpred, type="response")

# Produce spatial plot of the results on link scale
zlimits <- c(min(PRSPredl), max(PRSPredr))
par(mfrow=c(2,2))
for (impact in c(0, 1)){
  bWant <- dfpred$Impact %in% impact # a boolean specifying impact
  fields::quilt.plot(dfpred$XPos[bWant], 
                     dfpred$YPos[bWant], 
                     PRSPredl[bWant], 
                     nrow=25, ncol=60,
                     zlim=zlimits,
                     main=paste0("Impact", impact))}

# Produce spatial plot of the results on response scale
for (impact in c(0, 1)){
  bWant <- dfpred$Impact %in% impact # a boolean specifying impact
  fields::quilt.plot(dfpred$XPos[bWant], 
                     dfpred$YPos[bWant], 
                     PRSPredr[bWant], 
                     nrow=25, ncol=60,
                     zlim=zlimits,
                     main=paste0("Impact", impact))}

# 5. Fit a similar model to PRS_2D, but now add Impact as an interaction term in
# the two dimensional spatial smoother 
PRS_2DInt <- mgcv::gam(Nhat ~ s(XPos, YPos, by = as.factor(Impact)) + 
                         s(Depth) + Impact, data=df, family=quasipoisson, 
                       offset=log(Area))
summary(PRS_2DInt)

# 6. View partial plots for PRS_2DInt
plot(PRS_2DInt, shade=T, residuals=F)

# 7. Make pre/post-impact predictions on this fine grid using PRS_2DInt. 
# Link scale predictions
PRSPredlint<- predict(PRS_2DInt, newdata=dfpred, type="link")

# Response scale predictions
PRSPredrint<- predict(PRS_2DInt, newdata=dfpred, type="response")

# Produce spatial plot of the results on link scale
zlimits <- c(min(PRSPredlint), max(PRSPredrint))
par(mfrow=c(2,2))
for (impact in c(0, 1)){
  bWant <- dfpred$Impact %in% impact # a boolean specifying impact
  fields::quilt.plot(dfpred$XPos[bWant], 
                     dfpred$YPos[bWant], 
                     PRSPredlint[bWant], 
                     nrow=25, ncol=60,
                     zlim=zlimits,
                     main=paste0("Impact", impact))}

# Produce spatial plot of the results on response scale
for (impact in c(0, 1)){
  bWant <- dfpred$Impact %in% impact # a boolean specifying phase
  fields::quilt.plot(dfpred$XPos[bWant], 
                     dfpred$YPos[bWant], 
                     PRSPredrint[bWant], 
                     nrow=25, ncol=60,
                     zlim=zlimits,
                     main=paste0("Impact", impact))}

# 8. Plot the locations of the observation points and the fine prediction grid 
ggplot(dfpred)+geom_point(aes(x=XPos, y=YPos))

# 9. Try increasing (set k=15) and decreasing (set k=40) the limit on the 
# degrees of freedom for the two dimensional spatial smoother
PRS_2Dlim15 <- mgcv::gam(Nhat ~ s(XPos, YPos, by = as.factor(Impact), k=15) + 
                           s(Depth) + Impact, data=df, family=quasipoisson, 
                         offset=log(Area))
summary(PRS_2Dlim15)

PRS_2Dlim40 <- mgcv::gam(Nhat ~ s(XPos, YPos, by = as.factor(Impact), k=40) + 
                           s(Depth) + Impact, data=df, family=quasipoisson, 
                         offset=log(Area))
summary(PRS_2Dlim40)

# 10. Make pre/post-impact predictions on this fine grid using the models above
# PRS_2Dlim15:

# Link scale predictions
PRSPredllim15<- predict(PRS_2Dlim15, newdata=dfpred, type="link")

# Response scale predictions
PRSPredrlim15<- predict(PRS_2Dlim15, newdata=dfpred, type="response")

# Produce spatial plot of the results on link scale
zlimits <- c(min(PRSPredllim15), max(PRSPredrlim15))
par(mfrow=c(2,2))
for (impact in c(0, 1)){
  bWant <- dfpred$Impact %in% impact # a boolean specifying impact
  fields::quilt.plot(dfpred$XPos[bWant], 
                     dfpred$YPos[bWant], 
                     PRSPredllim15[bWant], 
                     nrow=25, ncol=60,
                     zlim=zlimits,
                     main=paste0("Impact", impact))}

# Produce spatial plot of the results on response scale
for (impact in c(0, 1)){
  bWant <- dfpred$Impact %in% impact # a boolean specifying phase
  fields::quilt.plot(dfpred$XPos[bWant], 
                     dfpred$YPos[bWant], 
                     PRSPredrlim15[bWant], 
                     nrow=25, ncol=60,
                     zlim=zlimits,
                     main=paste0("Impact", impact))}

# PRS_2Dlim40:

# Link scale predictions
PRSPredllim40<- predict(PRS_2Dlim40, newdata=dfpred, type="link")

# Response scale predictions
PRSPredrlim40 <- predict(PRS_2Dlim40, newdata=dfpred, type="response")

# Produce spatial plot of the results on link scale
zlimits <- c(min(PRSPredllim40), max(PRSPredrlim40))
par(mfrow=c(2,2))
for (impact in c(0, 1)){
  bWant <- dfpred$Impact %in% impact # a boolean specifying impact
  fields::quilt.plot(dfpred$XPos[bWant], 
                     dfpred$YPos[bWant], 
                     PRSPredllim40[bWant], 
                     nrow=25, ncol=60,
                     zlim=zlimits,
                     main=paste0("Impact", impact))}

# Produce spatial plot of the results on response scale
for (impact in c(0, 1)){
  bWant <- dfpred$Impact %in% impact # a boolean specifying phase
  fields::quilt.plot(dfpred$XPos[bWant], 
                     dfpred$YPos[bWant], 
                     PRSPredrlim40[bWant], 
                     nrow=25, ncol=60,
                     zlim=zlimits,
                     main=paste0("Impact", impact))}
