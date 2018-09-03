##################################################################
#
# crossValidationPerf.R
# Author: Amrit Singh
# Date: August 11, 2015
# Biomarker pipleine for microarry and rnaseq data 
#
##################################################################
## load libraries
require(snowfall, quietly=TRUE)  # required for parallelization
require(limma, quietly=TRUE)
require(caret, quietly=TRUE)   ## require to make cv folds
require(glmnet, quietly=TRUE)
require(randomForest, quietly = TRUE)
require(ROCR, quietly = TRUE)
require(mixOmics, quietly = TRUE)

#----------------------------------------------------
# Biomarker pipeline run in parallel (Main function)
#----------------------------------------------------
crossValidationPerf = function(X = X, Y = Y, folds = folds, nfolds = nfolds, topRanked = topRanked, covMat = NULL, alpha = alpha, cpus = cpus, file.path="~/Documents/RNA-SeqAnalysis/FinalAnalysis/src", file.name="crossValidationPerf.R", preSelectFolds = NULL){
  
  if(is.null(preSelectFolds)){
    folds.iter <- list()
    ## create cross-validation folds
    for(j in 1:nfolds){
      fsamples <- createFolds(Y, k=folds)
      folds.iter[[j]] <- fsamples
    } 
  } else {
    folds.iter <- preSelectFolds
  }
  
  sfInit(parallel = TRUE, cpus = cpus)
  sfLibrary(limma); sfLibrary(glmnet); sfLibrary(randomForest); 
  sfLibrary(ROCR); sfLibrary(mixOmics);
  setwd(file.path)
  sfSource(file.name)
  sfExport("X", "Y", "topRanked", "covMat", "alpha")
  result <- sfClusterApplyLB(folds.iter, classifier)
  sfStop()
  
  Perf <- list(Enet=lapply(result, function(i) i[[1]]$Enet), Rf=lapply(result, function(i) i[[1]]$Rf))
  Error <- list(Enet=lapply(result, function(i) i[[2]]$Enet), Rf=lapply(result, function(i) i[[2]]$Rf))
  PanelFeatures <- list(Enet=lapply(result, function(i) i[[3]]$Enet), Rf=lapply(result, function(i) i[[3]]$Rf))
  PanelLength <- list(Enet=lapply(result, function(i) i[[4]]$Enet), Rf=lapply(result, function(i) i[[4]]$Rf))
  SubjScore <- list(Enet=lapply(result, function(i) i[[5]]$Enet), Rf=lapply(result, function(i) i[[5]]$Rf))
  SubjPredictedClass <- list(Enet=lapply(result, function(i) i[[6]]$Enet), Rf=lapply(result, function(i) i[[6]]$Rf))
  SubjTrueClass <- list(Enet=lapply(result, function(i) i[[7]]$Enet), Rf=lapply(result, function(i) i[[7]]$Rf))

  return(list(Perf=Perf, Error=Error, PanelFeatures=PanelFeatures, PanelLength=PanelLength,
              SubjScore=SubjScore, SubjPredictedClass=SubjPredictedClass, SubjTrueClass=SubjTrueClass, folds.iter=folds.iter))
}

#x <-  classifier(folds.iter[[2]])
#iteration <- folds.iter[[2]]

#--------------------------------------------
# Deep cross-validation; limma + elastic net & randomForest
#--------------------------------------------
classifier = function(iteration){  
  result <- lapply(1 : length(iteration), function(z){
    test <- X[as.numeric(iteration[[z]]),, drop=FALSE]
    fsamples2 <- iteration
    fsamples2[[z]] <- NULL
    train <- X[as.numeric(unlist(fsamples2)),, drop=FALSE]
    
    ## Univariate filtering
    shortList <- filterFeatures(train=train, Y=Y, topRanked=topRanked, covMat=covMat)
    ## Classification
    classify <- classification(keepFeatures = shortList$sigGenes, train=train, Y=Y, alpha=alpha,
                               test=test, thres=thres, covMat=covMat)
    enetResultsList <- classify$fit.enet
    rfResultsList <- classify$fit.rf
    return(list(enetResultsList=enetResultsList, rfResultsList=rfResultsList))
  })
  enetResults <- lapply(result, function(w) w[[1]])
  rfResults <- lapply(result, function(w) w[[2]])
  
  Perf <- list(Enet=lapply(enetResults, function(y) y$Perf), Rf=lapply(rfResults, function(y) y$Perf))
  Error <- list(Enet=lapply(enetResults, function(y) y$Error), Rf=lapply(rfResults, function(y) y$Error))
  PanelFeatures <- list(Enet=lapply(enetResults, function(y) y$PanelFeatures), Rf=lapply(rfResults, function(y) y$PanelFeatures))
  PanelLength <- list(Enet=lapply(enetResults, function(y) y$PanelLength), Rf=lapply(rfResults, function(y) y$PanelLength))
  SubjScore <- list(Enet=lapply(enetResults, function(y) y$SubjScore), Rf=lapply(rfResults, function(y) y$SubjScore))
  SubjPredictedClass <- list(Enet=lapply(enetResults, function(y) y$SubjPredictedClass), Rf=lapply(rfResults, function(y) y$SubjPredictedClass))
  SubjTrueClass <- list(Enet=lapply(enetResults, function(y) y$SubjTrueClass), Rf=lapply(rfResults, function(y) y$SubjTrueClass))
  return(list(Perf=Perf, Error=Error, PanelFeatures=PanelFeatures, PanelLength=PanelLength,
              SubjScore=SubjScore, SubjPredictedClass=SubjPredictedClass, SubjTrueClass=SubjTrueClass))
}

filterFeatures = function(train=train, Y=Y, topRanked=topRanked, covMat=covMat){ 
  Y.fold <- Y[rownames(train)]
  if(is.matrix(covMat)){
    design.train <- cbind(model.matrix(~Y.fold), covMat[names(Y.fold), ])
  } else {
    design.train <- model.matrix(~Y.fold)
  }
  
  eset <- t(train)[, names(Y.fold), drop=FALSE]
  fit <- eBayes(lmFit(eset, design.train, method="ls"))
  top00 <- topTable(fit, coef=2, adjust.method='BH', n=nrow(fit))
  top0 <- top00[top00$P.Value < 0.05, ]
  top <- top0[order(abs(top0$logFC), decreasing = TRUE), ,drop=FALSE]
  top$ID <- rownames(top)
  if(nrow(top) < topRanked){
    sigGenes <- top$ID
  } else {
    sigGenes <- top$ID[1:topRanked] 
  }
return(list(sigGenes=sigGenes))
}

classification = function(keepFeatures=keepFeatures, train=train, Y=Y, alpha=alpha,
                          test=test, thres=thres, covMat=covMat){
  ## Training and Test dataset
  if(is.matrix(covMat)){
    X.train0 <- train[, keepFeatures, drop=FALSE]
    X.test0 <- test[, keepFeatures, drop=FALSE]
    X.train <- cbind(covMat[rownames(X.train0), ], X.train0)
    X.test <- cbind(covMat[rownames(X.test0), ], X.test0)
    y.train <- Y[rownames(X.train)]
    y.test <- Y[rownames(X.test)]
  } else {
    X.train <- train[, keepFeatures, drop=FALSE]
    y.train <- Y[rownames(X.train)]
    X.test <- test[, keepFeatures, drop=FALSE]
    y.test <- Y[rownames(X.test)]
  }
  
  ## Run Elastic net classifier
  fit.enet <- enetClassifier(X.train = X.train, y.train = y.train, alpha = alpha, X.test = X.test, y.test = y.test)
  ## Run RF classifier
  rfpanelLength <- fit.enet$PanelLength  
  fit.rf <- rfClassifier(X.train = X.train, y.train = y.train, rfpanelLength = rfpanelLength, X.test = X.test, y.test = y.test)

  return(list(fit.enet=fit.enet, fit.rf = fit.rf))  
}



enetClassifier = function(X.train = X.train, y.train = y.train, alpha = alpha, X.test = X.test, y.test = y.test){
  ## Run Elastic net classifier
  fit <- glmnet(X.train, y.train, family="binomial", alpha = alpha)
  cv.fit <- cv.glmnet(X.train, y.train, family="binomial")
  s = cv.fit$lambda.min
  Coefficients <- coef(fit, s = s)
  Active.Index <- which(Coefficients != 0)
  Active.Coefficients  <- Coefficients[Active.Index,]
  enet.panel <- names(Active.Coefficients)[-1]
  enet.panel.length <- length(enet.panel)
  probs0 <- predict(fit, newx=as.matrix(X.test), s=s, type='response')
  enet.predictResponse <- as.vector(predict(fit, newx=as.matrix(X.test), s=s, type='class'))
  enet.probs <- as.vector(probs0)
  names(enet.probs) <- names(enet.predictResponse) <- rownames(probs0)    
  # Calculate classification performance measurse
  validation.pred = prediction(enet.probs, y.test)
  ## error rate
  confusion = table(pred = enet.predictResponse, truth = y.test)
  error <- c(performance(validation.pred, 'auc')@y.values[[1]], diag(confusion)/colSums(confusion), sum(diag(confusion))/sum(confusion))
  names(error) <- c("AUC", paste(colnames(confusion), "ErrorRate", sep="."), "Overall.ErrorRate")
  enetPerf <- list()
  enetPerf[[1]] <- data.frame(threshold = performance(validation.pred, "sens")@x.values[[1]],
                          sensitivity = performance(validation.pred, "sens")@y.values[[1]],
                          specificity = performance(validation.pred, "spec")@y.values[[1]],
                          ppv = performance(validation.pred, "ppv")@y.values[[1]],
                          npv = performance(validation.pred, "npv")@y.values[[1]])
  enetPerf[[2]] <- error
  enetPerf[[3]] <- enet.panel
  enetPerf[[4]] <- enet.panel.length
  enetPerf[[5]] <- enet.probs
  enetPerf[[6]] <- enet.predictResponse
  enetPerf[[7]] <- y.test
  names(enetPerf) <- c("Perf", "Error", "PanelFeatures", "PanelLength", "SubjScore", "SubjPredictedClass", "SubjTrueClass")
  return(enetPerf)
}

rfClassifier = function(X.train = X.train, y.train = y.train, rfpanelLength = rfpanelLength,
                        X.test = X.test, y.test = y.test){
  ## Random Forest
  X.train0 <- X.train
  X.test0 <- X.test
  colnames(X.train0) <- colnames(X.test0) <-paste("Feature", 1:ncol(X.train0), sep=".")
  rf = randomForest(y.train ~ ., data=X.train0, importance=TRUE, proximity=TRUE)
  
  ## Select important variables
  rf.panel = c(rownames(rf$importance)[order(rf$importance[,4],decreasing=T)][1:rfpanelLength])
  refine.rf = randomForest(y.train ~ ., data = X.train0[, rf.panel, drop=FALSE], importance=T, proximity=T)
  rf.predictResponse <- predict(refine.rf, as.matrix(X.test0), type='response')
  rf.probs <- predict(refine.rf, as.matrix(X.test0), type='vote')[, 2]
  
  # Calculate classification performance measurse
  validation.pred = prediction(rf.probs, y.test)
  ## error rate
  confusion = table(pred = rf.predictResponse, truth = y.test)
  error <- c(performance(validation.pred, 'auc')@y.values[[1]], diag(confusion)/colSums(confusion), sum(diag(confusion))/sum(confusion))
  names(error) <- c("AUC", paste(colnames(confusion), "ErrorRate", sep="."), "Overall.ErrorRate")
  rfPerf <- list()
  rfPerf[[1]] <- data.frame(threshold = performance(validation.pred, "sens")@x.values[[1]],
                              sensitivity = performance(validation.pred, "sens")@y.values[[1]],
                              specificity = performance(validation.pred, "spec")@y.values[[1]],
                              ppv = performance(validation.pred, "ppv")@y.values[[1]],
                              npv = performance(validation.pred, "npv")@y.values[[1]])
  rfPerf[[2]] <- error
  rfPerf[[3]] <- colnames(X.train)[!is.na(match(colnames(X.train0), rf.panel))]
  rfPerf[[4]] <- rfpanelLength
  rfPerf[[5]] <- rf.probs
  rfPerf[[6]] <- rf.predictResponse
  rfPerf[[7]] <- y.test
  names(rfPerf) <- c("Perf", "Error", "PanelFeatures", "PanelLength", "SubjScore", "SubjPredictedClass", "SubjTrueClass")
  return(rfPerf)
}