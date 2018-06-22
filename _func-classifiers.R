enet = function (X, Y, alpha, lambda = NULL, family, X.test = NULL, 
  Y.test = NULL, filter = "p.value", topranked = 50, keepVar = NULL, pop.prev=0.5,
   cutoff=NULL){
  library(glmnet)
  library(limma)
  library(pROC)
  library(OptimalCutpoints)
  if (filter == "none") {
    X1 <- X
  }
  if (filter == "p.value") {
    design <- model.matrix(~Y)
    fit <- eBayes(lmFit(t(X), design))
    top <- topTable(fit, coef = 2, adjust.method = "BH", n = nrow(fit))
    if(sum(colnames(top) %in% "ID")){
      X1 <- X[, top$ID[1:topranked]]
    } else {
      X1 <- X[, rownames(top)[1:topranked]]
    }
  }
  if (is.null(keepVar)) {
    penalty.factor <- rep(1, ncol(X1))
    X2 <- X1
  } else {
    X1 <- X1[, setdiff(colnames(X1), keepVar)]
    X2 <- as.matrix(cbind(X1, X[, keepVar]))
    colnames(X2) <- c(colnames(X1), keepVar)
    penalty.factor <- c(rep(1, ncol(X1)), rep(0, length(keepVar)))
  }
  if (family == "binomial") {
    fit <- glmnet(X2, Y, family = "binomial", alpha = alpha, 
      penalty.factor = penalty.factor)
    cv.fit <- cv.glmnet(X2, Y, family = "binomial")
    if (is.null(lambda)) {
      lambda = cv.fit$lambda.min
    } else {
      lambda = lambda
    }
    Coefficients <- coef(fit, s = lambda)
    Active.Index <- which(Coefficients[, 1] != 0)
    Active.Coefficients <- Coefficients[Active.Index, ]
    enet.panel <- names(Active.Coefficients)[-1]
    enet.panel.length <- length(enet.panel)
  }
  if (family == "multinomial") {
    fit <- glmnet(X2, Y, family = "multinomial", alpha = alpha, 
      type.multinomial = "grouped", penalty.factor = penalty.factor)
    cv.fit <- cv.glmnet(X2, Y, family = "multinomial")
    if (is.null(lambda)) {
      lambda = cv.fit$lambda.min
    }
    else {
      lambda = lambda
    }
    Coefficients <- coef(fit, s = lambda)
    Active.Index <- which(Coefficients[[1]][, 1] != 0)
    Active.Coefficients <- Coefficients[[1]][Active.Index, 
      ]
    enet.panel <- names(Active.Coefficients)[-1]
    enet.panel.length <- length(enet.panel)
  }
  if (!is.null(X.test)) {
    library(pROC)
    library(OptimalCutpoints)
    probs <- predict(fit, newx = X.test[, colnames(X2)], 
      s = lambda, type = "response")
    predictResponse <- unlist(predict(fit, newx = X.test[, colnames(X2)], 
      s = lambda, type = "class"))
    if (family == "binomial") {
      if(is.null(cutoff)){
        perfTest <- tperformance(weights = as.numeric(as.matrix(probs)), 
          trueLabels = Y.test, pop.prev = pop.prev)
      } else {
        pred <- rep(NA, length(Y.test))
        pred[as.numeric(probs) >= 0.5] <- levels(Y.test)[2]
        pred[as.numeric(probs) < 0.5] <- levels(Y.test)[1]
        pred <- factor(pred, levels(Y.test))
        perfTest <- calcPerf(pred=pred, truth=Y.test, pop.prev)
      }
    }
    else {
      mat <- table(factor(as.character(predictResponse), 
        levels = levels(Y.test)), Y.test)
      mat2 <- mat
      diag(mat2) <- 0
      classError <- colSums(mat2)/colSums(mat)
      er <- sum(mat2)/sum(mat)
      ber <- mean(classError)
      perfTest <- c(classError, er, ber)
      names(perfTest) <- c(names(classError), "ER", "BER")
    }
  } else {
    perfTest <- predictResponse <- probs <- NA
  }
  return(list(X = X, Y = Y, fit = fit, enet.panel = enet.panel, 
    lambda = lambda, alpha = alpha, family = family, probs = probs, 
    Active.Coefficients = Active.Coefficients, perfTest = perfTest, 
    predictResponse = predictResponse, filter = filter, topranked = topranked, 
    keepVar = keepVar))
}


tperformance = function (weights, trueLabels, pop.prev){
  df = data.frame(prob = as.numeric(weights), status = model.matrix(~factor(as.character(trueLabels), 
    levels = levels(trueLabels)))[, 2])
  roc.score = roc(response = df$status, predictor = weights, 
    plot = FALSE, percent = TRUE, na.rm = TRUE, direction = "<")
  optimal.cutpoint.Youden <- optimal.cutpoints(X = "prob", 
    status = "status", tag.healthy = 0, methods = "Youden", 
    data = df, control = control.cutpoints(), ci.fit = FALSE, 
    conf.level = 0.95, trace = FALSE, pop.prev = pop.prev)
  optimalValues <- round(c(summary(optimal.cutpoint.Youden)$p.table$Global$Youden[[1]][1:5, 
    ], roc.score$auc/100), 3)
  names(optimalValues) <- c(names(optimalValues)[-length(names(optimalValues))], 
    "AUC")
  optimalValues
}




calcPerf = function(pred, truth, prev){
  mat <- table(pred, truth)
  if(rownames(mat)[1] != colnames(mat)[1])
    stop("check levels of inputs")
  
  spec = mat[1,1]/sum(mat[, 1])
  sens = mat[2,2]/sum(mat[, 2])
  npv = ((1-prev)*spec)/(prev*(1-sens)+(1-prev)*spec)
  ppv = (prev*sens)/((prev*sens)+(1-prev)*spec)
  accuracy <- sum(diag(mat))/sum(mat)
  perf = c(sens, spec, npv, ppv, accuracy)
  names(perf) <- c("sens", "spec", "npv", "ppv", "accuracy")
  perf
}












enetCV = function (X, Y, alpha, lambda, M, folds, progressBar, family, 
  filter, topranked, keepVar) 
{
  library(pROC)
  library(limma)
  library(glmnet)
  library(OptimalCutpoints)
  probs <- predictResponseList <- enet.panel <- list()
  if (progressBar == TRUE) 
    pb <- txtProgressBar(style = 3)
  for (i in 1:M) {
    if (progressBar == TRUE) 
      setTxtProgressBar(pb, i/M)
    omit = folds[[i]]
    X.train = X[-omit, , drop = FALSE]
    Y.train = Y[-omit]
    if (filter == "none") {
      X.train1 <- X.train
    }
    if (filter == "p.value") {
      design <- model.matrix(~Y.train)
      fit <- eBayes(lmFit(t(X.train), design))
      top <- topTable(fit, coef = 2, adjust.method = "BH", 
        n = nrow(fit))
      X.train1 <- X.train[, rownames(top)[1:topranked], 
        drop = FALSE]
    }
    if (is.null(keepVar)) {
      penalty.factor <- rep(1, ncol(X.train1))
      X.train2 <- X.train1
    }
    else {
      X.train1 <- X.train1[, setdiff(colnames(X.train1), 
        keepVar)]
      X.train2 <- as.matrix(cbind(X.train1, X.train[, keepVar]))
      colnames(X.train2) <- c(colnames(X.train1), keepVar)
      penalty.factor <- c(rep(1, ncol(X.train1)), rep(0, 
        length(keepVar)))
    }
    X.test1 = X[omit, colnames(X.train2), drop = FALSE]
    if (family == "binomial") {
      fit <- glmnet(X.train2, Y.train, family = "binomial", 
        alpha = alpha, penalty.factor = penalty.factor)
      cv.fit <- cv.glmnet(X.train2, Y.train, family = "binomial")
      if (is.null(lambda)) {
        lambda = cv.fit$lambda.min
      }
      else {
        lambda = lambda
      }
      Coefficients <- coef(fit, s = lambda)
      Active.Index <- which(Coefficients[, 1] != 0)
      Active.Coefficients <- Coefficients[Active.Index, 
        ]
      enet.panel[[i]] <- names(Active.Coefficients)[-1]
    }
    if (family == "multinomial") {
      fit <- glmnet(X.train2, Y.train, family = "multinomial", 
        alpha = alpha, type.multinomial = "grouped", 
        penalty.factor = penalty.factor)
      cv.fit <- cv.glmnet(X.train2, Y.train, family = "multinomial")
      if (is.null(lambda)) {
        lambda = cv.fit$lambda.min
      }
      else {
        lambda = lambda
      }
      Coefficients <- coef(fit, s = lambda)
      Active.Index <- which(Coefficients[[1]][, 1] != 0)
      Active.Coefficients <- Coefficients[[1]][Active.Index, 
        ]
      enet.panel[[i]] <- names(Active.Coefficients)[-1]
    }
    probs[[i]] <- predict(fit, newx = X.test1[, colnames(X.train2)], s = lambda, 
      type = "response")
    predictResponseList[[i]] <- predict(fit, newx = X.test1[, colnames(X.train2)], 
      s = lambda, type = "class")
  }
  predictResponse <- unlist(predictResponseList)
  if (family == "binomial") {
    probs <- unlist(probs)
    trueLabels = Y[unlist(folds)]
    library(pROC)
    library(OptimalCutpoints)
    perf <- amritr::tperformance(weights = probs, trueLabels = trueLabels)
  }
  else {
    trueLabels = Y[unlist(folds)]
    mat <- table(factor(trueLabels, levels(Y)), factor(predictResponse, 
      levels(Y)))
    mat2 <- mat
    diag(mat2) <- 0
    classError <- colSums(mat2)/colSums(mat)
    er <- sum(mat2)/sum(mat)
    ber <- mean(classError)
    perf <- c(classError, er, ber)
    names(perf) <- c(names(classError), "ER", "BER")
  }
  return(list(probs = probs, trueLabels = trueLabels, perf = perf, 
    enet.panel = enet.panel, predictResponse = predictResponse))
}

















perf.enet = function (object, validation = c("Mfold", "loo"), M = 5, iter = 10, 
  threads = 4, progressBar = TRUE) 
{
  library(dplyr)
  library(tidyr)
  X = object$X
  Y = object$Y
  n = nrow(X)
  alpha = object$alpha
  family = object$family
  lambda = object$lambda
  filter = object$filter
  topranked = object$topranked
  keepVar = object$keepVar
  if (validation == "Mfold") {
    folds <- lapply(1:iter, function(i) createFolds(Y, k = M))
    require(parallel)
    cl <- parallel::makeCluster(mc <- getOption("cl.cores", 
      threads))
    parallel::clusterExport(cl, varlist = c("enetCV", "enet", 
      "X", "Y", "alpha", "lambda", "M", "folds", "progressBar", 
      "family", "filter", "topranked", "keepVar"), envir = environment())
    cv <- parallel::parLapply(cl, folds, function(foldsi, 
      X, Y, alpha, lambda, M, progressBar, family, filter, 
      topranked, keepVar) {
      enetCV(X = X, Y = Y, alpha = alpha, lambda = lambda, 
        M = M, folds = foldsi, progressBar = progressBar, 
        family = family, filter = filter, topranked = topranked, 
        keepVar = keepVar)
    }, X, Y, alpha, lambda, M, progressBar, family, filter, 
      topranked, keepVar) %>% amritr::zip_nPure()
    parallel::stopCluster(cl)
    perf <- do.call(rbind, cv$perf) %>% as.data.frame %>% 
      gather(ErrName, Err) %>% dplyr::group_by(ErrName) %>% 
      dplyr::summarise(Mean = mean(Err), SD = sd(Err))
  }
  else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
    cv <- enetCV(X, Y, alpha, lambda, M, folds, progressBar, 
      family, filter, topranked, keepVar)
    perf <- data.frame(Mean = cv$perf) %>% mutate(ErrName = rownames(.))
    perf$SD <- NA
  }
  result = list()
  result$folds = folds
  result$probs = cv$probs
  result$trueLabels = cv$trueLabels
  result$panels = cv$enet.panel
  result$perf = perf
  method = "enet.mthd"
  result$meth = "enet.mthd"
  class(result) = c("perf", method)
  return(invisible(result))
}



performance = function (weights, trueLabels, pop.prev){
  df = data.frame(prob = as.numeric(weights), status = model.matrix(~factor(as.character(trueLabels), 
    levels = levels(trueLabels)))[, 2])
  roc.score = roc(response = df$status, predictor = weights, 
    plot = FALSE, percent = TRUE, na.rm = TRUE, direction = "<")
  optimal.cutpoint.Youden <- optimal.cutpoints(X = "prob", 
    status = "status", tag.healthy = 0, methods = "Youden", 
    data = df, control = control.cutpoints(), ci.fit = FALSE, 
    conf.level = 0.95, trace = FALSE, pop.prev = pop.prev)
  optimalValues <- round(c(summary(optimal.cutpoint.Youden)$p.table$Global$Youden[[1]][1:5, 
    ], roc.score$auc/100), 3)
  names(optimalValues) <- c(names(optimalValues)[-length(names(optimalValues))], 
    "AUC")
  optimalValues
}



rforest = function (X, Y, X.test = NULL, Y.test = NULL, family = "binomial", 
  filter = "p.value", topranked = 50){
  library(randomForest)
  library(limma)
  library(pROC)
  library(OptimalCutpoints)
  X1 <- X
  featIndex <- colnames(X1)
  names(featIndex) <- paste("f", 1:ncol(X1), sep = "_")
  colnames(X1) <- names(featIndex)
  if (filter == "none") {
    X1 <- X1
  }
  if (filter == "p.value") {
    design <- model.matrix(~Y)
    fit <- eBayes(lmFit(t(X1), design))
    top <- topTable(fit, coef = 2, adjust.method = "BH", 
      n = nrow(fit))
    X1 <- X1[, rownames(top)[1:topranked]]
  }
  rf.panel <- as.character(featIndex[colnames(X1)])
  y.train0 <- Y[rownames(X1)]
  fit = randomForest(y.train0 ~ ., data = X1, importance = TRUE, 
    proximity = TRUE)
  if (!is.null(X.test)) {
    colnames(X.test) <- paste("f", 1:ncol(X1), sep = "_")
    predictResponse <- predict(fit, as.matrix(X.test), type = "response")
    probs <- predict(fit, as.matrix(X.test), type = "vote")[, 
      levels(Y)[2]]
    if (family == "binomial") {
      perfTest <- amritr::tperformance(weights = as.numeric(as.matrix(probs)), 
        trueLabels = Y.test)
    }
    if (family == "multinomial") {
      mat <- table(Y.test, predictResponse)
      mat2 <- mat
      diag(mat2) <- 0
      classError <- colSums(mat2)/colSums(mat)
      er <- sum(mat2)/sum(mat)
      ber <- mean(classError)
      perfTest <- c(classError, er, ber)
      names(perfTest) <- c(names(classError), "ER", "BER")
    }
  }
  else {
    perfTest <- predictResponse <- NA
  }
  return(list(X = X, Y = Y, fit = fit, perfTest = perfTest, 
    rf.panel = rf.panel, probs=probs, predictResponse = predictResponse, 
    family = family, filter = filter, topranked = topranked))
}














