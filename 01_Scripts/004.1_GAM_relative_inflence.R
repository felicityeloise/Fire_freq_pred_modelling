# Edited by Felicity Charles 19/03/2025

# Code for gam.hp() edited lines 33:40 as the function was unable to handle the bam function used for model fitting. 



function (mod = fire_gam, iv = NULL, type = "dev", commonality = FALSE) 
{
  if (!inherits(mod, "gam")) 
    stop("gam.hp only supports gam objects at the moment")
  Formu <- strsplit(as.character(mod$call$formula)[3], "")[[1]]
  if ("*" %in% Formu) 
    stop("Please put the interaction term as a new variable (i.e. link variables by colon(:)) and  avoid the asterisk (*) and colon(:) in the original model")
  if (!is.null(attr(mod$terms, "offset"))) {
    ivname <- strsplit(as.character(mod$formula)[3], "[/+]")[[1]][-(attr(mod$terms, 
                                                                         "offset") + 1)]
  }
  if (is.null(attr(mod$terms, "offset"))) {
    ivname <- strsplit(gsub(" ", "", as.character(mod$formula)[3]), 
                       "[/+]")[[1]]
  }
  iv.name <- ivname
  if (type == "adjR2") 
    outr2 <- summary(mod)$r.sq
  if (type == "dev") 
    outr2 <- summary(mod)$dev.expl
  r2type <- row.names(outr2)
  nr2type <- length(r2type)
  if (nr2type == 0) {
    nr2type <- 1
    if (commonality) {
      r2type <- "commonality.analysis"
    }
    else {
      r2type <- "hierarchical.partitioning"
    }
  }
  dat <- na.omit(mod$model)
  
  mod_null <-  gam(Sentinel_ff ~ 1,data = Pres_back[trainSet,],family = "poisson")
  if (is.null(iv)) {
    nvar <- length(iv.name)
    if (nvar < 2) 
      stop("Analysis not conducted. Insufficient number of predictors.")
    totalN <- 2^nvar - 1
    binarymx <- matrix(0, nvar, totalN)
    for (i in 1:totalN) {
      binarymx <- creatbin(i, binarymx)
    }
    outputList <- list()
    outputList[[1]] <- outr2
    for (k in 1:nr2type) {
      commonM <- matrix(nrow = totalN, ncol = 3)
      for (i in 1:totalN) {
        tmp.name <- iv.name[as.logical(binarymx[, i])]
        to_add <- paste("~", paste(c(tmp.name, if (!is.null(mod$offset)) as.character(attr(mod$terms, 
                                                                                           "variables")[attr(mod$terms, "offset") + 1])), 
                                   collapse = " + "))
        modnew <- stats::update(object = mod_null, data = dat, 
                                formula = as.formula(to_add))
        if (type == "dev") 
          commonM[i, 2] <- summary(modnew)$dev.expl
        if (type == "adjR2") 
          commonM[i, 2] <- summary(modnew)$r.sq
      }
      commonlist <- vector("list", totalN)
      seqID <- vector()
      for (i in 1:nvar) {
        seqID[i] = 2^(i - 1)
      }
      for (i in 1:totalN) {
        bit <- binarymx[1, i]
        if (bit == 1) 
          ivname <- c(0, -seqID[1])
        else ivname <- seqID[1]
        for (j in 2:nvar) {
          bit <- binarymx[j, i]
          if (bit == 1) {
            alist <- ivname
            blist <- genList(ivname, -seqID[j])
            ivname <- c(alist, blist)
          }
          else ivname <- genList(ivname, seqID[j])
        }
        ivname <- ivname * -1
        commonlist[[i]] <- ivname
      }
      for (i in 1:totalN) {
        r2list <- unlist(commonlist[i])
        numlist <- length(r2list)
        ccsum <- 0
        for (j in 1:numlist) {
          indexs <- r2list[[j]]
          indexu <- abs(indexs)
          if (indexu != 0) {
            ccvalue <- commonM[indexu, 2]
            if (indexs < 0) 
              ccvalue <- ccvalue * -1
            ccsum <- ccsum + ccvalue
          }
        }
        commonM[i, 3] <- ccsum
      }
      orderList <- vector("list", totalN)
      index <- 0
      for (i in 1:nvar) {
        for (j in 1:totalN) {
          nbits <- sum(binarymx[, j])
          if (nbits == i) {
            index <- index + 1
            commonM[index, 1] <- j
          }
        }
      }
      outputcommonM <- matrix(nrow = totalN + 1, ncol = 2)
      totalRSquare <- sum(commonM[, 3])
      for (i in 1:totalN) {
        outputcommonM[i, 1] <- round(commonM[commonM[i, 
                                                     1], 3], digits = 4)
        outputcommonM[i, 2] <- round((commonM[commonM[i, 
                                                      1], 3]/totalRSquare) * 100, digits = 2)
      }
      outputcommonM[totalN + 1, 1] <- round(totalRSquare, 
                                            digits = 4)
      outputcommonM[totalN + 1, 2] <- round(100, digits = 4)
      rowNames <- NULL
      for (i in 1:totalN) {
        ii <- commonM[i, 1]
        nbits <- sum(binarymx[, ii])
        cbits <- 0
        if (nbits == 1) 
          rowName <- "Unique to "
        else rowName <- "Common to "
        for (j in 1:nvar) {
          if (binarymx[j, ii] == 1) {
            if (nbits == 1) 
              rowName <- paste(rowName, iv.name[j], 
                               sep = "")
            else {
              cbits <- cbits + 1
              if (cbits == nbits) {
                rowName <- paste(rowName, "and ", sep = "")
                rowName <- paste(rowName, iv.name[j], 
                                 sep = "")
              }
              else {
                rowName <- paste(rowName, iv.name[j], 
                                 sep = "")
                rowName <- paste(rowName, ", ", sep = "")
              }
            }
          }
        }
        rowNames <- c(rowNames, rowName)
      }
      rowNames <- c(rowNames, "Total")
      rowNames <- format.default(rowNames, justify = "left")
      colNames <- format.default(c("Fractions", " % Total"), 
                                 justify = "right")
      dimnames(outputcommonM) <- list(rowNames, colNames)
      VariableImportance <- matrix(nrow = nvar, ncol = 4)
      for (i in 1:nvar) {
        VariableImportance[i, 3] <- round(sum(binarymx[i, 
        ] * (commonM[, 3]/apply(binarymx, 2, sum))), 
        digits = 4)
      }
      VariableImportance[, 1] <- outputcommonM[1:nvar, 
                                               1]
      VariableImportance[, 2] <- VariableImportance[, 
                                                    3] - VariableImportance[, 1]
      total = round(sum(VariableImportance[, 3]), digits = 4)
      VariableImportance[, 4] <- round(100 * VariableImportance[, 
                                                                3]/total, 2)
      dimnames(VariableImportance) <- list(iv.name, c("Unique", 
                                                      "Average.share", "Individual", "I.perc(%)"))
      if (commonality) {
        outputList[[k + 1]] <- outputcommonM
      }
      else {
        outputList[[k + 1]] <- VariableImportance
      }
    }
  }
  else {
    nvar <- length(iv)
    ilist <- names(iv)
    if (is.null(ilist)) {
      names(iv) <- paste("X", 1:nvar, sep = "")
    }
    else {
      whichnoname <- which(ilist == "")
      names(iv)[whichnoname] <- paste("X", whichnoname, 
                                      sep = "")
    }
    ilist <- names(iv)
    ivlist <- ilist
    iv.name <- ilist
    ivID <- matrix(nrow = nvar, ncol = 1)
    for (i in 0:nvar - 1) {
      ivID[i + 1] <- 2^i
    }
    totalN <- 2^nvar - 1
    binarymx <- matrix(0, nvar, totalN)
    for (i in 1:totalN) {
      binarymx <- creatbin(i, binarymx)
    }
    outputList <- list()
    outputList[[1]] <- outr2
    for (k in 1:nr2type) {
      commonM <- matrix(nrow = totalN, ncol = 3)
      for (i in 1:totalN) {
        tmpname <- iv.name[as.logical(binarymx[, i])]
        tmp.name <- unlist(iv[names(iv) %in% tmpname])
        to_add <- paste("~", paste(c(tmp.name, if (!is.null(mod$offset)) as.character(attr(mod$terms, 
                                                                                           "variables")[attr(mod$terms, "offset") + 1])), 
                                   collapse = " + "))
        modnew <- stats::update(object = mod_null, data = dat, 
                                formula = as.formula(to_add))
        if (type == "dev") 
          commonM[i, 2] <- summary(modnew)$dev.expl
        if (type == "adjR2") 
          commonM[i, 2] <- summary(modnew)$r.sq
      }
      commonlist <- vector("list", totalN)
      seqID <- vector()
      for (i in 1:nvar) {
        seqID[i] = 2^(i - 1)
      }
      for (i in 1:totalN) {
        bit <- binarymx[1, i]
        if (bit == 1) 
          ivname <- c(0, -seqID[1])
        else ivname <- seqID[1]
        for (j in 2:nvar) {
          bit <- binarymx[j, i]
          if (bit == 1) {
            alist <- ivname
            blist <- genList(ivname, -seqID[j])
            ivname <- c(alist, blist)
          }
          else ivname <- genList(ivname, seqID[j])
        }
        ivname <- ivname * -1
        commonlist[[i]] <- ivname
      }
      for (i in 1:totalN) {
        r2list <- unlist(commonlist[i])
        numlist <- length(r2list)
        ccsum <- 0
        for (j in 1:numlist) {
          indexs <- r2list[[j]]
          indexu <- abs(indexs)
          if (indexu != 0) {
            ccvalue <- commonM[indexu, 2]
            if (indexs < 0) 
              ccvalue <- ccvalue * -1
            ccsum <- ccsum + ccvalue
          }
        }
        commonM[i, 3] <- ccsum
      }
      orderList <- vector("list", totalN)
      index <- 0
      for (i in 1:nvar) {
        for (j in 1:totalN) {
          nbits <- sum(binarymx[, j])
          if (nbits == i) {
            index <- index + 1
            commonM[index, 1] <- j
          }
        }
      }
      outputcommonM <- matrix(nrow = totalN + 1, ncol = 2)
      totalRSquare <- sum(commonM[, 3])
      for (i in 1:totalN) {
        outputcommonM[i, 1] <- round(commonM[commonM[i, 
                                                     1], 3], digits = 4)
        outputcommonM[i, 2] <- round((commonM[commonM[i, 
                                                      1], 3]/totalRSquare) * 100, digits = 2)
      }
      outputcommonM[totalN + 1, 1] <- round(totalRSquare, 
                                            digits = 4)
      outputcommonM[totalN + 1, 2] <- round(100, digits = 4)
      rowNames <- NULL
      for (i in 1:totalN) {
        ii <- commonM[i, 1]
        nbits <- sum(binarymx[, ii])
        cbits <- 0
        if (nbits == 1) 
          rowName <- "Unique to "
        else rowName <- "Common to "
        for (j in 1:nvar) {
          if (binarymx[j, ii] == 1) {
            if (nbits == 1) 
              rowName <- paste(rowName, iv.name[j], 
                               sep = "")
            else {
              cbits <- cbits + 1
              if (cbits == nbits) {
                rowName <- paste(rowName, "and ", sep = "")
                rowName <- paste(rowName, iv.name[j], 
                                 sep = "")
              }
              else {
                rowName <- paste(rowName, iv.name[j], 
                                 sep = "")
                rowName <- paste(rowName, ", ", sep = "")
              }
            }
          }
        }
        rowNames <- c(rowNames, rowName)
      }
      rowNames <- c(rowNames, "Total")
      rowNames <- format.default(rowNames, justify = "left")
      colNames <- format.default(c("Fractions", " % Total"), 
                                 justify = "right")
      dimnames(outputcommonM) <- list(rowNames, colNames)
      VariableImportance <- matrix(nrow = nvar, ncol = 4)
      for (i in 1:nvar) {
        VariableImportance[i, 3] <- round(sum(binarymx[i, 
        ] * (commonM[, 3]/apply(binarymx, 2, sum))), 
        digits = 4)
      }
      VariableImportance[, 1] <- outputcommonM[1:nvar, 
                                               1]
      VariableImportance[, 2] <- VariableImportance[, 
                                                    3] - VariableImportance[, 1]
      total = round(sum(VariableImportance[, 3]), digits = 4)
      VariableImportance[, 4] <- round(100 * VariableImportance[, 
                                                                3]/total, 2)
      dimnames(VariableImportance) <- list(iv.name, c("Unique", 
                                                      "Average.share", "Individual", "I.perc(%)"))
      if (commonality) {
        outputList[[k + 1]] <- outputcommonM
      }
      else {
        outputList[[k + 1]] <- VariableImportance
      }
    }
  }
  if (type == "adjR2") {
    names(outputList) <- c("adjusted.R2", r2type)
  }
  if (type == "dev") {
    names(outputList) <- c("Explained.deviance", r2type)
  }
  outputList$variables <- iv.name
  if (commonality) {
    outputList$type = "commonality.analysis"
  }
  if (!commonality) {
    outputList$type = "hierarchical.partitioning"
  }
  class(outputList) <- "gamhp"
  outputList

}
