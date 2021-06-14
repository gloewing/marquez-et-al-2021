###############
# GLMs
###############

robustSE <- function(fit, digits=3){  
    # for GLMs
    Xmat <- model.matrix(terms(fit), model.frame(fit)) 
    Umat <- residuals(fit, type="working") * fit$weights * Xmat
    modelV  <- summary(fit)$cov.unscaled  
    robustV <- modelV %*% t(Umat) %*% Umat %*% modelV
    value <- cbind(fit$coef, sqrt(diag(modelV)), sqrt(diag(robustV)),
                   sqrt(diag(robustV))/sqrt(diag(modelV)))
    colnames(value) <- c("Estimate", "Model SE", "Robust SE", "  Ratio")
    
    return(round(value, digits=digits))
}

robustVCov <- function(fit, digits=3){  
     # for GLMs
    Xmat <- model.matrix(terms(fit), model.frame(fit)) 
    Umat <- residuals(fit, type="working") * fit$weights * Xmat
    modelV  <- summary(fit)$cov.unscaled  
    robustV <- modelV %*% t(Umat) %*% Umat %*% modelV
    value <- cbind(fit$coef, sqrt(diag(modelV)), sqrt(diag(robustV)),
                   sqrt(diag(robustV))/sqrt(diag(modelV)))
    colnames(value) <- c("Estimate", "Model SE", "Robust SE", "  Ratio")
    
    return(robustV)
}

robustCIglm <- function(fit, alpha=0.05, digits=3, B=NULL, expo=TRUE)
{
    ## for GLMs
    covNaive <- vcov(fit) # also works
    seNaive  <- sqrt(diag(covNaive))
    ##
    X      <- model.matrix(fit)
    epsHat <- fit$residuals
    covRobust  <- robustVCov(fit, digits = digits)
    seRobust   <- sqrt(diag(covRobust))
    ##
    value <- cbind(coef(fit),
                   coef(fit) - qnorm(1-alpha/2)*seNaive,
                   coef(fit) + qnorm(1-alpha/2)*seNaive,
                   coef(fit) - qnorm(1-alpha/2)*seRobust,
                   coef(fit) + qnorm(1-alpha/2)*seRobust)
    if(expo){
        rowNames <- c("exp{betaHat}", "Naive Lo", "Naive Up", "Robust Lo", "Robust Up")
        value <- exp(value)
    }else{
        rowNames <- c("betaHat", "Naive Lo", "Naive Up", "Robust Lo", "Robust Up")
    }
    
    ##
    if(!is.null(B))
    {
        Y <- model.frame(fit)[,1]
        betaBoot <- matrix(NA, nrow=B, ncol=length(coef(fit)))
        for(b in 1:B)
        {
            bootSample   <- sample(1:length(Y), replace=TRUE)
            Yb           <- Y[bootSample]
            Xb           <- X[bootSample,]
            betaBoot[b,] <- solve(t(Xb)%*%Xb) %*% t(Xb) %*% Yb
        }
        covBoot <- cov(betaBoot)
        seBoot  <- sqrt(diag(covBoot))
        ##
        value <- cbind(value,
                       coef(fit) - qnorm(1-alpha/2)*seBoot,
                       coef(fit) + qnorm(1-alpha/2)*seBoot)
        ##
        rowNames <- c(rowNames, "Boot Lo", "Boot Up")
    }
    ##
    dimnames(value) <- list(names(coef(fit)), rowNames)
    return(round(value, digits=digits))
} 

##
expit <- function(x) exp(x) / (1 + exp(x))
logit <- function(p) log(p/(1-p))


###########
# Either
###########
getCI <- function(fit, alpha=0.05, digits=2, expo=TRUE)
{
    se <- sqrt(diag(summary(fit)$cov.scaled))
    value <- cbind(coef(fit),
                   coef(fit) - qnorm(1-alpha/2)*se,
                   coef(fit) + qnorm(1-alpha/2)*se)
    colNames <- c("beta", "lower", "upper")
    if(expo == TRUE)
    {
        value <- exp(value)
        colNames <- c("exp{beta}", "lower", "upper")
    }
    dimnames(value) <- list(names(coef(fit)), colNames)
    return(round(value, digits=digits))
}



################################
# for Linear Models (NOT GLMs)
################################



##
robustCI <- function(fit, alpha=0.05, digits=2, B=NULL)
{
	##
	covNaive <- summary(fit)$cov * summary(fit)$sigma^2 # vcov(fit) # also works
	seNaive  <- sqrt(diag(covNaive))
	##
	X      <- model.matrix(fit)
	epsHat <- fit$residuals
	covRobust  <- solve(t(X)%*%X) %*% (t(X)%*%diag(epsHat^2)%*%X) %*% solve(t(X)%*%X)
	seRobust   <- sqrt(diag(covRobust))
	##
	value <- cbind(coef(fit),
								 coef(fit) - qnorm(1-alpha/2)*seNaive,
								 coef(fit) + qnorm(1-alpha/2)*seNaive,
								 coef(fit) - qnorm(1-alpha/2)*seRobust,
								 coef(fit) + qnorm(1-alpha/2)*seRobust)
	rowNames <- c("betaHat", "Naive Lo", "Naive Up", "Robust Lo", "Robust Up")
	##
	if(!is.null(B))
	{
		Y <- model.frame(fit)[,1]
		betaBoot <- matrix(NA, nrow=B, ncol=length(coef(fit)))
		for(b in 1:B)
		{
			bootSample   <- sample(1:length(Y), replace=TRUE)
			Yb           <- Y[bootSample]
			Xb           <- X[bootSample,]
			betaBoot[b,] <- solve(t(Xb)%*%Xb) %*% t(Xb) %*% Yb
		}
		covBoot <- cov(betaBoot)
		seBoot  <- sqrt(diag(covBoot))
		##
		value <- cbind(value,
								   coef(fit) - qnorm(1-alpha/2)*seBoot,
								   coef(fit) + qnorm(1-alpha/2)*seBoot)
		##
		rowNames <- c(rowNames, "Boot Lo", "Boot Up")
	}
	##
	dimnames(value) <- list(names(coef(fit)), rowNames)
	return(round(value, digits=digits))
} 


##
waldTest <- function(fit, vec, digits=c(2, 4))
{
	beta     <- coef(fit)[vec]
	varMat   <- summary(fit)$cov.scaled[vec,vec]
	testStat <- t(beta) %*% solve(varMat) %*% beta
	pVal     <- 1 - pchisq(testStat, length(vec))
	return(c(round(testStat, digits=digits[1]), round(pVal, digits=digits[2])))
}

##
LRtest <- function(fit0, fit1, digits=2)
{
	testStat <- abs(deviance(fit0) - deviance(fit1))
	testDF   <- abs(length(fit0$coef) - length(fit1$coef))
	pVal     <- 1 - pchisq(testStat, testDF)
	cat("Test Statistic =", round(testStat, 1), "on", testDF, "df => p-value =", round(pVal, digits=digits), "\n")
	return(round(pVal, digits=digits))
}