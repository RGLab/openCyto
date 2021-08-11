fast_rlm <-
  function(x, y, weights, ...)
{
    w = rep(1, nrow(x))
    psi = psi.huber
    k2 = 1.345
    maxit = 20
    acc = 1e-4
    irls.delta <- function(old, new)
        sqrt(sum((old - new)^2)/max(1e-20, sum(old^2)))
   
    if(qr(x)$rank < ncol(x))
        stop("'x' is singular: singular fits are not implemented in 'rlm'")

    
    temp <-  lm.wfit(x, y, w, method="qr")
     coef <- temp$coefficients
    resid <- temp$residuals
   
    
    done <- FALSE
    conv <- NULL
    n1 <- nrow(x) - ncol(x)
    theta <- 2*pnorm(k2) - 1
    gamma <- theta + k2^2 * (1 - theta) - 2 * k2 * dnorm(k2)
    ## At this point the residuals are weighted for inv.var and
    ## unweighted for case weights.  Only Huber handles case weights
    ## correctly.
    scale <- mad(resid, 0)
    for(iiter in 1L:maxit) {
        testpv <- resid
        scale <- median(abs(resid))/0.6745
        if(scale == 0) {
            done <- TRUE
            break
        }
        
        w <- psi(resid/scale)
        temp <- lm.wfit(x, y, w, method="qr")
        coef <- temp$coefficients
        resid <- temp$residuals
        convi <- irls.delta(testpv, resid)
        conv <- c(conv, convi)
        done <- (convi <= acc)
        if(done) break
    }
    if(!done)
        warning(gettextf("'rlm' failed to converge in %d steps", maxit),
                domain = NA)
    fitted <- drop(x %*% coef)
    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl <- match.call()
    cl[[1L]] <- as.name("rlm")
    fit <- list(coefficients = coef, residuals = y - fitted, wresid = resid,
                effects = temp$effects,
                rank = temp$rank, fitted.values = fitted,
                assign = temp$assign,  qr = temp$qr, df.residual = NA, w = w,
                s = scale, psi = psi, k2 = k2,
                weights = NULL,
                conv = conv, converged = done, x = x, call = cl)
    class(fit) <- c("rlm", "lm")
    fit
}

psi.huber <- function(u, k = 1.345, deriv=0)
{
  if(!deriv) return(pmin(1, k / abs(u)))
  abs(u) <= k
}