fast_rlm <-
  function(x, y, maxit = 20, ...)
{
    
   
    temp <- rlm_cpp(x, y, maxit)
    fitted <- drop(x %*% temp$coefficients)
    fit <- list(coefficients = temp$coefficients
                , residuals = y - fitted
                , wresid = temp$residuals,
                effects = temp$effects,
                rank = temp$rank
                , fitted.values = fitted,
                assign = temp$assign,  qr = temp$qr
                , df.residual = NA
                , w = temp$w
                ,s = temp$scale
                # , psi = psi
                # , k2 = k2
                ,weights = NULL
                # ,conv = conv
                , converged = temp$done, x = x
                # , call = cl
                )
    class(fit) <- c("rlm", "lm")
    fit
}
