#' robust linear model using an M estimator 
#' 
#' rewritten in c++, till eval stats::lm.wfit r function in underlying cpp11 code 
#' It is internally used for singletGate, thus its output format may not be generic enough for common model fitting .
#' e.g. it doesn't take formula as input
#' 
#' @param x matrix with first column as weight (default can be 1s), the rest columns are predict variable
#' @param y numeric vector as response 
#' @param maxit maximum iterations
#' @examples 
#' @noRd
#' n <- 1e3
#' x <- seq_len(n)
#' y <- x * 2.5 - 1.3 + rnorm(n, sd = 30)
#' names(y) <- x
#' x <- cbind(1, x)
#' r2 <- fast_rlm(x, y)
fast_rlm <- function(x, y, maxit = 20)
{
    
    fit_result <- rlm_cpp(x, y, maxit)
    fitted <- drop(x %*% fit_result$coefficients)
    fit <- list(coefficients = fit_result$coefficients
                , residuals = y - fitted
                , wresid = fit_result$residuals,
                effects = fit_result$effects,
                rank = fit_result$rank
                , fitted.values = fitted,
                assign = fit_result$assign,  qr = fit_result$qr
                , df.residual = NA
                , w = fit_result$w
                ,s = fit_result$scale
                ,weights = NULL
                , converged = fit_result$done, x = x
                )
    class(fit) <- c("rlm", "lm")
    fit
}
