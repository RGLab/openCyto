#' rewrite huber estimator
robust_m_estimator <- function(x, sd)
{
  stopifnot(sd != 0) 
  Winsorizes_k = 1.5
  convergence_tol = 1e-6
  x <- x[!is.na(x)]
  n <- length(x)
  mu <- median(x)
  Winsorizes_k <- sd * Winsorizes_k
  convergence_tol <- sd * convergence_tol
  while(TRUE)
  {
    p <- mu - Winsorizes_k
    y <- sapply(x, function(i)max(p, i))
    q <- mu + Winsorizes_k
    y <- sapply(y, function(i)min(q, i))
    res <- mean(y)
    if(abs(mu - res) < convergence_tol)
      return(mu)
    mu = res
  }
  
}