solve_LSAP <- function (x) 
{
  storage.mode(x) <- "double"
  out <- solve_LSAP_cpp(x) + 1
  
  out
}