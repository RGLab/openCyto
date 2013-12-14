unlockNamespace <- function(env) {
  .Call('unlockNamespace', PACKAGE = 'openCyto', env)
}
