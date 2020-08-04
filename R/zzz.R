.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Predictions may take some time. Use multiple cores for a substantial speed increase.")
}
