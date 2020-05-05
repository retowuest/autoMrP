.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Depending on the number of context-level variables, the set of tuning
                        parameters and processing power of the computer, predictions may take
                        some time. Using the example data, 6 context-level variables and all
                        5 classifiers, it takes 30 minutes on average.")
}
