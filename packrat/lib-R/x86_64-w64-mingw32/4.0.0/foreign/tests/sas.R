library("foreign")
(dd <- data.frame(name  = c("apple", "banana", "carrot", NA),
                  gender= c("male", "female", "male", "female"), stringsAsFactors = FALSE))
##   name gender
## 1    a   male
## 2    b female
## 3    c   male
## 4 <NA> female
setwd(tempdir())
tfSi <- "temp_for_sas_import.txt"
ti   <- "temp_import.sas"
write.foreign(dd, datafile = tfSi, codefile = ti, package = "SAS")
file.show(tfSi) # the NA is shown as <empty>
writeLines(sasCodes <- readLines(ti))
## This failed in foreign <= 0.8-71 :
stopifnot(identical(" name $ 6",
                    grep(" name ", sasCodes, value=TRUE)))

## This site was unresponsive in Jan 2014
if(!nzchar(Sys.getenv("R_FOREIGN_FULL_TEST"))) q("no")
tfile <- "int1982ag.zip"
download.file("ftp://cusk.nmfs.noaa.gov/mrfss/intercept/ag/int1982ag.zip",
              tfile, quiet=TRUE, mode="wb")
zip.file.extract("int1982ag.xpt", tfile)
dfs <- read.xport("int1982ag.xpt")
foo <- dfs$I3_19822
nrow(foo)
stopifnot(nrow(foo) == 3650)
