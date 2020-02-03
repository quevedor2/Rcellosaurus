library(devtools)

pkg=file.path("~/git", "Rcellosaurus")

#### Assembling data ####
#usethis::use_data_raw()

#### Building ####
setwd(pkg)
# create("CCLid")
devtools::document(pkg)
devtools::check(pkg)
devtools::build(pkg)

