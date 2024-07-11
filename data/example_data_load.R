library("utils")
example_data = utils::read.csv("example_data.csv",row.names = 1)
usethis::use_data(example_data, overwrite = TRUE)
