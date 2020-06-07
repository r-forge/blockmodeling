## code to prepare `DATASET` dataset goes here

usethis::use_data("DATASET")
baker <-
  read.delim(file = "baker.txt",
             header = TRUE,
             row.names = 1)
# Transforming it to matrix format
baker <- as.matrix(baker)
#putting zeros on the diagonal
diag(baker) <- 0 
devtools::use_data(baker, overwrite = TRUE)
devtools::load_all()
