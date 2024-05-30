# load_all()

strsplit1 <- function(x, split) {
  strsplit(x, split = split)[[1]]
}


(x <- "alfa,bravo,charlie,delta")

strsplit1(x, split = ",")

exists("strsplit1", where = globalenv(), inherits = FALSE)
