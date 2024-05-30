# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
library(devtools)

(x <- "alfa,bravo,charlie,delta")

strsplit(x, split = ",")


str(strsplit(x, split = ","))

strsplit1 <- function(x, split) {
  strsplit(x, split = split)[[1]]
}

use_r("strsplit1")
