library(testthat)
library(mdatools)

test_check("mdatools")


#' Checks if all elements of two lists are equal
#'
#' @param l1
#' first list
#' @param l2
#' second list
testList <- function(l1, l2) {
   expect_equal(length(l1), length(l2))
   expect_equal(names(l1), names(l2))

   for (n in seq_len(length(l1))) {
      # print(names(l1)[n])
      v1 <- l1[[n]]
      v2 <- l2[[n]]

      if (is.list(v1)) testList(v1, v2)
      dim(v1) <- NULL
      dim(v2) <- NULL

      expect_equal(v1, v2)
   }
}
