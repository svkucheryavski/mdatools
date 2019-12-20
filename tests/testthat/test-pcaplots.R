#####################################################
# Tests for basic functionality of plot() class  #
#####################################################

pdf(file = "test_pca_plots.pdf")

# prepare cases

## 1. full data no test set
data(people)
x1 <- people
m1 <- pca(x1, 10, scale = T)

## 2. calibration and test sets
ind <- seq(1, 32, by = 4)
x2 <- people[-ind, ]
x2.test <- people[ind, ]
m2 <- pca(x2, 10, scale = T, x.test = x2.test)

## 3. calibration and test sets with excluded data
x3 <- people[-ind, ]
x3 <- mda.exclrows(x3, c(1, 10, 20))
x3 <- mda.exclcols(x3, c(3, 12))
x3.test <- people[ind, ]
m3 <- pca(x3, 10, scale = T, x.test = x3.test)

## combine all together
x <- list(x1, x2, x3)
m <- list("full" = m1, "full + test" = m2, "excluded + test" = m3)


#########################################
# Block 1: Variance plot                  #
#########################################

context("pca: variance plots")

tf <- function(model, name) {

   par(mfrow = c(1, 1))
   plot.new()
   text(0, 0, paste0("PCA - variance plot - ", name), pos = 4)

   # basic variance plots
   expect_silent({
      par(mfrow = c(2, 2))
      plotVariance(model)
      plotVariance(model, type = "h", show.labels = T)
      plotVariance(model, show.labels = T, col = c("red", "green"))
      plotVariance(model, res = list("cal" = model$calres))
   })

   # basic cumulative plots
   expect_silent({
      par(mfrow = c(2, 2))
      plotCumVariance(model)
      plotCumVariance(model, type = "h", show.labels = T)
      plotCumVariance(model, show.labels = T, col = c("red", "green"))
      plotCumVariance(model, res = list("cal" = model$calres))
   })

}

for (i in seq_len(length(m))) {
   tf(m[[i]], names(m)[[i]])
}

dev.off()
stop()

#########################################
# Block 2: Scores plot                  #
#########################################

context("pca: scores plots")

tf <- function(model, x.cal) {

   # basic plots (scatter)
   expect_silent({
      par(mfrow = c(2, 2))
      plotScores(model)
      plotScores(model, c(1, 3))
      plotScores(model, c(1, 3), show.labels = T, col = c("red", "green"))
      plotScores(model, c(1, 3), res = list("cal" = model$calres))
   })

   # basic plots (scatter, one components)
   expect_silent({
      par(mfrow = c(2, 2))
      plotScores(model, 1)
      plotScores(model, 2)
      plotScores(model, 2, show.labels = T, col = c("red", "green"))
      plotScores(model, 2, res = list("cal" = model$calres))
   })

   # basic plots (line plot, two components)
   expect_silent({
      par(mfrow = c(2, 2))
      plotScores(model, type = "l")
      plotScores(model, c(1, 3), type = "l")
      plotScores(model, c(1, 3), type = "l", show.labels = T, col = c("red", "green"))
      plotScores(model, c(1, 3), type = "l", res = list("cal" = model$calres))
   })

   # playing with legend
   expect_silent({
      par(mfrow = c(2, 2))
      plotScores(model, show.legend = FALSE)
      plotScores(model, c(1, 3), legend.position = "top")
      plotScores(model, c(1, 3), show.labels = T, legend.position = "topleft")
      plotScores(model, c(1, 3), legend.position = "bottom", res = list("cal" = model$calres))
   })

   # show hidden values
   expect_silent({
      par(mfrow = c(2, 2))
      plotScores(model, show.excluded = TRUE)
      plotScores(model, c(1, 3), show.excluded = TRUE)
      plotScores(model, c(1, 3), show.excluded = TRUE, show.labels = T, col = "red")
      plotScores(model, c(1, 3), show.excluded = TRUE, res = list("cal" = model$calres))
   })

   # show hidden values

}

for (i in seq_len(length(m))) {
   tf(m[[i]], x[[i]])
}

dev.off()

