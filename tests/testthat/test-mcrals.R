########################################################
# Block 1: testing methods implementing pca algorithms #
########################################################

#setup({
#   #pdf(file = "mdatools-test-mcrpure.pdf")
#   pdf(file = tempfile("mdatools-test-mcrpure-", fileext = ".pdf"))
#   sink(tempfile("mdatools-test-mcrpure-", fileext = ".txt"), append = FALSE, split = FALSE)
#})

#teardown({
#   dev.off()
#   sink()
#})
#

data(carbs)
D <- carbs$D
S <- carbs$S
C <- carbs$C

cc <- list(
   constraint("non-negativity")
)

cs <- list(
   constraint("non-negativity"),
   constraint("norm", params = list(type = "area"))
)

set.seed(6)
m <- mcrals(D, ncomp = 3, cont.constraints = cc, spec.constraints = cs, verbose = TRUE)
summary(m)

par(mfrow = c(2, 3))

# resolved spectra vs new
for (i in 1:min(c(m$ncomp, ncol(S)))) {
   mdaplotg(
      list(
         original = prep.norm(mda.subset(mda.t(S), i), "area"),
         resolved = prep.norm(mda.subset(mda.t(m$resspec), i), "area")
      ), type = "l", col = c("black", "red"), lwd = c(3, 1), opacity = c(0.75, 1),
      , xlim = c(1600, 200), xticks = seq(1600, 200, by = -200)
   )
}

# resolved contributions vs new
for (i in 1:min(c(m$ncomp, ncol(C)))) {
   mdaplotg(
      list(
         original = prep.norm(mda.subset(mda.t(C), i), "area"),
         resolved = prep.norm(mda.subset(mda.t(m$rescont), i), "area")
      ), type = "l", col = c("black", "red"), lwd = c(3, 1), opacity = c(0.75, 1)
   )
}

# check without constraints
# check default seeints verbode = FALSE
# check default seeints verbode = TRUE
