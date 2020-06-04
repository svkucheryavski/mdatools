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

data(simdata)
D <- simdata$spectra.c

cc <- list(constraint("non-negativity"))
cs <- list(constraint("non-negativity"))

m <- mcrals(D, ncomp = 3, cont.constraints = cc, spec.constraints = cs)
print(m)