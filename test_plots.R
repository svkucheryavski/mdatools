rm(list = ls())
roxygen2::roxygenize()
pdf(file = "~/Desktop/test_mdaplots.pdf")
devtools::test_file("tests/testthat/test-getcolors.R")
devtools::test_file("tests/testthat/test-plotseries.R")
devtools::test_file("tests/testthat/test-mdaplot1.R")
devtools::test_file("tests/testthat/test-mdaplot2.R")
devtools::test_file("tests/testthat/test-mdaplotg.R")
dev.off()
options(device = "quartz")


data(simdata)
x = simdata$spectra.c
y = simdata$conc.c[, 1]

model = pls(x, y, ncomp = 8, cv = 1)
summary(model)
plot.new()
plot(model)

plotRegcoeffs(model)

plot(model$coeffs)