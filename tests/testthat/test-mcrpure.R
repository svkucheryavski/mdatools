# temporary code for debugging while develop
library(mdatools)

data(simdata)
load("data/puredata.RData")
S <- simdata$S
attr(S, "yaxis.name") <- "Wavelength, nm"
attr(S, "yaxis.values") <- simdata$wavelength

# ceate concentration profiles
C1 <- dnorm(1:100, 20, 5)
C2 <- dnorm(1:100, 50, 20)
C3 <- dnorm(1:100, 70, 10)
C <- cbind(C1, C2, C3)
attr(C, "yaxis.name") <- "Time, c"
attr(C, "yaxis.values") <- (1:100)/10

D <- C %*% t(S) + matrix(rnorm(nrow(C) * nrow(S), 0, 0.001), nrow(C), nrow(S))
attr(D, "xaxis.name") <- "Wavelength, nm"
attr(D, "xaxis.values") <- simdata$wavelength
attr(D, "yaxis.name") <- "Time, c"
attr(D, "yaxis.values") <- (1:100)/10

par(mfrow = c(3, 1))
mdaplotg(mda.t(C), type = "l")
mdaplotg(mda.t(S), type = "l")
mdaplot(D, type = "l")


source("R/misc.R")
source("R/mcr.R")
source("R/mcrpure.R")


# ! test for user provided pure variables

# this should raise an error
#m <- mcrpure(D, ncomp = 3, purevars = c(50, 100, 200))

# this should be ok
m <- mcrpure(D, ncomp = 3, purevars = c(50, 100, 120))
show(m$purevars)

par(mfrow = c(2, 2))
plotPurity.mcrpure(m)
plotPurity.mcrpure(m, show.labels = TRUE)
plotPurity.mcrpure(m, type = "b", show.labels = TRUE, col = "red")

par(mfrow = c(2, 2))
plotPuritySpectra.mcrpure(m)
plotPuritySpectra.mcrpure(m, comp = 3)
plotPuritySpectra.mcrpure(m, comp = 1, type = "h", col = "black", show.lines = FALSE)
plotPuritySpectra.mcrpure(m, comp = c(2, 3))

par(mfrow = c(2, 2))
plotSpectra.mcr(m)
plotSpectra.mcr(m, comp = 3)
plotSpectra.mcr(m, comp = 1, type = "h", col = "black")
plotSpectra.mcr(m, comp = c(2, 3))

par(mfrow = c(2, 2))
plotContributions.mcr(m)
plotContributions.mcr(m, comp = 3)
plotContributions.mcr(m, comp = 1, type = "h", col = "black")
plotContributions.mcr(m, comp = c(2, 3))


# ! test for normal ones

m <- mcrpure(D, ncomp = 3)
show(m$purevars)

par(mfrow = c(2, 2))
plotPurity.mcrpure(m)
plotPurity.mcrpure(m, show.labels = TRUE)
plotPurity.mcrpure(m, type = "b", show.labels = TRUE, col = "red")

par(mfrow = c(2, 2))
plotPuritySpectra.mcrpure(m)
plotPuritySpectra.mcrpure(m, comp = 3)
plotPuritySpectra.mcrpure(m, comp = 1, type = "h", col = "black", show.lines = FALSE)
plotPuritySpectra.mcrpure(m, comp = c(2, 3))

par(mfrow = c(2, 2))
plotSpectra.mcr(m)
plotSpectra.mcr(m, comp = 3)
plotSpectra.mcr(m, comp = 1, type = "h", col = "black")
plotSpectra.mcr(m, comp = c(2, 3))

par(mfrow = c(2, 2))
plotContributions.mcr(m)
plotContributions.mcr(m, comp = 3)
plotContributions.mcr(m, comp = 1, type = "h", col = "black")
plotContributions.mcr(m, comp = c(2, 3))

