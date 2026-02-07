####################################
# Tests for all bugs found in 2022 #
####################################

# setup({
#    pdf(file = tempfile("mdatools-test-simca-", fileext = ".pdf"))
#    sink(tempfile("mdatools-test-test-simca-", fileext = ".txt"), append = FALSE, split = FALSE)
# })

# teardown({
#    dev.off()
#    sink()
# })

context("tests for model training:")

dc <- read.csv2("wines_baro.csv", row.names = 1)
dta <- read.csv2("wines_grigno.csv", row.names = 1)
dtf <- read.csv2("wines_full.csv", row.names = 1)
dtn <- read.csv2("wines_new.csv", row.names = 1)

Xc <- dc[, -1]
cc <- dc[, 1]

Xta <- dta[, -1]
cta <- dta[, 1]

Xtf <- dtf[, -1]
ctf <- dtf[, 1]

Xtn <- dtn

test_that("probabilities are computed correctly", {

   Xta <- mda.exclrows(Xta, 32)
   Xtn <- mda.exclrows(Xtn, 20)
   Xtf <- mda.exclrows(Xtf, c(32, 60, 120))

   m <- ddsimca(Xc, "Barolo", 10, scale = TRUE, exclrows = 32)
   rta <- predict(m, Xta, cta)
   rtf <- predict(m, Xtf, ctf)
   rtn <- predict(m, Xtn)

   # calres
   par(mfcol = c(2, 3))
   plotAcceptance(m$calres)
   plotAcceptance(m$calres, show.excluded = TRUE)
   plotAcceptance(m$calres, ncomp = 1, limType = "robust", log = TRUE)
   plotAcceptance(m$calres, ncomp = 1, limType = "robust", log = TRUE, show.excluded = TRUE)
   plotAcceptance(m$calres, ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE)
   plotAcceptance(m$calres, ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE, show.excluded = TRUE)

   # rta
   par(mfcol = c(2, 3))
   plotAcceptance(rta)
   plotAcceptance(rta, show.excluded = TRUE)
   plotAcceptance(rta, ncomp = 1, limType = "robust", log = TRUE)
   plotAcceptance(rta, ncomp = 1, limType = "robust", log = TRUE, show.excluded = TRUE)
   plotAcceptance(rta, ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE)
   plotAcceptance(rta, ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE, show.excluded = TRUE)

   # rtn
   par(mfcol = c(2, 3))
   plotAcceptance(rtn)
   plotAcceptance(rtn, show.excluded = TRUE)
   plotAcceptance(rtn, ncomp = 1, limType = "robust", log = TRUE)
   plotAcceptance(rtn, ncomp = 1, limType = "robust", log = TRUE, show.excluded = TRUE)
   plotAcceptance(rtn, ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE)
   plotAcceptance(rtn, ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE, show.excluded = TRUE)

   # rtf
   par(mfcol = c(2, 3))
   plotAcceptance(rtf)
   plotAcceptance(rtf, show.excluded = TRUE)
   plotAcceptance(rtf, ncomp = 1, limType = "robust", log = TRUE)
   plotAcceptance(rtf, ncomp = 1, limType = "robust", log = TRUE, show.excluded = TRUE)
   plotAcceptance(rtf, ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE)
   plotAcceptance(rtf, ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE, show.excluded = TRUE)

   par(mfcol = c(2, 3))
   plotAcceptance(rtf, show = "members")
   plotAcceptance(rtf, show = "members", show.excluded = TRUE)
   plotAcceptance(rtf, show = "members", ncomp = 1, limType = "robust", log = TRUE)
   plotAcceptance(rtf, show = "members", ncomp = 1, limType = "robust", log = TRUE, show.excluded = TRUE)
   plotAcceptance(rtf, show = "members", ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE)
   plotAcceptance(rtf, show = "members", ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE, show.excluded = TRUE)


   par(mfcol = c(2, 3))
   plotAcceptance(rtf, show = "strangers")
   plotAcceptance(rtf, show = "strangers", show.excluded = TRUE)
   plotAcceptance(rtf, show = "strangers", ncomp = 1, limType = "robust", log = TRUE)
   plotAcceptance(rtf, show = "strangers", ncomp = 1, limType = "robust", log = TRUE, show.excluded = TRUE)
   plotAcceptance(rtf, show = "strangers", ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE)
   plotAcceptance(rtf, show = "strangers", ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE, show.excluded = TRUE)

   par(mfcol = c(2, 3))
   plotAcceptance(rtf, show = "all")
   plotAcceptance(rtf, show = "all", show.excluded = TRUE)
   plotAcceptance(rtf, show = "all", ncomp = 1, limType = "robust", log = TRUE)
   plotAcceptance(rtf, show = "all", ncomp = 1, limType = "robust", log = TRUE, show.excluded = TRUE)
   plotAcceptance(rtf, show = "all", ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE)
   plotAcceptance(rtf, show = "all", ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE, show.excluded = TRUE)

   # as data frame
   print(as.data.frame(rta))
   print(as.data.frame(rtn))
   print(as.matrix(rta))
   print(as.matrix(rtf))
   print(as.matrix(rtn))

   writeCSV(rta, ncomp = 2, limType = "classic", fileName = "ddsimcares-alt.csv")
   writeCSV(rtf, ncomp = 2, limType = "classic", fileName = "ddsimcares-full.csv")
   writeCSV(rtn, ncomp = 2, limType = "robust", fileName = "ddsimcares-new.csv")


   par(mfcol = c(3, 4))
   plotDistances(m$calres)
   plotDistances(m$calres, ncomp = 4, distance = "q")
   plotDistances(m$calres, limType = "robust", ncomp = 2, distance = "f", log = TRUE)
   plotDistances(rta)
   plotDistances(rta, ncomp = 4, distance = "q")
   plotDistances(rta, limType = "robust", ncomp = 2, distance = "f", log = TRUE)
   plotDistances(rtf)
   plotDistances(rtf, ncomp = 4, distance = "q")
   plotDistances(rtf, limType = "robust", ncomp = 2, distance = "f", log = TRUE)
   plotDistances(rtn)
   plotDistances(rtn, ncomp = 4, distance = "q")
   plotDistances(rtn, limType = "robust", ncomp = 2, distance = "f", log = TRUE)


   m <- ddsimca(Xc, "Barolo", 10, scale = TRUE, exclrows = 32)
   rta <- predict(m, Xta, cta)
   rtf <- predict(m, Xtf, ctf)
   rtn <- predict(m, Xtn)

   writeJSON(m, "ddsimca-model.json")
})