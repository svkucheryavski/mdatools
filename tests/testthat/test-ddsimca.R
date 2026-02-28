####################################
# Tests for all bugs found in 2022 #
####################################

setup({
   pdf(file = "dump/mdatools-test-dsimca.pdf")
   sink("dump/mdatools-test-test-simca-.txt", append = FALSE, split = FALSE)
})

teardown({
   dev.off()
   sink()
})

compare_csv_files <- function(file1, file2, sep = ",", tolerance = 1e-3) {
  # Determine decimal separator based on field separator
  decimal_sep <- if (sep == ",") "." else ","

  # Read files as lines
  lines1 <- readLines(file1, warn = FALSE)
  lines2 <- readLines(file2, warn = FALSE)

  # Remove empty lines (trim whitespace first)
  lines1 <- lines1[trimws(lines1) != ""]
  lines2 <- lines2[trimws(lines2) != ""]

  # Check if number of lines match
  if (length(lines1) != length(lines2)) {
    stop(sprintf("Different number of non-empty lines: %d vs %d",
                 length(lines1), length(lines2)))
  }

  # Compare line by line
  differences <- list()

  for (i in seq_along(lines1)) {
    # Split by separator
    vals1 <- strsplit(lines1[i], sep, fixed = TRUE)[[1]]
    vals2 <- strsplit(lines2[i], sep, fixed = TRUE)[[1]]

    # Trim whitespace
    vals1 <- trimws(vals1)
    vals2 <- trimws(vals2)

    # Check if same number of values
    if (length(vals1) != length(vals2)) {
      differences[[paste0("Line_", i)]] <- "Different number of values"
      next
    }

    # Compare each value
    line_differs <- FALSE
    for (j in seq_along(vals1)) {
      # Prepare for numeric conversion
      num_val1 <- vals1[j]
      num_val2 <- vals2[j]

      # If decimal separator is comma, replace with dot for conversion
      if (decimal_sep == ",") {
        num_val1 <- gsub(",", ".", num_val1, fixed = TRUE)
        num_val2 <- gsub(",", ".", num_val2, fixed = TRUE)
      }

      # Try to convert to numeric
      num1 <- suppressWarnings(as.numeric(num_val1))
      num2 <- suppressWarnings(as.numeric(num_val2))

      # If both are numeric, compare with tolerance
      if (!is.na(num1) && !is.na(num2)) {
         if (abs(num1 - num2) > tolerance * max(abs(num1), abs(num2))) {
          line_differs <- TRUE
          break
        }
      } else {
        # Otherwise compare as strings
        if (vals1[j] != vals2[j]) {
          line_differs <- TRUE
          break
        }
      }
    }

    if (line_differs) {
      differences[[paste0("Line_", i)]] <- list(
        file1 = lines1[i],
        file2 = lines2[i]
      )
    }
  }

  # Return results
  list(
    equal = length(differences) == 0,
    total_lines = length(lines1),
    separator = sep,
    decimal_separator = decimal_sep,
    differences = if (length(differences) > 0) differences else NULL
  )
}

compareModels <- function(mr, mjs) {

   if (is.null(attr(mr$loadings, "yaxis.values"))) {
      attr(mjs$loadings, "yaxis.values") <- NULL
   }

   expect_equal(abs(mjs$loadings), abs(mda.purge(mr$loadings)), tolerance = 1e-4)
   expect_equal(mr$eigenvals, mjs$eigenvals)

   # we need to round DoF in R as they are rounded in web version
   mr$limParams$Q$moments$Nu <- round(mr$limParams$Q$moments$Nu)
   mr$limParams$Q$robust$Nu <- round(mr$limParams$Q$robust$Nu)
   mr$limParams$T2$moments$Nu <- round(mr$limParams$T2$moments$Nu)
   mr$limParams$T2$robust$Nu <- round(mr$limParams$T2$robust$Nu)
   testList(mr$limParams, mjs$limParams)

   # same for Qlim/T2lim
   expect_equal(mr$Qlim, mjs$Qlim)
   expect_equal(mr$T2lim, mjs$T2lim)

   # DDSIMCA related settings
   expect_equal(mr$alpha, mjs$alpha)
   expect_equal(mr$gamma, mjs$gamma)
   expect_equal(mr$classname, mjs$classname)
   expect_equal(mr$ncomp.selected, mjs$ncomp.selected)
   expect_equal(mr$nclasses, mjs$nclasses)
   expect_equal(mr$nrows, mjs$nrows)

   # prep
   ncr <- ncol(mr$loadings)
   ncjs <- ncol(mjs$loadings)
   if (!is.null(mr$prep) && ncr == ncjs) {
      for (n in seq_along(mr$prep)) {
         if (!is.list(mr$prep[[n]])) next
         testList(mr$prep[[n]]$params, mjs$prep[[n]]$params)
      }
   }
}

context("tests for model training:")

dc <- read.csv2("datasets/wines_baro.csv", row.names = 1, check.names = FALSE)
dta <- read.csv2("datasets/wines_grigno.csv", row.names = 1)
dtf <- read.csv2("datasets/wines_full.csv", row.names = 1)
dtn <- read.csv2("datasets/wines_new.csv", row.names = 1)

src <- "https://mda.tools/ddsimca/"

check_case <- function(m, Xta, cta, Xtf, ctf, Xtn, limType = "classic", caseNum = 1) {

   print(paste0("Testing case ", caseNum))


   ra <- predict(m, Xta, cta)
   rf <- predict(m, Xtf, ctf)
   rn <- predict(m, Xtn)
   rc <- m$res$cal

   train.file <- sprintf('ddsimca-case%d-train.csv', caseNum)
   alt.file <- sprintf('ddsimca-case%d-alt.csv', caseNum)
   full.file <- sprintf('ddsimca-case%d-full.csv', caseNum)
   new.file <- sprintf('ddsimca-case%d-new.csv', caseNum)

   writeCSV(rc, paste0('dump/', train.file), limType = limType, src = src, name = 'train', sep = ';', dataFile = 'wines_baro')
   writeCSV(ra, paste0('dump/', alt.file), limType = limType, src = src, name = 'test', sep = ';', dataFile = 'wines_grigno')
   writeCSV(rf, paste0('dump/', full.file), limType = limType, src = src, name = 'test', sep = ';', dataFile = 'wines_full')
   writeCSV(rn, paste0('dump/', new.file), limType = limType, src = src, name = 'new', sep = ';', dataFile = 'wines_new')

   res1 <- compare_csv_files(paste0('dump/', train.file), paste0('auxfiles/', train.file), sep = ";")
   expect_null(res1$differences)

   res2 <- compare_csv_files(paste0('dump/', alt.file), paste0('auxfiles/', alt.file), sep = ";")
   expect_null(res2$differences)

   res3 <- compare_csv_files(paste0('dump/', full.file), paste0('auxfiles/', full.file), sep = ";")
   expect_null(res3$differences)

   res4 <- compare_csv_files(paste0('dump/', new.file), paste0('auxfiles/', new.file), sep = ";")
   expect_null(res4$differences)

   model.file <- sprintf("ddsimca-model-%d.json", caseNum)
   writeJSON(m, paste0("./dump/", model.file))
   m1 <- readJSON(paste0("./dump/", model.file))
   m2 <- readJSON(paste0("./jsonfiles/", model.file))
   compareModels(m, m1)
   compareModels(m, m2)

   par(mfcol = c(1, 1))
   plot(0, type = "n", xlim = c(0, 1), ylim = c(0, 1))
   text(0.5, 0.5, paste0("Case ", caseNum), cex = 3)

   par(mfcol = c(2, 3))
   plotAcceptance(rc)
   plotAcceptance(rc, pch = 5, show = "all")
   plotAcceptance(rc, ncomp = 5, log = TRUE, limType = "robust", pch = 5)
   plotAcceptance(rc, ncomp = 5, log = TRUE, limType = "robust", pch = 5, show = "all", show.labels = TRUE)
   plotAcceptance(rc, main = "training set", pch = 5, show.excluded = TRUE, show.labels = TRUE)
   plotAcceptance(rc, main = "training set", pch = 5, show = "all", show.excluded = TRUE, show.labels = TRUE)

   par(mfcol = c(2, 3))
   plotAcceptance(m)
   plotAcceptance(m, pch = 5, show = "all")
   plotAcceptance(m, ncomp = 5, log = TRUE, limType = "robust", pch = 5)
   plotAcceptance(m, ncomp = 5, log = TRUE, limType = "robust", pch = 5, show = "all", show.labels = TRUE)
   plotAcceptance(m, main = "training set", pch = 5, show.excluded = TRUE, show.labels = TRUE)
   plotAcceptance(m, main = "training set", pch = 5, show = "all", show.excluded = TRUE, show.labels = TRUE)

   par(mfcol = c(2, 3))
   plotAcceptance(m, res = "pv", show.labels = TRUE)
   plotAcceptance(m, res = "pv", pch = 5, show = "all")
   plotAcceptance(m, res = "pv", ncomp = 5, log = TRUE, limType = "robust", pch = 5)
   plotAcceptance(m, res = "pv", ncomp = 5, log = TRUE, limType = "robust", pch = 5, show = "all", show.labels = TRUE)
   plotAcceptance(m, res = "pv", main = "training set", pch = 5, show.excluded = TRUE, show.labels = TRUE)
   plotAcceptance(m, res = "pv", main = "training set", pch = 5, show = "all", show.excluded = TRUE, show.labels = TRUE)

   par(mfrow = c(2, 3))
   plotAcceptance(rf)
   plotAcceptance(rf, pch = 5, res.name = "test full", show = "strangers")
   plotAcceptance(rf, pch = 5, show = "all")

   plotAcceptance(rf, log = TRUE, pch = 5)
   plotAcceptance(rf, log = TRUE, pch = 5, show = "strangers")
   plotAcceptance(rf, log = TRUE, pch = 5, show = "all")

   par(mfrow = c(2, 3))
   plotAcceptance(ra)
   plotAcceptance(ra, log = TRUE, pch = 5, show = "strangers", res.name = "test alt")
   plotAcceptance(ra, log = TRUE, ncomp = 10, limType = "robust", pch = 5)
   plotAcceptance(rn)
   plotAcceptance(rn, ncomp = 5, res.name = "new")
   plotAcceptance(rn, show.labels = TRUE, show.excluded = TRUE, log = TRUE, pch = 5, show = "all")


   par(mfcol = c(3, 2))
   plotDistances(rc)
   plotDistances(rc, distance = "q", log = TRUE, show.excluded = TRUE, show.labels = TRUE)
   plotDistances(rc, distance = "f", ncomp = 5, show.labels = TRUE)

   plotDistances(rf)
   plotDistances(rf, distance = "q", log = TRUE, show.excluded = TRUE, show.labels = TRUE)
   plotDistances(rf, distance = "f", log = TRUE, ncomp = 5, show.labels = TRUE)

   par(mfcol = c(3, 2))
   plotDistances(m)
   plotDistances(m, distance = "q", log = TRUE, show.excluded = TRUE, show.labels = TRUE)
   plotDistances(m, distance = "f", ncomp = 5, show.labels = TRUE)

   plotDistances(m, res = "cal")
   plotDistances(m, res = "pv", distance = "q", log = TRUE, show.excluded = TRUE, show.labels = TRUE)
   plotDistances(m, res = "diff", col = "darkgray", distance = "f", ncomp = 5, show.labels = TRUE)

   plotDistDoF(m)

   # extremes plot
   expect_error(plotExtreme(rn))
   expect_error(plotExtreme(ra))
   expect_error(plotExtreme(m, col = "blue"))
   expect_error(plotExtreme(m, pch = 1))

   par(mfrow = c(2, 3))
   plotExtremes(rc)
   plotExtremes(rf)
   plotExtremes(rf, ncomp = 5, limType = "robust", col = "red", pch = 2, cex = 0.8)

   plotExtremes(m)
   plotExtremes(m, ncomp = 5, limType = "robust")
   plotExtremes(m, col = c("green", "pink"), pch = c(21, 22), bg = "yellow", cex = 0.8)

   # sensitivity and other FomS
   expect_error(plotSensitivity(rn))
   expect_error(plotSensitivity(ra))

   par(mfrow = c(2, 3))
   plotSensitivity(rc)
   plotSensitivity(rc, ylim = c(0.8, 1.1), show.ci = TRUE)
   plotSensitivity(rc, ylim = c(0.8, 1.1), limType = "robust", col = "red", pch = 1)
   plotSensitivity(rf)
   plotSensitivity(rf, ylim = c(0.8, 1.1), show.ci = TRUE)
   plotSensitivity(rf, ylim = c(0.8, 1.1), limType = "robust", col = "red", pch = 1)

   par(mfrow = c(2, 3))
   plotSensitivity(m)
   plotSensitivity(m, ylim = c(0.8, 1.1), show.ci = TRUE, ci.col = "#e0e0e0", ci.lty = c(3, 1, 3))
   plotSensitivity(m, ylim = c(0.7, 1.1), limType = "robust", legend.position = "topleft", col = c("green", "orange"))

   plotLoadings(m)
   plotScores(m, pch = 1)

   par(mfrow = c(2, 4))
   plotFoM(rf, "sens")
   plotFoM(rf, "spec")
   plotFoM(rf, "eff")
   plotFoM(rf, "sel", col = "red", pch = 3, ylim = c(0.5, 1.1))

   plotFoM(rc, "sens")
   plotFoM(ra, "spec")
   plotFoM(ra, "sel")
   plotFoM(rf, "acc")

   plotFoMs(rc)
   plotFoMs(ra)
   plotFoMs(rf, show.labels = TRUE, labels = "values", type = "h")
   plotFoMs(rf, c("sens", "spec", "eff"), legend.position = "top")


   # selectivity area

   expect_error(plotSelectivityArea(rc))
   expect_error(plotSelectivityArea(rn))

   par(mfrow = c(2, 2))
   plotSelectivityArea(rf)
   plotSelectivityArea(rf, ncomp = 5, limType = "robust")
   plotSelectivityArea(ra, ncomp = 1, col = "blue")
   plotSelectivityArea(ra, ncomp = 5, col = "blue")

   # aliens plot
   expect_error(plotAliens(rn))
   expect_error(plotAliens(rc))

   par(mfrow = c(2, 3))
   plotAliens(ra)
   plotAliens(ra, ncomp = 5)
   plotAliens(ra, ncomp = 5, limType = "robust", col = "red", pch = 2, cex = 0.8)
   plotAliens(rf)
   plotAliens(rf, ncomp = 5)
   plotAliens(rf, ncomp = 5, limType = "robust", col = "red", pch = 2, cex = 0.8)

   # plot summary
   plot(m)
   plot(rc)
   plot(rf)
   plot(ra)
   plot(rn)

   # print methods
   summary(m)
   summary(rc)
   summary(rf, limType = "robust", ncomp = 5)
   summary(ra)
   summary(rn)

   print(m)
   print(rc)
   print(rf)
   print(rn)
   print(ra)

   print(as.matrix(rc))
   print(as.matrix(rf))
   print(as.matrix(rn))
   print(as.matrix(ra))

   print(as.data.frame(rc))
   print(as.data.frame(rf))
   print(as.data.frame(rn))
   print(as.data.frame(ra))
}

test_that("case1: no prep, no outliers, no exclvars, A = 6, default settings", {

   Xc <- dc[, -1]
   Xta <- dta[, -1]
   Xtf <- dtf[, -1]
   Xtn <- dtn

   cc <- dc[, 1]
   cta <- dta[, 1]
   ctf <- dtf[, 1]

   m <- ddsimca(Xc, "Barolo", 10, scale = TRUE)
   m <- selectCompNum(m, 6)

   check_case(m, Xta, cta, Xtf, ctf, Xtn, caseNum = 1)
});


test_that("case2: no prep, no outliers, no exclvars, A = 4, limType = 'robust', alpha = 1, gamma = 0.1", {

   Xc <- dc[, -1]
   Xta <- dta[, -1]
   Xtf <- dtf[, -1]
   Xtn <- dtn

   cc <- dc[, 1]
   cta <- dta[, 1]
   ctf <- dtf[, 1]

   m <- ddsimca(Xc, "Barolo", 10, scale = TRUE)
   m <- selectCompNum(m, 4)
   m <- setParams(m, alpha = 0.01, gamma = 0.001)
   check_case(m, Xta, cta, Xtf, ctf, Xtn, limType = "robust", caseNum = 2)
});

test_that("case3: no prep, with outliers, no exclvars, A = 4, limType = 'robust', alpha = 1, gamma = 0.1", {

   Xc <- dc[, -1]
   Xta <- dta[, -1]
   Xtf <- dtf[, -1]
   Xtn <- dtn

   cc <- dc[, 1]
   cta <- dta[, 1]
   ctf <- dtf[, 1]

   m <- ddsimca(Xc, "Barolo", 10, scale = TRUE, exclrows = c(32, 40))
   m <- selectCompNum(m, 2)
   m <- setParams(m, alpha = 0.01, gamma = 0.01)
   check_case(m, Xta, cta, Xtf, ctf, Xtn, limType = "classic", caseNum = 3)
});

test_that("case4: with prep, with outliers, no exclvars, A = 4, limType = 'robust', alpha = 1, gamma = 0.1", {

   Xc <- dc[, -1]
   Xta <- dta[, -1]
   Xtf <- dtf[, -1]
   Xtn <- dtn

   cc <- dc[, 1]
   cta <- dta[, 1]
   ctf <- dtf[, 1]

   p <- list(
      prep("center", type = "median"),
      prep("scale", type = "pareto")
   )
   m <- ddsimca(Xc, "Barolo", 10, scale = FALSE, exclrows = c(32, 40), prep = p)
   m <- selectCompNum(m, 2)
   m <- setParams(m, alpha = 0.01, gamma = 0.01)
   check_case(m, Xta, cta, Xtf, ctf, Xtn, limType = "classic", caseNum = 4)
});


# test_that("case1: no prep, with outliers, no exclvars", {
#    Xta <- mda.exclrows(Xta, 32)
#    Xtn <- mda.exclrows(Xtn, 20)
#    Xtf <- mda.exclrows(Xtf, c(32, 60, 120))

#    m <- ddsimca(Xc, "Barolo", 10, scale = TRUE, exclrows = 32)
#    rta <- predict(m, Xta, cta)
#    rtf <- predict(m, Xtf, ctf)
#    rtn <- predict(m, Xtn)

#    writeCSV(m$calres, 'res_exp_train_noprep.csv',  ncomp = 2, src = src, name = 'train', sep = ';', dataFile = 'wines_baro')
#    writeCSV(ra,       'res_exp_alt_noprep.csv',    ncomp = 2, src = src, name = 'test', sep = ';', dataFile = 'wines_grigno')
#    writeCSV(rf,       'res_exp_full_noprep.csv',   ncomp = 2, src = src, name = 'test', sep = ';', dataFile = 'wines_full')
#    writeCSV(rn,       'res_exp_new_noprep.csv',    ncomp = 2, src = src, name = 'new', sep = ';', dataFile = 'wines_new')
# })


# #    # calres
# #    par(mfcol = c(2, 3))
# #    plotAcceptance(m$calres)
# #    plotAcceptance(m$calres, show.excluded = TRUE)
# #    plotAcceptance(m$calres, ncomp = 1, limType = "robust", log = TRUE)
# #    plotAcceptance(m$calres, ncomp = 1, limType = "robust", log = TRUE, show.excluded = TRUE)
# #    plotAcceptance(m$calres, ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE)
# #    plotAcceptance(m$calres, ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE, show.excluded = TRUE)

# #    # rta
# #    par(mfcol = c(2, 3))
# #    plotAcceptance(rta)
# #    plotAcceptance(rta, show.excluded = TRUE)
# #    plotAcceptance(rta, ncomp = 1, limType = "robust", log = TRUE)
# #    plotAcceptance(rta, ncomp = 1, limType = "robust", log = TRUE, show.excluded = TRUE)
# #    plotAcceptance(rta, ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE)
# #    plotAcceptance(rta, ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE, show.excluded = TRUE)

# #    # rtn
# #    par(mfcol = c(2, 3))
# #    plotAcceptance(rtn)
# #    plotAcceptance(rtn, show.excluded = TRUE)
# #    plotAcceptance(rtn, ncomp = 1, limType = "robust", log = TRUE)
# #    plotAcceptance(rtn, ncomp = 1, limType = "robust", log = TRUE, show.excluded = TRUE)
# #    plotAcceptance(rtn, ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE)
# #    plotAcceptance(rtn, ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE, show.excluded = TRUE)

# #    # rtf
# #    par(mfcol = c(2, 3))
# #    plotAcceptance(rtf)
# #    plotAcceptance(rtf, show.excluded = TRUE)
# #    plotAcceptance(rtf, ncomp = 1, limType = "robust", log = TRUE)
# #    plotAcceptance(rtf, ncomp = 1, limType = "robust", log = TRUE, show.excluded = TRUE)
# #    plotAcceptance(rtf, ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE)
# #    plotAcceptance(rtf, ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE, show.excluded = TRUE)

# #    par(mfcol = c(2, 3))
# #    plotAcceptance(rtf, show = "members")
# #    plotAcceptance(rtf, show = "members", show.excluded = TRUE)
# #    plotAcceptance(rtf, show = "members", ncomp = 1, limType = "robust", log = TRUE)
# #    plotAcceptance(rtf, show = "members", ncomp = 1, limType = "robust", log = TRUE, show.excluded = TRUE)
# #    plotAcceptance(rtf, show = "members", ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE)
# #    plotAcceptance(rtf, show = "members", ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE, show.excluded = TRUE)


# #    par(mfcol = c(2, 3))
# #    plotAcceptance(rtf, show = "strangers")
# #    plotAcceptance(rtf, show = "strangers", show.excluded = TRUE)
# #    plotAcceptance(rtf, show = "strangers", ncomp = 1, limType = "robust", log = TRUE)
# #    plotAcceptance(rtf, show = "strangers", ncomp = 1, limType = "robust", log = TRUE, show.excluded = TRUE)
# #    plotAcceptance(rtf, show = "strangers", ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE)
# #    plotAcceptance(rtf, show = "strangers", ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE, show.excluded = TRUE)

# #    par(mfcol = c(2, 3))
# #    plotAcceptance(rtf, show = "all")
# #    plotAcceptance(rtf, show = "all", show.excluded = TRUE)
# #    plotAcceptance(rtf, show = "all", ncomp = 1, limType = "robust", log = TRUE)
# #    plotAcceptance(rtf, show = "all", ncomp = 1, limType = "robust", log = TRUE, show.excluded = TRUE)
# #    plotAcceptance(rtf, show = "all", ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE)
# #    plotAcceptance(rtf, show = "all", ncomp = 4, limType = "robust", log = TRUE, show.labels = TRUE, show.excluded = TRUE)

# #    # as data frame
# #    print(as.data.frame(rta))
# #    print(as.data.frame(rtn))
# #    print(as.matrix(rta))
# #    print(as.matrix(rtf))
# #    print(as.matrix(rtn))

# #    writeCSV(rta, ncomp = 2, limType = "classic", fileName = "ddsimcares-alt.csv")
# #    writeCSV(rtf, ncomp = 2, limType = "classic", fileName = "ddsimcares-full.csv")
# #    writeCSV(rtn, ncomp = 2, limType = "robust", fileName = "ddsimcares-new.csv")


# #    par(mfcol = c(3, 4))
# #    plotDistances(m$calres)
# #    plotDistances(m$calres, ncomp = 4, distance = "q")
# #    plotDistances(m$calres, limType = "robust", ncomp = 2, distance = "f", log = TRUE)
# #    plotDistances(rta)
# #    plotDistances(rta, ncomp = 4, distance = "q")
# #    plotDistances(rta, limType = "robust", ncomp = 2, distance = "f", log = TRUE)
# #    plotDistances(rtf)
# #    plotDistances(rtf, ncomp = 4, distance = "q")
# #    plotDistances(rtf, limType = "robust", ncomp = 2, distance = "f", log = TRUE)
# #    plotDistances(rtn)
# #    plotDistances(rtn, ncomp = 4, distance = "q")
# #    plotDistances(rtn, limType = "robust", ncomp = 2, distance = "f", log = TRUE)


# #    m <- ddsimca(Xc, "Barolo", 10, scale = TRUE, exclrows = 32)
# #    rta <- predict(m, Xta, cta)
# #    rtf <- predict(m, Xtf, ctf)
# #    rtn <- predict(m, Xtn)

# #    writeJSON(m, "ddsimca-model-r.json")
# # })