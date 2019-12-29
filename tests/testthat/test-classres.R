#####################################
# Tests for classres class methods  #
#####################################

# mock data

## reference values
c.ref <- as.factor(c(rep("C2", 10), rep("C1", 10), rep("C3", 10)))

## predicted values
nrows <- 30
ncomp <- 3
nclasses <- 2
n <- nrows * ncomp * nclasses
classnames <- paste0("C", seq_len(nclasses))

set.seed(42)
c.pred <- array((rnorm(n) > 0) * 2 - 1, dim = c(nrows, ncomp, nclasses))
dimnames(c.pred)[[2]] <- paste0("Comp ", seq_len(ncomp))
dimnames(c.pred)[[3]] <- classnames

## excluded values
excluded_rows <- c(1, 7, 15, 25)


######################################################
# Block 1: tests for getClassificationPerformance()  #
######################################################

## testing function
tf <- function(res, tp, tn, fp, fn, classnames, ncomp) {
   sn <- tp / (tp + fn)
   sn <- rbind(sn, colSums(tp) / colSums(tp + fn))

   sp <- tn / (tn + fp)
   sp <- rbind(sp, colSums(tn) / colSums(tn + fp))

   ms <- (fp + fn) / (tn + fp + tp + fn)
   ms <- rbind(ms, colSums(fn + fp) / colSums(tn + tp + fn + fp))


   # 1. check main outcomes

   ## check dimension
   expect_equal(dim(res$tp), c(length(classnames), ncomp))
   expect_equal(dim(res$tn), c(length(classnames), ncomp))
   expect_equal(dim(res$fp), c(length(classnames), ncomp))
   expect_equal(dim(res$fn), c(length(classnames), ncomp))

   ## check names
   expect_equal(rownames(res$tp), classnames)
   expect_equal(rownames(res$tn), classnames)
   expect_equal(rownames(res$fp), classnames)
   expect_equal(rownames(res$fn), classnames)

   ## check values
   expect_equivalent(res$tp, tp)
   expect_equivalent(res$tn, tn)
   expect_equivalent(res$fp, fp)
   expect_equivalent(res$fn, fn)

   # 2. check statistics

   ## check dimension
   expect_equal(dim(res$sensitivity), c(length(classnames) + 1, ncomp))
   expect_equal(dim(res$specificity), c(length(classnames) + 1, ncomp))
   expect_equal(dim(res$misclassified), c(length(classnames) + 1, ncomp))

   ## check names
   expect_equal(rownames(res$sensitivity), c(classnames, "Total"))
   expect_equal(rownames(res$specificity), c(classnames, "Total"))
   expect_equal(rownames(res$misclassified), c(classnames, "Total"))

   ## check values
   expect_equivalent(res$sensitivity, sn)
   expect_equivalent(res$specificity, sp)
   expect_equivalent(res$misclassified, ms)

   ## 3.
}

# 1. only c.pred, no c.ref
context("classres: no reference data")
expect_error(classres.getPerformance(NULL, c.pred))

# 2. c.ref and c.pred have different classes

context("classres: predicted and reference have different classes")

tp <- rbind(c(4, 5, 4), c(5, 2, 2))
fn <- rbind(c(6, 5, 6), c(5, 8, 8))
tn <- rbind(c(11, 7, 6), c(11, 13, 9))
fp <- rbind(c(9, 13, 14), c(9, 7, 11))

res <- classres.getPerformance(c.ref, c.pred)
tf(res, tp, tn, fp, fn, classnames, ncomp)

# 3. c.ref and c.pred have same classes
context("classres: predicted and reference have same classes")

tp <- rbind(c(4, 5, 4), c(5, 2, 2))
fn <- rbind(c(6, 5, 6), c(5, 8, 8))
tn <- rbind(c(4, 5, 2), c(4, 5, 4))
fp <- rbind(c(6, 5, 8), c(6, 5, 6))

res <- classres.getPerformance(c.ref[1:20], c.pred[1:20, , ])
tf(res, tp, tn, fp, fn, classnames, ncomp)

# 4. same as #2 but with excluded rows

context("classres: different classes and exclrows")

## exclude some values
c.predexcl <- c.pred
attr(c.predexcl, "exclrows") <- excluded_rows

tp <- rbind(c(4, 5, 4), c(4, 2, 2))
fn <- rbind(c(5, 4, 5), c(4, 6, 6))
tn <- rbind(c(11, 6, 4), c(9, 12, 9))
fp <- rbind(c(6, 11, 13), c(9, 6, 9))

res <- classres.getPerformance(c.ref, c.predexcl)
tf(res, tp, tn, fp, fn, classnames, ncomp)

# 5. same as #4 but with only one component

context("classres: different classes and exclrows (one component)")

## subset and exclude some values
c.predexcl <- c.pred[, 1, , drop = FALSE]
attr(c.predexcl, "exclrows") <- excluded_rows

tp <- rbind(c(4), c(4))
fn <- rbind(c(5), c(4))
tn <- rbind(c(11), c(9))
fp <- rbind(c(6), c(9))

res <- classres.getPerformance(c.ref, c.predexcl)
tf(res, tp, tn, fp, fn, classnames, ncomp = 1)

# 6. same as #4 but with only one class

context("classres: different classes and exclrows (one class)")

## subset and exclude some values
c.predexcl <- c.pred[, , 1, drop = FALSE]
attr(c.predexcl, "exclrows") <- excluded_rows

tp <- rbind(c(4, 5, 4))
fn <- rbind(c(5, 4, 5))
tn <- rbind(c(11, 6, 4))
fp <- rbind(c(6, 11, 13))

res <- classres.getPerformance(c.ref, c.predexcl)
tf(res, tp, tn, fp, fn, classnames[[1]], ncomp)

# 7. same as #4 but with only one class and one component

context("classres: different classes and exclrows (one class + one component)")

## subset and exclude some values
c.predexcl <- c.pred[, 1, 1, drop = FALSE]
attr(c.predexcl, "exclrows") <- excluded_rows

tp <- rbind(c(4))
fn <- rbind(c(5))
tn <- rbind(c(11))
fp <- rbind(c(6))

res <- classres.getPerformance(c.ref, c.predexcl)
tf(res, tp, tn, fp, fn, classnames[[1]], ncomp = 1)

# 8. test classres() constructor

context("classres: constructor")

## exclude some values
c.predexcl <- c.pred
attr(c.predexcl, "exclrows") <- excluded_rows

res <- classres.getPerformance(c.ref, c.predexcl)
classres <- classres(c.predexcl, c.ref)
tf(classres, res$tp, res$tn, res$fp, res$fn, classnames, ncomp)

## make sure wrong number of dimensions in c.pref lead to error
expect_error(classress(c.predexcl[, 1], c.ref))

## if no reference - no statistics
classres <- classres(c.predexcl, NULL)
expect_null(classres$tp)
expect_null(classres$tn)
expect_null(classres$fp)
expect_null(classres$fn)
expect_null(classres$sensitivity)
expect_null(classres$specificity)
expect_null(classres$misclassified)
expect_null(classres$c.ref)
expect_equal(classres$c.pred, c.predexcl)

############################################
# Block 2: tests for getConfusionMatrix()  #
############################################

context("classres: getConfusionMatrix()")

## for all classes

c.predexcl <- c.pred
attr(c.predexcl, "exclrows") <- excluded_rows
res <- classres(c.predexcl, c.ref)

## expected values
cm.ref1 <- cbind(c(4, 4, 2), c(6, 4, 3), c(3, 3, 4))
cm.ref2 <- cbind(c(5, 4, 7), c(4, 2, 2), c(3, 4, 1))
cm.ref3 <- cbind(c(4, 7, 6), c(5, 2, 4), c(3, 1, 1))

expect_equivalent(cm.ref1, getConfusionMatrix(res, 1))
expect_equivalent(cm.ref2, getConfusionMatrix(res, 2))
expect_equivalent(cm.ref3, getConfusionMatrix(res, 3))

## for one class

c.predexcl <- c.pred[, , 1, drop = FALSE]
attr(c.predexcl, "exclrows") <- excluded_rows
res <- classres(c.predexcl, c.ref)

## expected values
cm.ref1 <- cbind(c(4, 4, 2), c(5, 4, 7))
cm.ref2 <- cbind(c(5, 4, 7), c(4, 4, 2))
cm.ref3 <- cbind(c(4, 7, 6), c(5, 1, 3))

expect_equivalent(cm.ref1, getConfusionMatrix(res, 1))
expect_equivalent(cm.ref2, getConfusionMatrix(res, 2))
expect_equivalent(cm.ref3, getConfusionMatrix(res, 3))
