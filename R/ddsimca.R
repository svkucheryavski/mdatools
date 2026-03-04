#' Data Driven SIMCA
#'
#' @description
#' \code{ddsimca} is used to develop DD-SIMCA (Data Driven SIMCA) model for
#' one-class classification.
#'
#' @param x
#' a numerical matrix with data values.
#' @param classname
#' short text (up to 20 symbols) with class name.
#' @param ncomp
#' maximum number of components to calculate.
#' @param center
#' logical, do mean centering of data or not.
#' @param scale
#' logical, do standardization of data or not.
#' @param pcv
#' Procrustes cross-validation settings (see details).
#' @param alpha
#' significance level for making the predictions (can be also adjusted when model is applied to data).
#' @param gamma
#' significance level for detection of outliers (can be also adjusted when model is applied to data).
#' @param exclrows
#' rows to be excluded from calculations (numbers, names or vector with logical values)
#' @param exclcols
#' columns to be excluded from calculations (numbers, names or vector with logical values)
#' @param prep
#' optional list with preprocessing methods created using `\code{\link{prep}}` function.
#' @param do.round
#' logical, round or not DoF for distances.
#' @param ...
#' any other parameters suitable for \code{\link{pca}} method.
#'
#' @details
#' DD-SIMCA is based on PCA model with additional functionality, so \code{ddsimca} class inherits most
#' of the functionality of \code{\link{pca}} class.
#'
#' In order to make a decision, DDSIMCA uses score and orthogonal distances to PCA model. It combines
#' the two distances to joint full distance and uses chi-distribution for finding a critical value
#' which is employed as decision rule. More details about DD-SIMCA can be found in [1] (open access).
#'
#' Procrustes cross-validation (PCV) is used to generate a validation set in order to find optimal model
#' complexity (number of components). The PCV settings are similar to the ones used for conventional
#' cross-validation. The best way is to set `pcv` value to a list, for example:  \code{pcv = list('ven', nseg)}
#' for systematic splits or \code{pcv = list('rand', nseg)} for random splits. In case if full cross-validation
#' must be employed, use \code{pcv = list('loo')}.
#'
#' @return
#' Returns an object of \code{ddsimca} class with following fields:
#' \item{classname }{a short text with class name.}
#' \item{calres }{an object of class \code{\link{simcares}} with classification results for a
#' calibration data.}
#' \item{pvres }{an object of class \code{\link{simcares}} with classification results for a Procrustes
#' validation set.}
#'
#' Fields, inherited from \code{\link{pca}} class:
#' \item{ncomp }{number of components included to the model.}
#' \item{ncomp.selected }{selected (optimal) number of components.}
#' \item{loadings }{matrix with loading values (nvar x ncomp).}
#' \item{eigenvals }{vector with eigenvalues for all existent components.}
#' \item{expvar }{vector with explained variance for each component (in percent).}
#' \item{cumexpvar }{vector with cumulative explained variance for each component (in percent).}
#' \item{info }{information about the model, provided by user when build the model.}
#'
#' @references
#' 1. Kucheryavskiy S, Rodionova O, Pomerantsev A. A comprehensive tutorial on Data-Driven SIMCA:
#' Theory and implementation in web. Journal of Chemometrics. 2024; 38(7):e3556. doi:10.1002/cem.3556
#'
#' 2. S. Kucheryavskiy, O. Rodionova, A. Pomerantsev, Procrustes cross-validation of multivariate
#' regression models. Analytica Chimica Acta. 2023; 1255:341096. doi:10.1016/j.aca.2023.341096.
#'
#' @seealso
#' Methods for \code{ddsimca} objects:
#' \tabular{ll}{
#'  \code{print.ddsimca} \tab shows information about the object.\cr
#'  \code{summary.ddsimca} \tab shows summary statistics for the model.\cr
#'  \code{plot.ddsimca} \tab makes an overview of DD-SIMCA model with four plots.\cr
#'  \code{\link{predict.ddsimca}} \tab applies DD-SIMCA model to a new data.\cr
#' }
#'
#' Methods, inherited from \code{classmodel} class:
#' \tabular{ll}{
#'  \code{\link{plotPredictions.classmodel}} \tab shows plot with predicted values.\cr
#'  \code{\link{plotSensitivity.classmodel}} \tab shows sensitivity plot.\cr
#'  \code{\link{plotSpecificity.classmodel}} \tab shows specificity plot.\cr
#'  \code{\link{plotMisclassified.classmodel}} \tab shows misclassified ratio plot.\cr
#' }
#'
#' Methods, inherited from \code{\link{pca}} class:
#' \tabular{ll}{
#'  \code{\link{selectCompNum.pca}} \tab set number of optimal components in the model\cr
#'  \code{\link{plotScores.pca}} \tab shows scores plot.\cr
#'  \code{\link{plotLoadings.pca}} \tab shows loadings plot.\cr
#'  \code{\link{plotVariance.pca}} \tab shows explained variance plot.\cr
#'  \code{\link{plotCumVariance.pca}} \tab shows cumulative explained variance plot.\cr
#' }
#'
#' @author
#' Sergey Kucheryavskiy (svkucheryavski@@gmail.com)
#'
#' @examples
#' ## make a SIMCA model for Iris setosa class with full cross-validation
#' library(mdatools)
#'
#' data = iris[, 1:4]
#' class = iris[, 5]
#'
#' # take first 20 objects of setosa as calibration set
#' se = data[1:20, ]
#'
#' # make SIMCA model and apply to test set
#' model = ddsimca(se, "setosa", pcv = list("ven", 10))
#' model = selectCompNum(model, 1)
#'
#' # show information, summary and plot overview
#' print(model)
#' summary(model)
#' plot(model)
#'
#'
#' @importFrom pcv pcvpca
#'
#' @export
ddsimca <- function(x, classname, ncomp = min(nrow(x) - 1, ncol(x) - 1, 20), center = TRUE,
   scale = FALSE, pcv = list("ven", 10), alpha = 0.05, gamma = 0.01,
   exclrows = NULL, exclcols = NULL, prep = NULL, do.round = TRUE, ...) {

   if (!is.character(classname)) {
      stop("Argument 'classname' must be a text.", call. = FALSE)
   }

   if (nchar(classname) > 20) {
      stop("Argument 'classname' must have up to 20 symbols.", call. = FALSE)
   }

   # check calibration data and process excluded rows and columns
   x <- prepCalData(x, exclrows = exclrows, exclcols = exclcols, min.nrows = 2, min.ncols = 2)
   x <- mda.purgeCols(x)

   # correct number of components
   ncomp <- min(nrow(x) - 1, ncol(x) - 1 - length(attr(x, "exclcols")), ncomp)

   # calibrate model
   model <- pca(x, ncomp = ncomp, center = center, scale = scale, prep = prep, do.round = do.round, ...)
   model$nrows <- nrow(x) - length(attr(x, "exclrows"))
   model$nclasses <- 1
   model$classname <- classname
   model$alpha <- alpha
   model$gamma <- gamma
   model$limType <- "moments"
   model$call <- match.call()
   class(model) <- c("ddsimca", "pca")

   # apply model to calibration set
   c.ref <- rep(classname, nrow(x))
   model$res[["cal"]] <- predict(model, x, c.ref)
   model$calres <- model$res[["cal"]]

   # do Procrustes cross-validation if needed
   if (!is.null(pcv)) {

      # remove excluded rows (capture processed indices before purging clears them)
      exclrows_idx <- attr(x, "exclrows")
      x <- mda.purgeRows(x)
      c.ref.pv <- if (length(exclrows_idx) > 0) c.ref[-exclrows_idx] else c.ref

      # if preprocessing is available apply it
      if (!is.null(model$prep)) x <- prep.apply(model$prep, x)

      # generate xpv set
      xpv <- pcv::pcvpca(x, ncomp = ncomp, cv = pcv, center = center, scale = scale)

      # make a copy of PCA model and remove preprocessing - already applied
      # then apply the model
      m <- model
      m$prep <- NULL
      model$res[["pv"]] <- predict(m, xpv, c.ref.pv)
      model$pvres <- model$res[["pv"]]
   }

   return(model)
}


#' DD-SIMCA predictions
#'
#' @description
#' Applies DD-SIMCA model to a new data set
#'
#' @param object
#' a DD-SIMCA model (object of class \code{ddsimca})
#' @param x
#' a matrix with x values (predictors)
#' @param c.ref
#' a vector with reference class names (same as class names for models)
#' @param alpha
#' significance level for making the predictions.
#' @param gamma
#' significance level for detection of outliers.
#' @param ...
#' other optional parameters.
#'
#' @return
#' DD-SIMCA results (an object of class \code{ddsimcares})
#'
#' @details
#' See examples in help for \code{\link{ddsimca}} function.
#'
#' @export
predict.ddsimca <- function(object, x, c.ref = NULL, alpha = object$alpha, gamma = object$gamma, ...) {

   if (alpha > 0.5 || alpha < 0.00001) {
      stop("Wrong value for parameter 'alpha'.")
   }

   if (gamma > 0.5 || gamma < 0.00001) {
      stop("Wrong value for parameter 'gamma'.")
   }

   # get PCA results
   pcares <- predict.pca(object, x)
   nobj <- nrow(x)
   nobj.cal <- object$nrows
   classname <- object$classname

   indIncluded <- rep(TRUE, nobj)
   if (length(attr(pcares$scores, "exclrows"))> 0) {
      indIncluded[attr(pcares$scores, "exclrows")] <- FALSE
   }
   indExcluded <- !indIncluded
   nExcluded <- sum(indExcluded)

   # identify indices for strangers, members and object whose class is unknown
   if (is.null(c.ref)) {
      indMembers <- NULL
      indStrangers <- NULL
      indUnknowns <- indIncluded
      nMembers <- 0
      nStrangers <- 0
      nUnknowns <- sum(indUnknowns)
   } else {
      indMembers <- (c.ref == object$classname) & indIncluded
      indStrangers <- (!indMembers) & indIncluded
      indUnknowns <- NULL
      nMembers <- sum(indMembers)
      nStrangers <- sum(indStrangers)
      nUnknowns <- 0
   }

   indices <- list(members = indMembers, strangers = indStrangers, unknown = indUnknowns, excluded = indExcluded)
   numbers <- list(members = nMembers,   strangers = nStrangers,   unknown = nUnknowns, excluded = nExcluded)

   # compute full distances and probabilities
   lp <- object$limParams
   outcomes <- list(
      "moments" = classify(pcares, indices, numbers, lp$Q[["moments"]], lp$T2[["moments"]], nobj.cal, alpha, gamma, classname, c.ref),
      "robust" = classify(pcares, indices, numbers, lp$Q[["robust"]],  lp$T2[["robust"]],  nobj.cal, alpha, gamma, classname, c.ref)
   )

   return(ddsimcares(pcares, outcomes, indices, numbers, alpha = alpha, classname = classname, c.ref = c.ref))
}

#' Set default parameters for the DD-SIMCA model.
#'
#' @param obj
#' a DD-SIMCA model (object of class \code{ddsimca})
#' @param alpha
#' significance level for making the predictions.
#' @param gamma
#' significance level for detection of outliers.
#' @param ...
#' other arguments
#'
#' @return DD-SIMCA model with the redefined parameters.
#'
#' @export
setParams.ddsimca <- function(obj, alpha = obj$alpha, gamma = obj$gamma, ...) {
   obj$alpha <- alpha
   obj$gamma <- gamma

   nobj <- obj$nrows
   lp <- obj$limParams
   classname <- obj$classname
   if (!is.null(obj[["res"]]) && !is.null(obj$res[["cal"]])) {

      r <- obj[["res"]][["cal"]]
      indices <- r$simca$indices
      numbers <- r$simca$numbers
      c.ref <- r$simca$c.ref

      outcomes <- list(
         "moments" = classify(r, indices, numbers, lp$Q[["moments"]], lp$T2[["moments"]], nobj, alpha, gamma, classname, c.ref),
         "robust" = classify(r, indices, numbers, lp$Q[["robust"]],  lp$T2[["robust"]],  nobj, alpha, gamma, classname, c.ref)
      )

      obj$res$cal <- ddsimcares(r, outcomes, indices, numbers, alpha = alpha, classname = classname, c.ref = c.ref)
      obj$calres <- obj$res$cal
   }

   if (!is.null(obj[["res"]]) && !is.null(obj$res[["pv"]])) {

      r <- obj[["res"]][["pv"]]
      indices <- r$simca$indices
      numbers <- r$simca$numbers
      c.ref <- r$simca$c.ref

      outcomes <- list(
         "moments" = classify(r, indices, numbers, lp$Q[["moments"]], lp$T2[["moments"]], nobj, alpha, gamma, classname, c.ref),
         "robust" = classify(r, indices, numbers, lp$Q[["robust"]],  lp$T2[["robust"]],  nobj, alpha, gamma, classname, c.ref)
      )

      obj$res$pv <- ddsimcares(r, outcomes, indices, numbers, alpha = alpha, classname = classname, c.ref = c.ref)
      obj$pvres <- obj$res$pv
   }

   return(obj)
}


#' Model overview plot for DD-SIMCA
#'
#' @description
#' Shows a set of plots for DD-SIMCA model.
#'
#' @param x
#' a DD-SIMCA model (object of class \code{ddsimca})
#' @param comp
#' which components to show on scores and loadings plot
#' @param ncomp
#' how many components to use for residuals plot
#' @param ...
#' other arguments
#'
#' @details
#' See examples in help for \code{\link{ddsimca}} function.
#'
#' @export
plot.ddsimca <- function(x, comp = c(1, 2), ncomp = x$ncomp.selected, ...) {
   op <- par(mfrow = c(2, 2))
   on.exit(par(op))
   plotScores(x, comp, ...)
   plotLoadings(x, comp = comp, ...)
   plotAcceptance(x, ncomp = ncomp, ...)
   plotSensitivity(x, ...)
}

#' Summary method for DD-SIMCA model object
#'
#' @description
#' Shows performance statistics for the model.
#'
#' @param object
#' a DD-SIMCA model (object of class \code{ddsimca})
#' @param ncomp
#' number of components to show summary for
#' @param res
#' list of result objects to show summary for
#' @param ...
#' other arguments
#'
#' @export
summary.ddsimca <- function(object, ncomp = object$ncomp.selected, res = object$res, ...) {

   fprintf("\nDD-SIMCA model for class '%s'\n\n", object$classname)

   if (!is.null(object$info) && nchar(object$info)) {
      fprintf("Info: %s\n", object$info)
   }

   if (!is.null(object$rand)) {
      fprintf("Parameters for randomized algorithm: q = %d, p = %d\n",
         object$rand[1], object$rand[2])
   }

   if (length(object$exclrows) > 0) {
      fprintf("Excluded rows: %d\n", length(object$exclrows))
   }

   if (length(object$exclcols) > 0) {
      fprintf("Excluded columns: %d\n", length(object$exclcols))
   }

   fprintf("Number of components: %d\n", object$ncomp)
   fprintf("Number of selected components: %d\n", ncomp)
   fprintf("Type of limits: %s\n", if (object$limType == "robust" || object$limType == "ddrobust") "robust" else "classic")
   fprintf("Alpha: %s\n", object$alpha)
   fprintf("Gamma: %s\n", object$gamma)

   if (!is.null(object$prep)) {
      fprintf("\nPreprocessing methods:\n")
      print(object$prep)
   }

   cat("\n")

   if (is.null(object[["res"]]) || is.null(object$res[["cal"]])) {
      message("Calibration results not found (most probably this model is loaded from web-application).")
      return(invisible(object))
   }


   sum_data <- do.call(rbind, lapply(res, function(x) as.matrix(x)[ncomp, , drop = FALSE]))
   rownames(sum_data) <- capitalize(names(res))
   print(sum_data, 4)
   cat("\n")

   return(invisible(object))
}


#' Print method for DD-SIMCA model object
#'
#' @description
#' Prints information about the object structure
#'
#' @param x
#' a DD-SIMCA model (object of class \code{ddsimca})
#' @param ...
#' other arguments
#'
#' @export
print.ddsimca <- function(x, ...) {
   cat("\nDD-SIMCA model (class 'ddsimca')\n")

   cat("\nCall:\n")
   print(x$call)

   cat("\nFields inherited from 'pca' class:\n")
   cat(" $info - information about the model\n")
   cat(" $classname - name of the class\n")
   cat(" $ncomp - number of calculated components\n")
   cat(" $ncomp.selected - number of selected components\n")
   cat(" $loadings - matrix with loadings\n")
   cat(" $eigenvals - eigenvalues for components\n")
   cat(" $center - values for centering data\n")
   cat(" $scale - values for scaling data\n")
   cat(" $limParams - parameters for distribution of distances\n")

   cat("\nFields related to 'ddsimca' class:\n")
   cat(" $alpha - significance level for critical limits\n")
   cat(" $gamma - significance level for outlier limits\n")
   cat(" $classname - name of the target class\n")
   cat(" $nrows - number of objects in the training set\n")

   if (!is.null(x$res))
      cat(" $res - list with model results (calibration, PV-set)\n")
   if (!is.null(x$prep))
      cat(" $prep - preprocessing model\n")

   return(invisible(x))
}


################################
#  Static and local methods    #
################################

#' Creates classification outcomes for given PCA result objects and distance parameters.
#'
#' @param pcares
#' object of class \code{\link{pcares}}.
#' @param indices
#' list with indices of members, strangers, and unknown samples as well as information about excluded objects.
#' @param numbers
#' list with numbers of members, strangers, unknown and excluded samples.
#' @param qp
#' distribution parameters for q-distances.
#' @param hp
#' distribution parameters for h-distances.
#' @param nobj.cal
#' number of objects in calibration set the model was built on (needed for outliers detection).
#' @param alpha
#' significance level for decision boundary.
#' @param gamma
#' significance level for outlier detection boundary,
#' @param classname
#' name of the target class.
#' @param c.ref
#' vector with names of the reference classes.
#'
#' @return a list with classification outcomes for each number of components.
classify <- function(pcares, indices, numbers, qp, hp, nobj.cal, alpha = 0.05, gamma = 0.01, classname, c.ref = NULL) {

   # get distances and compute matrix with full distance
   H <- pcares$T2
   Q <- pcares$Q


   # get distribution parameters
   q0 <- qp$u0
   Nq <- qp$Nu
   h0 <- hp$u0
   Nh <- hp$Nu
   Nf <- Nh + Nq

   # compute critical limits
   fCritE <- qchisq(1 - alpha, Nf)
   fCritO <- qchisq((1 - gamma)^(1/nobj.cal), Nf)

   ncomp <- ncol(H)
   nobj <- nrow(H)

   # vectors for main FoMs
   sns <- rep(0, ncomp)
   spc <- rep(0, ncomp)
   sel <- rep(0, ncomp)
   eff <- rep(0, ncomp)
   acc <- rep(0, ncomp)

   # vectors for number of objects accepted/rejected
   nin <- rep(0, ncomp)
   nout <- rep(0, ncomp)

   # vectors with number of negatives/positives
   TN <- rep(0, ncomp)
   TP <- rep(0, ncomp)
   FN <- rep(0, ncomp)
   FP <- rep(0, ncomp)

   # list to save the classification outcomes for each component
   classres <- list()

   for (a in seq_len(ncomp)) {

      ha <- H[, a]
      qa <- Q[, a]
      fa <- ha / h0[a] * Nh[a] + qa / q0[a] * Nq[a]

      Nfa <- Nf[a]
      Nha <- Nh[a]
      Nqa <- Nq[a]

      fcea <- fCritE[a]
      fcoa <- fCritO[a]
      names(fa) <- NULL

      res <- list(
         f = fa,
         h = ha / h0[a],
         q = qa / q0[a],
         Nh = Nha,
         Nq = Nqa,
         Nf = Nfa,
         fce = fcea,
         fco = fcoa,
         roles = rep("", nobj),
         decisions = fa < fcea,
         TP = 0,
         FN = 0,
         TN = 0,
         FP = 0,
         beta = 0,
         s = 0,
         f0 = 0,
         hz = 0,
         Mz = 0,
         Sz = 0,
         k = 0,
         m = 0
      )

      res$roles[indices$excluded] <- "excluded"

      res <- processMembers(res, indices$members)
      res <- processStrangers(res, indices$strangers)

      if (!is.null(indices$unknown)) {
         res <- processStrangers(res, indices$unknown)
         res$TN <- 0
         res$FP <- 0
      }

      TP[a] <- res$TP
      TN[a] <- res$TN
      FP[a] <- res$FP
      FN[a] <- res$FN

      nin[a] <- sum(res$decisions)
      nout[a] <- sum(!res$decisions)

      sns[a] <- if (res$TP > 0) res$TP / (res$TP + res$FN) else 0
      spc[a] <- if (res$TN > 0) res$TN / (res$TN + res$FP) else 0
      sel[a] <- if (res$beta > 0) 1 - res$beta else 0
      eff[a] <- sqrt(sns[a] * spc[a])
      acc[a] <- if(res$TP > 0 && res$TN > 0) (res$TP + res$TN) / (res$TP + res$TN + res$FP + res$FN) else 0

      classres[[a]] <- res
   }

   return(list(values = classres, sns = sns, spc = spc, sel = sel, eff = eff, acc = acc, nin = nin, nout = nout, TP = TP, TN = TN, FP = FP, FN = FN))
}

#' Computes classification outcomes for target class members.
#'
#' @param res
#' a list with classification outcomes created by method \code{\link{classify}} (part of it will be filled by this method).
#' @param indMembers
#' a vector with logical values pointing on data items corresponding to class members.
#'
#' @return
#' the \code{res} list where roles vector is filled for class members, plus values for TP and FN.
processMembers <- function(res, indMembers) {
   if (is.null(indMembers) || sum(indMembers) < 1) return(res)

   regular <- indMembers & res$f < res$fce
   outlier <- indMembers & res$f > res$fco
   extreme <- indMembers & !regular & !outlier

   if (any(regular)) res$roles[regular] <- "regular"
   if (any(extreme)) res$roles[extreme] <- "extreme"
   if (any(outlier)) res$roles[outlier] <- "outlier"

   res$TP <- sum(regular)
   res$FN <- sum(outlier) + sum(extreme)

   return(res)
}

#' Computes classification outcomes for members of non-target classes.
#'
#' @param res
#' a list with classification outcomes created by method \code{\link{classify}} (part of it will be filled by this method).
#' @param indStrangers
#' a vector with logical values pointing on data items corresponding to members of non-target classes.
#'
#' @return
#' the \code{res} list where roles vector is filled for class members, values for TN and FP, as
#' well as statistics for computing Type II error and related things (beta, s, f0, hz, Mz, Sz, k, m).
processStrangers <- function(res, indStrangers) {

   if (is.null(indStrangers) || sum(indStrangers) < 1) return(res)

   res$roles[indStrangers] <- "alien"
   res$TN <- sum(indStrangers & !res$decisions)
   res$FP <- sum(indStrangers & res$decisions)

   # not enough strangers to fit non-central chi-square — skip beta computation
   if (sum(indStrangers) < 2) return(res)

   # Step 1. Sort all distances
   # make a copy of strangers indices as we are going to change it
   indv <- which(indStrangers)

   # sort the indices so the largest will be at the end
   indv <- indv[order(res$f[indv])]
   I <- length(indv)

   # sort the distance values as well and save into temporary variable
   fp <- res$f[indv]

   # Step 2. Try to fit the non-central chi-square
   Disc <- -1
   n <- 0
   m <- 0
   d <- 0
   M1 <- 0
   k <- res$Nf

   while (Disc < 0) {
      if (n > 0) {
         # sample with largest f does not fit, so we change its role to "external"
         # and amend the number of aliens and externals
         res$roles[indv[I]] <- "external"
         # then we remove this sample from the temporary vector and
         # assess the next biggest
         indv <- indv[-I]
         fp <- fp[-I]
         I <- I - 1
      }

      # not enough strangers left to estimate distribution parameters
      if (I < 2) break

      # compute parameters for moments squared equation and
      # its discriminant
      m <- mean(fp)
      d <- var(fp)
      M1 <- d / (m * m)
      Disc <- 4 - 2 * k * M1
      n <- n + 1
   }

   # if loop exhausted strangers or all distances are identical (M1 = 0), skip beta
   if (I < 2 || M1 < .Machine$double.eps) return(res)

   # Step 3. Calculate  x  by Eq. (18)
   x <- (2 + sqrt(Disc)) / M1

   # Calculate  s and  f'0 by Eq. (19)
   s <- (x - k)
   f0 <- m / x

   # Calculate: z, hz, r, p, Mz, Sz by Eq. (20)
   z <- res$fce / f0
   hz <- 1 - 2 * (k + s) * (k + 3 * s) / (3 * (k + 2 * s)^2)
   r <- (hz - 1) * (1 - 3 * hz)
   p <- (k + 2 * s) / (k + s)^2
   Mz <- 1 + hz * p * (hz - 1 - 0.5 * (2 - hz) * r * p)
   Sz <- hz * sqrt(2 * p) * (1 + 0.5 * r * p)

   # Step 4.	If α is given then calculate β by Eq. (21)
   beta <- pnorm((((z / (k + s))^hz) - Mz) / Sz)

   res$beta <- beta
   res$s <- s
   res$f0 <- f0
   res$hz <- hz
   res$Mz <- Mz
   res$Sz <- Sz
   res$k <- k
   res$m <- m

   return(res)
}


################################
# JSON methods                 #
################################

#' Converts JSON string created in mda.tools/ddsimca app to \code{ddsimca} object
#'
#' @param str
#' stringified JSON (from model file)
#'
#' @return
#' object of \code{\link{ddsimca}} class
#'
ddsimca.fromjson <- function(str) {
   m <- pca.fromjson(str)
   simca.str <- extractBlock(str, "simca")

   m$nrows <- extractValue(str, "nrows")
   m$nclasses <- 1
   m$classname <- extractValue(str, "className")
   m$alpha <- extractValue(str, "alpha") / 100
   m$gamma <- extractValue(str, "gamma") / 100
   m$ncomp.selected <- extractValue(simca.str, "ncomp")
   m$limType <- extractValue(simca.str, "limType")
   class(m) <- c("ddsimca", "pca")
   return(m)
}

#' Reads DD-SIMCA model from JSON file made in web-application (mda.tools/ddsimca).
#'
#' @param fileName
#' file name (or full path) to JSON file.
#'
#' @return list with DD-SIMCA model similar to what \code{ddsimca()} creates.
#'
#' @export
ddsimca.readJSON <- function(fileName) {
   fileConn <- file(fileName)
   str <- readLines(fileConn, warn = FALSE)
   close(fileConn)
   return(ddsimca.fromjson(str))
}

#' Converts object with DD-SIMCA model to JSON string compatible with web-application.
#'
#' @param obj
#' Object with DD-SIMCA model (from \code{\link{ddsimca}}).
#' @param limType
#' type of critical limits to use ('classic' or 'robust').
#' @param alpha
#' significance level for decision boundary.
#' @param gamma
#' significance level for outlier detection boundary.
#' @param ncomp
#' number of selected components.
#' @param ...
#' other arguments
#'
#' @return stringified JSON
#'
#' @export
asjson.ddsimca <- function(obj, limType = "classic", alpha = obj$alpha, gamma = obj$gamma, ncomp = obj$ncomp.selected, ...) {

   limType <- processLimType(limType)
   hash <- paste0("'", genhash(), "'")

   Nf <- rep(0, obj$ncomp)
   eCrit <- rep(0, obj$ncomp)
   oCrit <- rep(0, obj$ncomp)

   for (i in seq_len(obj$ncomp)) {
      v <- obj$calres$simca$outcomes[[limType]]$values[[i]]
      Nf[i] <- v$Nf
      eCrit[i] <- v$fce
      oCrit[i] <- v$fco
   }

   if (limType == "moments") limType <- "classic"
   m <- paste0(
      "{'pca':", asjson.pca(obj), ", 'simca': {",
         "'ncomp':", ncomp, ",",
         "'alpha':", alpha * 100, ",",
         "'gamma':", gamma * 100, ",",
         "'limType': '", limType, "',",
         "'Nf':{'__type':'Float64Array','data':[",paste0(Nf, collapse = ","), "]},",
         "'eCrit':{'__type':'Float64Array','data':[",paste0(eCrit, collapse = ","), "]},",
         "'oCrit':{'__type':'Float64Array','data':[",paste0(oCrit, collapse = ","), "]},",
         "'className':'", obj$classname, "',",
         "'hash':", hash,
      "}}"
   )

   m <- gsub("\'", "\"", m)
}


#' Saves DD-SIMCA model as JSON file compatible with web-application (https://mda.tools/ddsimca).
#'
#' @description
#' You can load created JSON file to web-app and use it for prediction.
#'
#' @param obj
#' Object with DD-SIMCA model (from \code{\link{ddsimca}}).
#' @param fileName
#' Name or full path to JSON file to be created.
#' @param ...
#' other arguments (passed to \code{\link{asjson.ddsimca}}).
#'
#' @export
writeJSON.ddsimca <- function(obj, fileName, ...) {
   m <- asjson(obj, ...)
   fileConn <- file(fileName)
   writeLines(m, fileConn)
   close(fileConn)
}


################################
#  Plotting methods            #
################################

#' Sensitivity plot.
#'
#' @param obj
#' DD-SIMCA model (object of class \code{ddsimca}).
#' @param limType
#' limit type to show the plot for ('classic' or 'robust').
#' @param col
#' vector with two colors (for calibration and PV-set sensitivity).
#' @param type
#' type of the plot ("b", "l", or "h").
#' @param pch
#' vector with two markers (for calibration and PV-set sensitivity).
#' @param legend.position
#' position of legend on the plot.
#' @param ...
#' any parameters suitable for the \code{\link{plotSensitivity.ddsimcares}} method.
#'
#' @details
#' The method shows sensitivity vs. number of components for calibration and PV-set results.
#'
#' @export
plotSensitivity.ddsimca <- function(obj,
   limType = "classic",
   col = mdaplot.getColors(2),
   type = "b",
   pch = c(16, 16),
   legend.position = "bottomright",
   ...) {

   if (is.null(obj[["res"]]) || is.null(obj$res[["cal"]])) {
      message("Calibration results not found (most probably this model is loaded from web-application).")
      return(invisible(NULL))
   }

   p <- plotSensitivity.ddsimcares(obj$res[["cal"]], pch = pch[1], col = col[1], limType = limType, ...)
   if (!is.null(obj$res[["pv"]])) {
      stopifnot("Number of values for 'col' parameter should be 2." = length(col) == 2)
      stopifnot("Number of values for 'pch' parameter should be 2." = length(pch) == 2)

      pd <- plotSensitivity.ddsimcares(obj$res[["pv"]], limType = limType, show.plot = FALSE, ...)
      mdaplot(pd, show.axes = FALSE, type = type, pch = pch[2], col = col[2], ...)
      mdaplotg.showLegend(c("cal", "pv"), pch = pch, col = col, lty = NA, position = legend.position)
   } else {
      mdaplotg.showLegend("cal", pch = pch[1], col = col[1], lty = NA, position = legend.position)
   }
}


#' A shortcut to \code{\link{plotExtremes.ddsimca}}.
#'
#' @param obj
#' DD-SIMCA model (object of class \code{ddsimca}).
#' @param ...
#' any parameters suitable for \code{\link{plotExtremes.ddsimca}}.
#'
#' @export
plotExtreme.ddsimca <- function(obj, ...) {
   return(plotExtremes.ddsimca(obj, ...))
}

#' Extreme plot
#'
#' @description
#' Shows a plot with number of expected vs. number of observed extreme objects for different
#' significance levels (alpha values) for calibration and PV-set results.
#'
#' @param obj
#' a DD-SIMCA model (object of class \code{ddsimca})
#' @param ncomp
#' number of components to show the plot for
#' @param limType
#' limit type to show the plot for ('classic' or 'robust').
#' @param col
#' vector with two colors (for calibration and PV-set sensitivity).
#' @param legend.position
#' position of the legend
#' @param ...
#' other arguments, suitable for \code{\link{plotExtremes.ddsimcares}}
#'
#' @export
plotExtremes.ddsimca <- function(obj,
   ncomp = obj$ncomp.selected,
   limType = "classic",
   col = mdaplot.getColors(2),
   legend.position = "topleft",
   ...) {

   if (is.null(obj[["res"]]) || is.null(obj$res[["cal"]])) {
      message("Calibration results not found (most probably this model is loaded from web-application).")
      return(invisible(NULL))
   }

   p <- plotExtremes(obj$res[["cal"]], col = col[1], ncomp = ncomp, limType = limType, ...)

   if (!is.null(obj$res[["pv"]])) {
      stopifnot("Number of values for 'col' parameter should be 2." = length(col) == 2)

      pd <- plotExtremes(obj$res[["pv"]], limType = limType, ncomp = ncomp, show.plot = FALSE, ...)
      mdaplot(pd, show.axes = FALSE, col = col[2], ...)
      mdaplotg.showLegend(c("cal", "pv"), col = col, lty = NA, position = legend.position)
   } else {
      mdaplotg.showLegend("cal", col = col[1], lty = NA, position = legend.position)
   }
}

#' Acceptance plot for DD-SIMCA model.
#'
#' @description
#' Shows Acceptance plot for calibration of Procrustes validation results (if available).
#'
#' @param obj
#' DD-SIMCA model  (object of class \code{ddsimca}).
#' @param res
#' name of the results (either 'cal' or 'pv')
#' @param ...
#' any parameters suitable for \code{\link{plotAcceptance.ddsimcares}}.
#'
#' @export
plotAcceptance.ddsimca <- function(obj, res = "cal", ...) {

   res <- match.arg(res, c("cal", "pv"))

   if (is.null(obj[["res"]]) || is.null(obj$res[["cal"]])) {
      message("Calibration results not found (most probably this model is loaded from web-application).")
      return(invisible(NULL))
   }

   if (is.null(obj$res[[res]])) {
      message(sprintf("Result object '%s' not found.", res))
      return(invisible(NULL))
   }

   return(plotAcceptance(obj$res[[res]], res.name = res, ...))
}


#' Show with distance values (score, orthogonal or full) vs object indices for calibration
#' and PV-set results.
#'
#' @param obj
#' DD-SIMCA model (object of class \code{ddsimca})
#'
#' @param ncomp
#' number of components to show the plot for.
#' @param res
#' name of the results (can be 'cal', 'pv', 'both', or 'diff', in the latter case shows difference
#' of the distance values between the training and the PV-set).
#' @param limType
#' limit type to show the plot for ('classic' or 'robust')
#' @param distance
#' which distance to show (\code{"h"} for score, \code{"q"} for orthogonal or \code{"f"} for full distance).
#' @param log
#' logical, apply log transformation or not.
#' @param lim.lty
#' vector with two values - types of the lines showing critical limits (extreme, outlier).
#' @param lim.col
#' vector with two values - colors of the lines showing critical limits (extreme, outlier).
#' @param lim.lwd
#' vector with two values - thickness of the lines showing critical limits (extreme, outlier).
#' @param ylim
#' limits for y-axis, if NULL will be estimated automatically.
#' @param xlim
#' limits for x-axis, if NULL will be estimated automatically.
#' @param ...
#' other parameters compatible with \code{\link{plotDistances.ddsimcares}} method.
#'
#' @export
plotDistances.ddsimca <- function(obj, res = "both",
   ncomp = obj$ncomp.selected,
   limType = "classic",
   distance = "h",
   log = FALSE,
   lim.lty = c(2, 3),
   lim.col = c("darkgray", "darkgray"),
   lim.lwd = c(1, 1),
   ylim = NULL, xlim = NULL,
   ...) {

   if (is.null(obj[["res"]])) {
      stop("Object with results not found (most probably this model is loaded from web-application).", call. = FALSE)
   }

   res <- match.arg(res, c(names(obj$res), "both", "diff"))

   if (res == "both" || res == "diff") {
      plot_data <- lapply(obj$res, plotDistances.ddsimcares, limType = limType, ncomp = ncomp,
         distance = distance, log = log, show.plot = FALSE, ...)
   } else {
      plot_data <- list(plotDistances.ddsimcares(obj$res[[res]], limType = limType, ncomp = ncomp,
         distance = distance, log = log, show.plot = FALSE, ...))
      names(plot_data) <- res
   }

   if (distance != "f") {
      return(invisible(mdaplotg(plot_data, type = "p", ylim = ylim, ...)))
   }

   if (res == "diff") {
      x <- plot_data[[1]][, 1]
      y <- matrix(abs(plot_data[[1]][, 2] - plot_data[[2]][, 2]), nrow = 1)
      attr(y, "xaxis.values") <- x
      attr(y, "xaxis.name") <- colnames(plot_data[[1]])[1]
      attr(y, "yaxis.name") <- colnames(plot_data[[1]])[2]
      attr(y, "name") <- attr(plot_data[[1]], "name")
      rownames(y) <- "abs(cal - pv)"
      colnames(y) <- rownames(plot_data[[2]])
      return(invisible(mdaplotg(y, type = "h", ylab = colnames(plot_data[[1]])[2], ...)))
   }

   limType <- processLimType(limType)
   res <- obj$res$cal$simca
   v <- res$outcomes[[limType]]$values[[ncomp]]

   yle <- v$fce / v$Nf
   ylo <- v$fco / v$Nf

   if (log == TRUE) {
      yle <- log(1 + yle)
      ylo <- log(1 + ylo)
   }

   if (is.null(ylim)) {
      ymax <- max(vapply(plot_data, function(x) attr(x, "ylim")[2], 0))
      ylim <- c(0, ymax * 1.2)
   }

   p <- mdaplotg(plot_data, type = "p", ylim = ylim, ...)

   if (yle > 0 && !is.null(lim.lty) && !is.na(lim.lty[1])) {
      abline(h = yle, col = lim.col[1], lty = lim.lty[1], lwd = lim.lwd[1])
   }

   if (ylo > 0 && !is.null(lim.lty) && !is.na(lim.lty[2])) {
      abline(h = ylo, col = lim.col[2], lty = lim.lty[2], lwd = lim.lwd[2])
   }

   return(invisible(p))
}

