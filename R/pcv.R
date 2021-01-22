#' Compute matrix with pseudo-validation set
#'
#' @param x
#' matrix with calibration set (IxJ)
#' @param ncomp
#' number of components for PCA decomposition
#' @param nseg
#' number of segments in cross-validation
#' @param scale
#' logical, standardize columns of X prior to decompositon or not
#'
#' @description
#' The method computes pseudo-validation matrix Xpv, based on PCA decomposition of calibration set X
#' and systematic (venetian blinds) cross-validation. It is assumed that data rows are ordered
#' correctly, so systematic cross-validation can be applied.
#'
#' All details can be found in [1]
#'
#' @references
#' 1. Kucheryavskiy, S., Zhilin, S., Rodionova, O., & Pomerantsev, A.
#' Procrustes Cross-Validation—A Bridge between Cross-Validation and Independent Validation Sets.
#' Analytical Chemistry, 92 (17), 2020. pp.11842–11850. DOI: 10.1021/acs.analchem.0c02175
#'
#' @return
#' Pseudo-validation matrix (IxJ)
#'
#' @export
pcv <- function(x, ncomp = min(round(nrow(x)/nseg) - 1, col(x), 20), nseg = 4, scale = FALSE) {

   # keep names if any
   attrs <- attributes(x)

   # convert to matrix if necessary
   x <- mda.df2mat(x)

   mx <- apply(x, 2, mean)
   sx <- if (scale) apply(x, 2, sd) else rep(1, ncol(x))

   # autoscale the calibration set
   x <- scale(x, center = mx, scale = sx)

   # create a global model
   P <- svd(x)$v[, seq_len(ncomp), drop = FALSE]

   # create matrix with indices for cross-validation, so
   # each column is number of rows to be taken as local validation set
   # in corresponding segment
   ind <- matrix(seq_len(nrow(x)), ncol = nseg, byrow = TRUE)

   # prepare empty matrix for pseudo-validation set
   x.pv <- matrix(0, nrow(x), ncol(x))
   a <- NULL

   # cv-loop
   for (k in seq_len(nseg)) {

      # split data to calibration and validation
      x.c <- x[-ind[, k], , drop = FALSE]
      x.k <- x[ ind[, k], , drop = FALSE]

      # get loadings for local model and rotation matrix between global and local models
      P.k <- svd(x.c, nv = ncomp)$v[, seq_len(ncomp), drop = FALSE]

      # correct direction of loadings for local model
      a <- acos(colSums(P * P.k)) < pi / 2
      P.k <- P.k %*% diag(a * 2 - 1, ncol(P), ncol(P))

      # get rotation matrix between the PC spaces
      R <- getR(P.k, P)

      # rotate the local validation set and save as a part of Xpv
      x.pv[ind[, k], ] <- tcrossprod(x.k, R)
   }

   # uscenter and unscale the data
   x.pv <- sweep(x.pv, 2, sx, "*")
   x.pv <- sweep(x.pv, 2, mx, "+")

   attributes(x.pv) <- attrs
   return(x.pv)
}

#' Creates rotation matrix to map a set vectors \code{base1} to a set of vectors \code{base2}.
#'
#' @param base1
#' Matrix (JxA) with A orthonormal vectors as columns to be rotated (A <= J)
#' @param base2
#' Matrix (JxA) with A orthonormal vectors as columns, \code{base1} should be aligned with
#'
#' @description
#' In both sets vectors should be orthonormal.
#'
#' @return
#' Rotation matrix (JxJ)
getR <- function(base1, base2) {
   base1 <- as.matrix(base1);
   base2 <- as.matrix(base2);

   R1 <- rotationMatrixToX1(base1[, 1])
   R2 <- rotationMatrixToX1(base2[, 1])

   if (ncol(base1) == 1) {
      R <- t(R2) %*% R1
   } else {
      # Compute bases rotated to match their first vectors to [1 0 0 ... 0]'
      base1_r <- as.matrix(R1 %*% base1)
      base2_r <- as.matrix(R2 %*% base2)

      # Get bases of subspaces of dimension n-1 (forget x1)
      nr <- nrow(base1_r) # equal to nrow(base2_r)
      nc <- ncol(base1_r) # equal to ncol(base2_r)
      base1_rs <- base1_r[2:nr, 2:nc]
      base2_rs <- base2_r[2:nr, 2:nc]

      # Recursevely compute rotation matrix to map subspaces
      Rs <- getR(base1_rs, base2_rs)

      # Construct rotation matrix of the whole space (recall x1)
      M <- eye(nr)
      M[2:nr, 2:nr] <- Rs

      R <- crossprod(R2, (M %*% R1))
   }

   return(R);
}

#' Creates a rotation matrix to map a vector x to [1 0 0 ... 0]
#'
#' @param x
#' Vector (sequence with J coordinates)
#'
#' @return
#' Rotation matrix (JxJ)
rotationMatrixToX1 <- function(x) {
   N <- length(x)
   R <- eye(N)
   step <- 1
   while (step < N) {
      A <- eye(N)
      n <- 1
      while (n <= N - step) {
         r2 <- x[n]^2 + x[n + step]^2
         if (r2 > 0) {
            r <- sqrt(r2)
            pcos <- x[n] / r
            psin <- -x[n + step] / r
            A[n, n] <- pcos
            A[n, n + step] <- -psin
            A[n + step, n] <- psin
            A[n + step, n + step] <- pcos
         }
         n <- n + 2 * step
      }
      step <- 2 * step
      x <- A %*% x
      R <- A %*% R
   }
   return(R)
}


#' Create the identity matrix
#'
#' @param n
#' Size of the matrix
#'
#' @return
#' The identity matrix (n x n)
#'
#' @export
eye <- function(n) {
   X <- matrix(0, n, n)
   diag(X) <- 1
   return(X)
}

