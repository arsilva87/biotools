mvpaircomp <- 
function(model, factor1, nesting.factor = NULL,
   test = "Pillai", adjust = "none", SSPerror = NULL, DFerror = NULL) 
{
   stopifnot(inherits(model, "mlm"))
   test <- match.arg(test, c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"))
   adjust <- match.arg(adjust, p.adjust.methods)
   W <- if(is.null(SSPerror)) SSD(model)$SSD else SSPerror
   dfe <- if(is.null(DFerror)) SSD(model)$df else DFerror
   y <- model$model[[1]]
   y.names <- colnames(y)
   fac <- model$model
   lev <- model$xlevels[[factor1]]
   pares <- combn(lev, 2)
   nc <- ncol(pares)
   if( is.null(nesting.factor) ) {
      f2 <- ""
      nreps <- cbind(xtabs(~ fac[[factor1]]))
   } else {
      f2 <- model$model[nesting.factor]
      nreps <- xtabs(~ fac[[factor1]] + fac[[nesting.factor]])
   }
   dat <- data.frame(y, f1 = fac[[factor1]])
   dl <- split(dat, f = f2) 
   # loop ------------
   st <- array(dim = c(nc, 5, length(dl)), 
      dimnames = list(paste(pares[1, ], "-", pares[2, ]), 
         c(test, "approx F", "num DF", "den DF", "Pr(>F)"), 
         names(dl)))
   B <- array(dim = c(nrow(W), ncol(W), nc, length(dl)),
      dimnames = list(rownames(W), colnames(W), 
         paste(pares[1, ], "-", pares[2, ]), names(dl)))
   for(j in 1:length(dl)) {
      med <- aggregate( .~ f1, data = dl[[j]], FUN = mean)
      rownames(med) <- med[, "f1"]
      med <- med[, y.names]
      nrep <- nreps[, j]
      for(i in 1:nc) {
         d.matrix <- apply(med[pares[, i], ], 2, diff)
         a <- 1/sum(1/nrep[pares[, i]])
         B[,, i,j] <- crossprod(t(d.matrix)) * a
         eigs <- eigen(solve(W, B[,, i,j]), only.values = TRUE)$values
         if(is.complex(eigs)) eigs <- makeReal(eigs)
         st[i, 1:4, j] <- switch(test, 
             Pillai = bioPillai(eigs, 1, dfe), 
             Wilks = bioWilks(eigs, 1, dfe),
             'Hotelling-Lawley' = bioHL(eigs, 1, dfe), 
             Roy = bioRoy(eigs, 1, dfe))
         pval <- pf(st[i, 2, j], df1 = st[i, 3, j], 
             df2 = st[i, 4, j], lower.tail = FALSE)
         st[i, 5, j] <- ifelse(adjust == "none", pval, 
             p.adjust(pval, adjust, n = nc))
      }
   }
   out <- list(st = st, SSPcontrast = B, adjust = adjust, 
      fac1 = factor1, fac2 = nesting.factor)
   class(out) <- "mvpaircomp"
   return(out)
}

# -----------------------------------------
makeReal <- function(x) if(all(Im(z <- zapsmall(x))==0)) as.numeric(z) else x

fstr <- function(x) {
   structure(as.data.frame(x), 
      heading = names(x), class = c("anova", "data.frame"))
}

print.mvpaircomp <- function(x, ...) {
   cat("\n            Multivariate Pairwise Comparisons\n\n")
   if(!is.null(x$fac2))
      cat("---\nComparing levels of", x$fac1, "nested in", x$fac2, "\n\n")
   if(dim(x$st)[3L] > 1) {
      k <- apply(x$st, 3, fstr)
      print(k, ...)
   } else {
       k <- fstr(x$st)
       colnames(k) <- colnames(x$st)
       print(k, ...)
   }
   cat("With", x$adjust, "p-value adjustment for multiple comparisons\n")
}

bioPillai <- function (eig, q, df.res) 
{
    test <- sum(eig/(1 + eig))
    p <- length(eig)
    s <- min(p, q)
    n <- 0.5 * (df.res - p - 1)
    m <- 0.5 * (abs(p - q) - 1)
    tmp1 <- 2 * m + s + 1
    tmp2 <- 2 * n + s + 1
    c(test, (tmp2/tmp1 * test)/(s - test), s * tmp1, s * tmp2)
}

bioWilks <- function (eig, q, df.res) 
{
    test <- prod(1/(1 + eig))
    p <- length(eig)
    tmp1 <- df.res - 0.5 * (p - q + 1)
    tmp2 <- (p * q - 2)/4
    tmp3 <- p^2 + q^2 - 5
    tmp3 <- if (tmp3 > 0) 
        sqrt(((p * q)^2 - 4)/tmp3)
    else 1
    c(test, ((test^(-1/tmp3) - 1) * (tmp1 * tmp3 - 2 * tmp2))/p/q, 
        p * q, tmp1 * tmp3 - 2 * tmp2)
}

bioHL <- function (eig, q, df.res) 
{
    test <- sum(eig)
    p <- length(eig)
    m <- 0.5 * (abs(p - q) - 1)
    n <- 0.5 * (df.res - p - 1)
    s <- min(p, q)
    tmp1 <- 2 * m + s + 1
    tmp2 <- 2 * (s * n + 1)
    c(test, (tmp2 * test)/s/s/tmp1, s * tmp1, tmp2)
}

bioRoy <- function (eig, q, df.res) 
{
    p <- length(eig)
    test <- max(eig)
    tmp1 <- max(p, q)
    tmp2 <- df.res - tmp1 + q
    c(test, (tmp2 * test)/tmp1, tmp1, tmp2)
}
