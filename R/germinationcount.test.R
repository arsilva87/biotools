germinationcount.test <- 
function(r, nsamples, n, N, K) {
   p.val <- 1 - pRangeHyper(r, nsamples, N, K, n)
   offrange <- r
   germrate <- K/N
   out <- list(R.value = offrange, p.value = p.val,
       germination.rate = germrate)
   class(out) <- "germinationcount.test"
   return(out)
}

# -----------------------------------------
# print method
print.germinationcount.test <- function(x, digits = 4, ...) {
    cat("\n            Germination count range test for seed sample heterogeneity",
      "\n\nGermination rate of the seed lot:", 100*x$germination.rate, "%",
	 "\nR-value (count difference):", x$R.value, "seeds",
	 "\np-value: ", x$p.value)
    cat("\nNull hypothesis: germination homogeneity", "\n")
    invisible(x)
}

# ------------------------------------------
# theoretical distribution based of count range
dRangeHyper <- function(r, size, N, K, n) 
{
   sapply(r, function(r) {
      if(r < 0) {
         prob <- 0
      } else if(r == 0) {
         prob <- sum(dhyper(0:n, K, N-K, n)^size)
      } else {
         prob <- 0
         for(x in 0:n) {
            prob <- prob +  
               (phyper(x+r, K, N-K, n) - phyper(x-1, K, N-K, n))^size - 
               (phyper(x+r, K, N-K, n) - phyper(x, K, N-K, n))^size -
               (phyper(x+r-1, K, N-K, n) - phyper(x-1, K, N-K, n))^size +
               (phyper(x+r-1, K, N-K, n) - phyper(x, K, N-K, n))^size
         }
      }
      return(prob)
   })
}

# CDF --------------------------------------
pRangeHyper <- function(r, size, N, K, n)
{
   sapply(r, function(r) {
      sum(dRangeHyper(0:r, size, N, K, n))  
   })
}

