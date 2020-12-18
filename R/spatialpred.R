spatialpred <- function(coords, data, grid)
{
   if (!requireNamespace("tcltk", quietly = TRUE)) 
        stop("package tcltk is required")
   y <- data
   n <- nrow(coords)
   N <- nrow(grid)
   pred <- data.frame(grid, pred = NA, SE = NA, radius = NA, n = NA)

   # progress bar... for the loop structure
   pb <- tcltk::tkProgressBar(title = "Spatial Interpolation", 
       label = "CALCULATION PROGRESS", min = 0, max = N, width = 400L)
   
   # loop for interpolation
   i = 1
   repeat{
      s <- sweep(coords, 2, as.matrix(grid[i, ]), FUN = "-")
      d <- sqrt(apply(s^2, 1, sum))
      w <- 1/d^2
      max.w <- max(w[is.finite(w)])
      w[!is.finite(w)] <- max.w

      # selecting the radius for the i-th grid point
      # criterion: data causing the least variability (IV)
      seq.d <- seq(min(d), max(d), length.out = 200)
      vari <- c()
      for(j in 1:length(seq.d)) {
         vari[j] <- iv( y[which(d <= seq.d[j])] )
      }
      pred$radius[i] <- seq.d[which.min(vari)]

      id <- which(d <= pred$radius[i])
      pred$pred[i] <- weighted.mean(y[id], w[id])
      pred$n[i] <- length(y[id])
      pred$SE[i] <- sd(y[id])/sqrt(pred$n[i])

      tcltk::setTkProgressBar(pb, i, 
         label = sprintf("CALCULATION PROGRESS (%.0f%%)", 
         100 * i/N))
      i = i + 1
      if (i > N) break()
   }

   # percentual absolut mean error
   if (all(coords == grid)) {
      pame <- 100*mean(abs(data - pred$pred)/data)
      cat("PAME:", round(pame, 2), "%\n")
   }
   Sys.sleep(0.5)
   close(pb)
   return(pred)
}

# aux functions
# index of variation: CV%/sqrt(n)
iv <- function(x) sd(x)/( mean(x) * sqrt(length(x)) ) 


