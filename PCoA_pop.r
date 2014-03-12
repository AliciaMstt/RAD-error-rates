## PCoA display by POPULATION and defined colors
PCoA_pop = function (dat.d, vect.grps = 0, cols.pop.key, pile= FALSE) {
  # dat.d is the distance matrix
  # cols.pop.key is a list of colors specific to each pop, defined externally
  # set pile to TRUE to cancel par(mfrow = c(1, 2) and allow to pile the results of the function run with different data
  vect.grps = as.factor(as.character(vect.grps))
  num.grps <- length(table(vect.grps))
  names.grps <- levels(vect.grps)
  cols.pop.key = cols.pop.key #cols.pop.key is a list of colors specific to each pop, defined externally
  topo.col <- cols.pop.key[match(levels(vect.grps), cols.pop.key[, 1]), 2] 
  # perform classical multidimensional sclaing (PCoA) of the dist matrix  
  # k is the number of samples - 5 and eigenvalues are returned
  acop <- cmdscale(dat.d, k = nrow(as.matrix(dat.d)) - 1, eig = TRUE) # keep n-1 eigenvalues
  axes.tot <- acop$eig # eig are the n eigenvalues computed by cmdscale. Axes are ranked by their eigenvalues, so the first axis has the highest eigenvalue, the second axis has the second highest eigenvalue, etc.
  inertia <- sum(acop$eig[acop$eig > 0]) 
  percents <- round(as.numeric(100 * axes.tot/inertia), digits = 0) # The eigenvalues represent the variance extracted by each axis, here they are expressed as a percentage of the sum of all eigenvalues (i.e. total variance).
  par()
  
  if(pile == TRUE){
    par(pty = "s") # and dont forget to define mfrow=(x,2) externally
  } else {
    par(mfrow = c(1, 2), pty = "s") 
  }
  coord1 <- acop$points[, c(1, 2)] # points is a matrix whose rows give the coordinates of the points chosen to represent the dissimilarities
  col.grps <- data.frame(vect.grps, coord1)
  # plot so that the maximum variance is projected along the first axis, then on the second and so on
  plot(coord1, asp = 1, cex = 0.1, xlab = paste("Axis 1 ", 
    "(", percents[1], " % variance explained)", sep = ""), 
    ylab = paste("Axis 2 ", "(", percents[2], " % variance explained)", 
      sep = ""), main = "", type = "n", bty = "n")
  abline(h = 0, lty = 2, col = "grey")
  abline(v = 0, lty = 2, col = "grey")
  if (length(vect.grps) == nrow(as.matrix(dat.d))) {
    for (g in 1:length(names.grps)) {
      text(x = coord1[col.grps[, 1] == names.grps[g], 1], 
        y = coord1[col.grps[, 1] == names.grps[g], 2], 
        labels = names.grps[g], col = topo.col[g], cex = 0.7)
    }
  }
  else {
    points(coord1, pch = 19, col = "blue", cex = 0.5)
  }
  coord1 <- acop$points[, c(3, 4)]
  col.grps <- data.frame(vect.grps, coord1)
  plot(coord1, asp = 1, cex = 0.1, xlab = paste("Axis 3 ", 
    "(", percents[3], " % variance explained)", sep = ""), 
    ylab = paste("Axis 4 ", "(", percents[4], " % variance explained)", 
      sep = ""), main = "", type = "n", bty = "n")
  abline(h = 0, lty = 2, col = "grey")
  abline(v = 0, lty = 2, col = "grey")
  if (length(vect.grps) == nrow(as.matrix(dat.d))) {
    for (g in 1:length(names.grps)) {
      text(x = coord1[col.grps[, 1] == names.grps[g], 1], 
        y = coord1[col.grps[, 1] == names.grps[g], 2], 
        labels = names.grps[g], col = topo.col[g], cex = 0.7)
    }
  }
  else {
    points(coord1, pch = 19, col = "blue", cex = 0.5)
  }
}
