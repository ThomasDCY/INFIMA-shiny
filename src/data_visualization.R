###### This file contains functions to plot INFIMA input data ######


# x should be an element in input_queries_small
plot_input <- function(x, option = NULL, ...) {
  
  if (is.null(option)) {
    stop('Please input the option argument!')
  }
  
  if (!option %in% c(1, 2, 3, 4)) {
    stop('Invalid value for the option argument!')
  }
  
  strains <- c('129', 'AJ', 'B6', 'CAST',
               'NOD', 'NZO', 'PWK', 'WSB')
  
  colors <- c(
    rgb(240, 128, 128, maxColorValue = 255, alpha = 255),
    rgb(218, 165, 32, maxColorValue = 255, alpha = 255),
    rgb(128, 128, 128, maxColorValue = 255, alpha = 255),
    rgb(0, 160, 0, maxColorValue = 255, alpha = 255),
    rgb(16, 16, 240, maxColorValue = 255, alpha = 255),
    rgb(0, 160, 240, maxColorValue = 255, alpha = 255),
    rgb(240, 0, 0, maxColorValue = 255, alpha = 255),
    rgb(144, 0, 224, maxColorValue = 255, alpha = 255)
  )
  
  if (option == 1) {
    grid <- par(mfrow = c(1, 2),
                mar = c(1, 1, 1, 1) + 3)
    plot(
      x$Y,
      main = '',
      xaxt = 'n',
      xlab = 'Strain',
      ylab = 'DO allele effect',
      col = colors,
      pch = 19
    )
    axis(1, at = 1:8, labels = strains, las = 2)
    
    plot(
      x$Y.t,
      main = '',
      xaxt = 'n',
      xlab = 'Strain',
      ylab = 'DO allele effect (trinarized)',
      col = colors,
      pch = 19,
      yaxt = 'n',
      ylim = c(-1,1)
    )
    axis(1, at = 1:8, labels = strains, las = 2)
    axis(2, at = c(-1,0,1), labels = c('low', 'middle', 'high'))
    par(grid)
    
  }
  
  if (option == 2) {
    grid <- par(mfrow = c(1, 2),
                mar = c(1, 1, 1, 1) + 3)
    plot(
      x$A,
      main = '',
      xaxt = 'n',
      xlab = 'Strain',
      ylab = 'Local ATAC-seq signal',
      col = colors,
      pch = 19
    )
    axis(1, at = 1:8, labels = strains, las = 2)
    
    plot(
      x$A.t,
      main = '',
      xaxt = 'n',
      xlab = 'Strain',
      ylab = 'Local ATAC-seq signal (trinarized)',
      col = colors,
      pch = 19,
      yaxt = 'n',
      ylim = c(-1,1)
    )
    axis(1, at = 1:8, labels = strains, las = 2)
    axis(2, at = c(-1,0,1), labels = c('low', 'middle', 'high'))
    par(grid)
  }
  
  if (option == 3) {
    grid <- par(mfrow = c(1, 2),
                mar = c(1, 1, 1, 1) + 3)
    plot(
      x$B.avg,
      main = '',
      xaxt = 'n',
      xlab = 'Strain',
      ylab = 'Founder gene expression',
      col = colors,
      pch = 19
    )
    axis(1, at = 1:8, labels = strains, las = 2)
    
    plot(
      x$B.t,
      main = '',
      xaxt = 'n',
      xlab = 'Strain',
      ylab = 'Founder gene expression (trinarized)',
      col = colors,
      pch = 19,
      yaxt = 'n',
      ylim = c(-1,1)
    )
    axis(1, at = 1:8, labels = strains, las = 2)
    axis(2, at = c(-1,0,1), labels = c('low', 'middle', 'high'))
    par(grid)
  }
  
  if (option == 4) {
    grid <- par(mfrow = c(1, 2),
                mar = c(1, 1, 1, 1) + 3)
    plot(
      x$E.t,
      main = '',
      xaxt = 'n',
      xlab = 'Strain',
      ylab = 'Founder allele effect (trinarized)',
      col = colors,
      pch = 19,
      yaxt = 'n',
      ylim = c(-1,1)
    )
    axis(1, at = 1:8, labels = strains, las = 2)
    axis(2, at = c(-1,0,1), labels = c('low', 'middle', 'high'))
    
    plot(
      x$D,
      main = '',
      xaxt = 'n',
      xlab = 'Strain',
      ylab = 'Edit distance',
      col = colors,
      pch = 19,
      yaxt = 'n'
    )
    axis(1, at = 1:8, labels = strains, las = 2)
    axis(2, at = -2:2)
    par(grid)
  }
}


# plot the Dirichlet prior information
# x should be a row from the prior data frame
plot_prior <- function(x, ...){
  rankscore <- as.numeric(x[, c('cor.A.E.rs', 'cor.A.B.rs', 'footprint.rs', 'dist.rs')])
  realscore <- as.numeric(x[, c('cor.A.E', 'cor.A.B', 'footprint', 'dist')])
  label <- c('cor(A,E)', 'cor(A,B)', 'footprint', 'distance')
  
  df <- data.frame(rankscore = rankscore,
                   realscore = realscore,
                   label = label)
  df <- df[order(df$rankscore),]
  
  rankscore <- df$rankscore
  names(rankscore) <- df$label
  
  bp <- barplot(rankscore, 
                xlab = '', 
                ylab = 'Rank score (compared to other candidate SNPs)', 
                main = 'Contribution of prior components\n (original values in brackets)',
                ylim = c(0,1),
                xaxt = 'n')
  realscore <- df$realscore
  text(bp, 0.1, labels = paste0('(', format(round(realscore, 2), nsmall = 2), ')'))
  axis(1, at = bp, labels = df$label, las = 2)
}

