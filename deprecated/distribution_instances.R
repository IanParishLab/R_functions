distribution_instances <- function(df,col_name, nbreaks=10){
  y <- hist(df[[col_name]], breaks = nbreaks)
  plot(y, ylim=c(0, max(y$counts)+5))
  text(y$mids, y$counts+3, y$counts, cex=0.75)
  
}
