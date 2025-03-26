logfc <- function(x, y, base = 2){ return(mean(log(x, base)) - mean(log(y, base))) }

chi.test <- function(a, b) { return(chisq.test(cbind(a, b))) }
