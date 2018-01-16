################################################################################
# Truncated Piecewised Function (tp)
# Author: Yu Liu
# Date: June 2017
################################################################################
tp <- function(x, name, knots = NULL, nknots = NULL, raw = FALSE)
{
  # x atomic vector return data frame
  # Check if there is missing data
  if (any(is.na(x))) stop('Missing data in Variables')
  # Create Knots
  if(is.null(knots)){
    if(!is.null(nknots)){
      if(nknots > min(as.integer(length(unique(x))/4),30)){
        stop("Number of Knots is larger than thumb of rules ")
      }
    } else{
      nknots <- min(as.integer(length(unique(x))/4),30)
    }
    knots <- quantile(
      x = x, probs = seq(0, 1, length.out = nknots + 2)
    )[-c(1, nknots + 2)]
  }

  basis_matrix = sapply(knots, function(knot) pmax(x - knot, 0))
  if(raw){
    basis_matrix <- cbind(x, basis_matrix)
    colnames(basis_matrix) <- sapply(0: nknots, function(i) paste0(name,i))
  } else{
    colnames(basis_matrix) <- sapply(1: nknots, function(i) paste0(name,i))
  }
  return(as.data.frame(basis_matrix))
}
