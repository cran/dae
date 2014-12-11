fac.match <- function(x, table, col.names)
#Function to find, for each combination of col.names in x, the row that has the same combination in table
#It can be viewed as a generalization of the match function from a single vector to multiple vectors
{ if (class(col.names) != "character") 
    stop("Must supply a character vector of column names")
  ncols <- length(col.names)
  if (any(!(col.names %in% names(x))))
    stop("All column names must be in x")
  if (any(!(col.names %in% names(table))))
    stop("All column names must be in x")
  if (length(dim(x)) != 2 | length(dim(table)) != 2)
    stop("Both x and table must have two dimensions")
  col.list <- as.list(x[col.names])
  index <- as.data.frame.table(with(x, by(x, col.list, 
                                          function(x, table, col.names)
                                          { k <- rep(TRUE, nrow(table))
                                            i <- 1
                                            for (i in 1:ncols)
                                            { k <- k & (table[[col.names[i]]] == x[[col.names[i]]])
                                            }
                                            if (sum(k) != 1)
                                              stop("Some combinations in x have more than one combination in table")
                                            else
                                              k <- c(1:nrow(table))[k]
                                            return(k)
                                          },
                                          table = table, col.names = col.names)))$Freq
  return(index)
}
