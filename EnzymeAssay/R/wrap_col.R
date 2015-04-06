#' Replicate Vector to N
#'
#' This function extends a vector of colors to specified length.
#' @param col A vector of colors.
#' @param n The desired length for the returned vector of colors.
#' @export
#' @examples
#' wrap_col(c('red','blue'), 8)

wrap_col <- function(col, n){
        newcol <- col
        while(length(newcol) < n){
                newcol <- c(newcol, col)
        }
        if (length(newcol) != n){
                warning('n not divisible by col, truncating\n')
        }
        return(newcol)
}

