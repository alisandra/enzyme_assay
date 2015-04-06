#' Convert to Seconds
#'
#' This function converts time in the format exported by the plate reader (H:MM:SS) to seconds
#' @param time a value or vector of strings for time points in the format "H:MM:SS"
#' @export
#' @examples
#' plate_reader_seconds(enzyme_assay.time)

plate_reader_seconds <- function(time){
        time <- sapply(time,as.character)
        total_s <- sapply(time,function(x){
                splittime <- as.numeric(unlist(strsplit(x, split=":")))
                h <- splittime[1]
                m <- splittime[2]
                s <- splittime[3]
                out <- s + 60 * m + 60^2 * h
                return(out)
        })
        return (total_s)
}

