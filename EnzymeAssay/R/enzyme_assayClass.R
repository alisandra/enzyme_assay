#  Enzyme_AssayClass.R
#
#  Copyright (C) 2010-2014 Alisandra Kaye Denton,
#  Institute for Plant Biochemistry, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: adenton@bio1.rwth-aachen.de
#
#  This file is part of EnzymeAssay.
#
#  EnzymeAssay is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  EnzymeAssay is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with EnzymeAssay.  If not, see <http://www.gnu.org/licenses/>.


# Enzyme_Assay


#------------------------------------------------------------------------------#
#                            class definitions                                 #
#------------------------------------------------------------------------------#

#' Organizing Enzyme Assay Data
#'
#' This class holds, and facilitates the curation of enzyme assay data
#' @slot data A matrix with measurements in rows for different wells in columns
#' @slot seconds The second values corresponding to the rows in slot 'data'
#' @slot auto_slope The automatically calculated maximum slopes from \code{\link{enzyme_assay}}
#' @slot auto_range The ranges corresponding to slopes in slot 'auto_slope'
#' @slot curated_slope The user curated slopes from \code{\link{enzyme_assay.curate}} or \code{\link{enzyme_assay.curate_batch}}
#' @slot curated_range The ranges corresponding to slot 'curated_slope'
#' @export 
#' @examples 
#' ea <- enzyme_assay(enzyme_assay.data, enzyme_assay.time)
#' eaCurated <- enzyme_assay.curate(ea)
#' new_range <- ea@@auto_range
#' new_range[1,] <- 1
#' eaBatchCurated <- enzyme_assay.curate_batch(ea, new_range)

setClass('Enzyme_Assay',
        representation(
                        data = "matrix", 		# raw assay data
                        seconds = "vector",		# time for measurements
                        auto_slope = "vector",		# slope of automatically selected ranges
                        auto_range = "matrix",		# automatically selected ranges
                        curated_slope = "vector",	# slope of curated ranges
                        curated_range = "matrix"),	# curated ranges
        prototype(
                        data = matrix(),
                        seconds = vector(),
                        auto_slope = vector(),
                        auto_range = matrix(),
                        curated_slope = vector(),
                        curated_range = matrix())
)

#------------------------------------------------------------------------------#
#                              user constructors                               #
#------------------------------------------------------------------------------#

#' Getting Maximal Slopes From Enzyme Assay
#'
#' This functions generates an object of Enzyme_Assay class from raw data
#' @param raw_data A matrix of measurements (rows) and wells (columns).
#' @param x The time measurements corresponding to each row, either in seconds, or formatted "H:MM:SS"
#' @param seconds Whether or not x is given in seconds. Defaults to FALSE
#' @param decreasing Whether or not assay is decreasing. Defaults to FALSE
#' @param windowsize Length of ranges to consider.
#' @param method Method for calculating range of maximum slope, "individual" for well by well, "together" for assay mean
#' @export
#' @examples
#' enzyme_assay(enzyme_assay.data, enzyme_assay.time)


# init function
enzyme_assay <- function(raw_data, x, seconds=FALSE, decreasing = FALSE, windowsize = 8, method = c("individual", "together")){
        raw_data<-apply(raw_data,c(1,2),as.numeric)
        method <- match.arg(method)
        if (! seconds){
                s <- plate_reader_seconds(x)
        }else{
                s <- sapply(x,as.numeric)
        }
        # set whether we will look for the highest or lowest slopes
        if (decreasing){
                direction <- -1
        }else{
                direction <- 1
        }

        # calculate auto ranges and slope for each well individually
        if (method == "individual"){
                auto_slope <- rep(NA, ncol(raw_data))
                high_coordinates <- matrix(data=NA, nrow=2, ncol=ncol(raw_data))
                for (i in 1:ncol(raw_data)){
                        temp.slope <- 0
                        temp.coord <- c(1,nrow(raw_data))
                        # go through each position
                        for (j in 1:(nrow(raw_data) - windowsize)){
                                linreg <- lm(raw_data[j:(j + windowsize), i] ~ s[j:(j + windowsize)])
                                # [[1]][2]] parses the slope out of the regression object
                                new.slope <- linreg[[1]][[2]]
                                if (is.na(new.slope)){
                                        new.slope <- 0
                               }
                                if (new.slope * direction > temp.slope * direction){  # tests whether it's the highest slope yet observed
                                        temp.slope <- new.slope
                                        temp.coord <- c(j, j + windowsize)
                                }
                        }
                        auto_slope[i] <- temp.slope
                        high_coordinates[,i] <- temp.coord
                }
	# calculate auto ranges and slope for assay mean
        }else if (method == "together"){
                high_coordinates <- matrix(data=NA, nrow=2, ncol=ncol(raw_data))
                temp.slope <- 0
                for (j in 1:(nrow(raw_data) - windowsize)){
                        slopes <- c()
                        for (i in 1:ncol(raw_data)){
                                linreg <- lm(raw_data[j:(j + windowsize), i] ~ s[j:(j + windowsize)])
                                slopes <- c(slopes, linreg[[1]][[2]])
                        }
                        new.slope <- mean(slopes)
                        if (new.slope * direction > temp.slope * direction){  # tests whether it's the highest slope yet observed
                                temp.slope <- new.slope
                                temp.coord <- c(j, j + windowsize)
                                auto_slope <- slopes
                        }
                }
                high_coordinates[1,] <- temp.coord[1]
                high_coordinates[2,] <- temp.coord[2]
        }
        #match names
        names(auto_slope) <- colnames(raw_data)
        colnames(high_coordinates) <- names(auto_slope)
        rownames(high_coordinates) <- c('start','stop')
        out <- new('Enzyme_Assay', data=as.matrix(raw_data), seconds=s, auto_slope=auto_slope, auto_range=as.matrix(high_coordinates))
        return(out)
}

#------------------------------------------------------------------------------#
#                              curation methods                                #
#------------------------------------------------------------------------------#

# interpret user input: seconds to measurement
#' Finding Measurement Closest to Seconds
#'
#' This function returns the measurement closest to a selected seconds value
#' @param this An object of the Enzyme_Assay class
#' @param s A number in seconds
#' @export
#' @examples
#' ea <- enzyme_assay(enzyme_assay.data, enzyme_assay.time)
#' enzyme_assay.s_to_measurement(ea, 111)

setGeneric("enzyme_assay.s_to_measurement", function(this, s) standardGeneric("enzyme_assay.s_to_measurement"))
setMethod('enzyme_assay.s_to_measurement','Enzyme_Assay',
        function(this, s){
                #absolute differences between every second entry and 's' value of interest
                second_diff <- abs(this@seconds - s) 
                #where was the lowest difference
                index_min_s <- which(second_diff == min(second_diff))
                #take first (lower) on tie
                if (length(index_min_s) > 1){
                        index_min_s <- index_min_s[1]
                }
                return(index_min_s)
        }
)

# get and handle user input
#' Handling User Input
#'
#' This internal function asks the user to curate start or stop measurement for determining slope in an enzyme assay
#' @param object An object of the Enzyme_Assay class
#' @param oldvalue Previous start or stop value
#' @export
#' @examples
#' ea <- enzyme_assay(enzyme_assay.data, enzyme_assay.time)
#' eaCurated <- enzyme_assay.curate(ea)
setGeneric("enzyme_assay.getinput", function(object, ...) standardGeneric("enzyme_assay.getinput"))
setMethod('enzyme_assay.getinput','Enzyme_Assay',
        function(object, oldvalue){
                input <- readline()
                numinput <- suppressWarnings(as.numeric(input))
                continue <- 'continue'
                if (input == 'c'){
                        continue <- 'cancel'
                        out <- oldvalue
                }else if(input == 'k'){
                        out <- oldvalue
                }else if (! is.na(numinput)){
                        out <- enzyme_assay.s_to_measurement(object, numinput)
                        cat('rounded to nearest measurement:', out, '\n')
                }else{
                        out <- oldvalue
                        print('non-numeric value given, returning previous value untouched')
                }
                return(list(out,continue))
        }
)

# by column curation
#' Curate One Column
#'
#' This internal function handles curation of one well (column) of an enzyme assay
#' @param object An object of the Enzyme_Assay class
#' @param column The curent column
#' @param status 'curated' or 'auto', whether the initial values are taken from slot 'auto_range' or 'curated_range'
#' @param at_once The number of columns/wells to plot at once
#' @export
#' @examples
#' ea <- enzyme_assay(enzyme_assay.data, enzyme_assay.time)
#' eaCurated <- enzyme_assay.curate(ea)

setGeneric("enzyme_assay.subcurate", function(object, ...) standardGeneric("enzyme_assay.subcurate"))
setMethod('enzyme_assay.subcurate','Enzyme_Assay',
        function(object, column, status = c("curated", "auto"), at_once=12){
                #method <- match.arg(method) #and in the parameters 'method = c("individual", "together"),'
		method <- "individual" #relic, of the intent to implement a together option.
                status <- match.arg(status)
                if (status == "curated"){
                        range <- object@curated_range
                }else if (status =="auto"){
                        range <- object@auto_range
                }
                if (method == "individual"){
                        if (column <= ncol(range)){
                                out <- range[,column]
                                cname <- colnames(range)[column] # the column name
                                startsec <- object@seconds[range[1,column]] # the second value where previous slope measurement starts
                                endsec <- object@seconds[range[2,column]] # the second value where previous slope measurement ends
                                # get and handle input for starting position
                                cat("for", cname, status, "start =", startsec, "seconds. Please enter 'k' to keep, 'c' to cancel, or a new start in seconds: ")
                                processed_input <- enzyme_assay.getinput(object, input, range[1,column])
                                out[1] <- processed_input[[1]]
                                continue <- processed_input[[2]]
                                if (continue != 'cancel'){
                                        # get and handle input for stopping position
                                        cat("for", cname, status, "end =", endsec, "seconds. Please enter 'k' to keep, 'c' to cancel, or a new end in seconds: ")
                                        processed_input <- enzyme_assay.getinput(object, input, range[2,column])
                                        out[2] <- processed_input[[1]]
                                        continue <- processed_input[[2]]
                                }
                        }
                }
                return(list(out,continue))
        }
)


#calculate slope from curated range
#' Getting Slope from Range
#'
#' This internal function calculates the slopes associated with curated ranges
#' @param object An object of the Enzyme_Assay class
#' @export
#' @examples
#' ea <- enzyme_assay(enzyme_assay.data, enzyme_assay.time)
#' eaCurated <- enzyme_assay.curate(ea)

setGeneric("enzyme_assay.calc_slope", function(object) standardGeneric("enzyme_assay.calc_slope"))
setMethod('enzyme_assay.calc_slope','Enzyme_Assay',
        function(object){
                d.cur <- dim(object@curated_range)
                d.auto <- dim(object@auto_range)
                same_size <- all(d.cur == d.auto)
                if (same_size){
                        newslopes <- object@auto_slope
                        for (i in 1:ncol(object@curated_range)){
                                start <- object@curated_range[1,i]
                                end <- object@curated_range[2,i]
                                linreg <- lm(object@data[start:end, i] ~ object@seconds[start:end])
                                newslopes[i] <- linreg[[1]][[2]]
                        }
                        object@curated_slope <- newslopes

                }else{
                        warning('size of curated_range does not match auto_range in object, please curate properly')
                        out <- matrix()
                }
                return(object)

        }
)
 
#calculate slope from curated range
#' Well-By-Well Curation
#'
#' This function plots and asks the user to curate an Enzyme_Assay object well by well
#' @param x An object of the Enzyme_Assay class
#' @param at_once Number of wells/columns to plot at once
#' @param status Whether to start from 'auto' or 'curated' ranges
#' @export
#' @examples
#' ea <- enzyme_assay(enzyme_assay.data, enzyme_assay.time)
#' eaCurated <- enzyme_assay.curate(ea)


#curation
setGeneric("enzyme_assay.curate", function(x, ...) standardGeneric("enzyme_assay.curate"))
setMethod('enzyme_assay.curate','Enzyme_Assay',
        function(x,
                   ylim, xlim,
                   xlab = "Seconds",
                   ylab = "A 340",
                   pch = 20,
                   col = "black",
                   at_once = 12,
                   status = c("auto","curated")
                   #method = c("individual","together")
		){
                status = match.arg(status)
                #method = match.arg(method)
                if (missing(ylim)) {
                        ylim <- range(x@data)
                        yadj <- 0.05 * max(abs(ylim))
                        ylim <- c(ylim[1] - yadj, ylim[2] + yadj)
                }
                if (missing(xlim)){
                        xlim <- c(1,max(x@seconds))
                }
                if (status == "auto"){
                        ranges <- x@auto_range
                }else if (status == "curated"){
                        ranges <- x@curated_range
                }
                newcol = wrap_col(col=col, n=at_once)

                x@curated_range <- ranges

                par(ask=T)
                nplots <- ceiling(ncol(x@data)/at_once)
                # plot data in sets of reasonable size (default, 12, for one row on plate reader)
                my_grey = 'grey60'
                i <- 1
                while(i <= nplots){
                        plot(NA, ylim=ylim, xlab=xlab, ylab=ylab, pch=pch, xlim=xlim)
                        start <- (i - 1) * at_once
                        if (i < nplots){
                                all_j <- 1:at_once + start
                        }else{
                                all_j <- (start+1):ncol(x@data)
                        }
                        #preplot everything in light grey for perspective
                        for (j in all_j){
                                lines(x@seconds, x@data[,j], col = my_grey)
                                highlight <- ranges[1,j]:ranges[2,j]
                                lines(x@seconds[highlight], x@data[highlight,j],col = my_grey, lwd=3)
                        }
                        for (j in all_j){
                                # run only if the condition (i < nplots) of the outer for loop is met
                                if (i <= nplots){
                                        lines(x@seconds, x@data[,j], col = newcol[j-start])
                                        highlight <- ranges[1,j]:ranges[2,j]
                                        lines(x@seconds[highlight], x@data[highlight,j],col = newcol[j-start], lwd=3)
                                        subcuration <- enzyme_assay.subcurate(object=x, column=j, status=status, at_once=at_once)
                                        newranges <- subcuration[[1]]
                                        #stop the while and for loops (by setting i) if the user canceled in the subcurationfunction
                                        if (subcuration[[2]] == 'cancel'){
                                                i <- nplots + 1
                                        }
                                        x@curated_range[,j] <- newranges
                                        # 'erase' old highlight (cover with white)
                                        lines(x@seconds[highlight], x@data[highlight,j],col = 'white', lwd=3)
                                        # plot new in grey 
                                        lines(x@seconds, x@data[,j], col = my_grey)
                                        lines(x@seconds[newranges[1]:newranges[2]], x@data[newranges[1]:newranges[2],j],col = my_grey, lwd=3)
                                }
                        }
                        i <- i + 1

                }
                par(ask=F)
                return(x)
        }
)


#batch curation
#' Batch Curation
#'
#' This function replaces slot 'curated_range' with user supplied ranges, and recalculates 'curated_slope'
#' @param object An object of the Enzyme_Assay class
#' @param new_Range Matrix with start values for each column in slot 'data' in row one and end values in row two
#' @export
#' @examples
#' ea <- enzyme_assay(enzyme_assay.data, enzyme_assay.time)
#' new_range <- ea@@auto_range
#' new_range[1,] <- 1
#' eaCurated <- enzyme_assay.curate_batch(ea, new_range)

setGeneric("enzyme_assay.curate_batch", function(object, new_range) standardGeneric("enzyme_assay.curate_batch"))
setMethod('enzyme_assay.curate_batch','Enzyme_Assay',
        function(object, new_range){
                d.new <- dim(new_range)
                d.old <- dim(object@auto_range)
                same_size <- all(d.new == d.old)
                if (same_size){
                        object@curated_range <- new_range
                        object <- enzyme_assay.calc_slope(object)
                }else{
                        warning('size of new_range does not match slot auto_range of object')
                }
                return(object)
        }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#


setMethod('plot', signature(x = 'Enzyme_Assay', y = 'missing'),
        function(x,
                   ylim, xlim,
                   xlab = "Seconds",
                   ylab = "A 340",
                   pch = 20,
                   col = "black",
                   at_once = 12,
                   status = c("auto","curated")
                   ) {
                status = match.arg(status)
                if (missing(ylim)) {
                        ylim <- range(x@data)
                        yadj <- 0.05 * max(abs(ylim))
                        ylim <- c(ylim[1] - yadj, ylim[2] + yadj)
                }
                if (missing(xlim)){
                        xlim <- c(1,max(x@seconds))
                }
                if (status == "auto"){
                        ranges <- x@auto_range
                }else if (status == "curated"){
                        ranges <- x@curated_range
                }
                newcol = wrap_col(col=col, n=at_once)

                par(ask=T)
                nplots <- ceiling(ncol(x@data)/at_once)
                # plot data in sets of reasonable size (default, 12, for one row on plate reader)
                for (i in 1:(nplots - 1)){
                        plot(NA, ylim=ylim, xlab=xlab, ylab=ylab, pch=pch, xlim=xlim)
                        start = (i - 1) * at_once
                        for (j in 1:at_once + start){
                                lines(x@seconds, x@data[,j], col = newcol[j-start])
                                highlight <- ranges[1,j]:ranges[2,j]
                                lines(x@seconds[highlight], x@data[highlight,j],col = newcol[j-start], lwd=3)
                        }
                }
                plot(NA, ylim=ylim, xlab=xlab, ylab=ylab, pch=pch, xlim=xlim)
                start = i * at_once
                for (j in (start+1):ncol(x@data)){
                        lines(x@seconds, x@data[,j], col = newcol[j-start])
                        highlight <- ranges[1,j]:ranges[2,j]
                        lines(x@seconds[highlight], x@data[highlight,j],col = newcol[j-start], lwd=3)

                }
                par(ask=F)
        }
)


#show
setMethod("show", signature(object = "Enzyme_Assay"),
        function(object){
                h <- dim(object@data)[1]
                l <- dim(object@data)[2]
                if (length(object@curated_slope) > 0){
                        status <- 'curated'
                }else{
                        status <- 'not curated'
                }
                cat('\nEnzyme_Assay object')
                cat('\n\tStatus:', status)
                cat('\n\tSamples:', l)
                cat('\n\tMeasurements:', h)
                cat('\n\tAutomatic max slopes:\n\t\t', object@auto_slope[1:5], '...')
                cat('\n\tAutomatic ranges:\n\t\t', object@auto_range[1,1:5], '...')
                cat('\n\t\t', object@auto_range[2,1:5], '...', '\n')
                if (length(object@curated_slope) > 0){
                        cat('\n\tCurated max slopes:\n\t\t', object@curated_slope[1:5], '...')
                        cat('\n\tCurated ranges:\n\t\t', object@curated_range[1,1:5], '...')
                        cat('\n\t\t', object@curated_range[2,1:5], '...', '\n')
                }
        }
)
#########
#datasets
#########


#' Example Enzyme Assay Output
#'
#' An example enzyme assay run with absorbance values.
#' Rows for different measurement time points and columns
#' for different wells. 
#' 
#' @format A data frame with 69 rows and 42 columns
#' @name enzyme_assay.data
NULL


#' Example Plate Reader Timepoints
#'
#' An example vector of time values produced by the plate
#' reader that corresponds the the measurements in 
#' enzyme_assay.time 
#' 
#' @format A data frame with 69 rows and 42 columns
#' @name enzyme_assay.time
NULL
