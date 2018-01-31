#' Class whippetDataSet
#'
#' Class \code{whippetDataSet} contains information read from whippet output files
#'
#' @name whippetDataSet-class
#' @rdname whippetDataSet-class
#' @exportClass whippetDataSet
setClass("whippetDataSet", slots=list(coordinates="GRanges",
                                      diffSplicingResults="data.frame",
                                      comparisons="character",
                                      junctions="GRanges",
                                      readCounts="data.frame",
                                      filePath="character"))

#' Method diffSplicingResults
#' @name diffSplicingResults
#' @rdname diffSplicingResults-methods
#' @exportMethod diffSplicingResults
#' @param whippetDataSet whippetDataSet generated from \code{readWhippetDataSet()}
setGeneric("diffSplicingResults",
           def=function(whippetDataSet)
           {
               standardGeneric("diffSplicingResults")
           }
)

#' @rdname diffSplicingResults-methods
setMethod("diffSplicingResults", signature="whippetDataSet",
          definition=function(whippetDataSet)
          {
              return(whippetDataSet@diffSplicingResults)
          }
)

#' Method readCounts
#' @name readCounts
#' @rdname readCounts-methods
#' @exportMethod readCounts
#' @param whippetDataSet whippetDataSet generated from \code{readWhippetDataSet()}
setGeneric("readCounts",
           def=function(whippetDataSet)
           {
               standardGeneric("readCounts")
           }
)
#' @rdname readCounts-methods
setMethod("readCounts", signature="whippetDataSet",
          definition=function(whippetDataSet)
          {
              return(whippetDataSet@readCounts)
          }
)

#' Method junctions
#' @name junctions
#' @rdname junctions-methods
#' @exportMethod junctions
#' @param whippetDataSet whippetDataSet generated from \code{readWhippetDataSet()}
setGeneric("junctions",
           def=function(whippetDataSet)
           {
               standardGeneric("junctions")
           }
)
#' @rdname junctions-methods
setMethod("junctions", signature="whippetDataSet",
          definition=function(whippetDataSet)
          {
              return(whippetDataSet@junctions)
          }
)

#' Method coordinates
#' @name coordinates
#' @rdname coordinates-methods
#' @exportMethod coordinates
#' @param whippetDataSet whippetDataSet generated from \code{readWhippetDataSet()}
setGeneric("coordinates",
           def=function(whippetDataSet)
           {
               standardGeneric("coordinates")
           }
)
#' @rdname coordinates-methods
setMethod("coordinates", signature="whippetDataSet",
          definition=function(whippetDataSet)
          {
              return(whippetDataSet@coordinates)
          }
)
