#' Class whippetDataSet
#'
#' Class \code{whippetDataSet} contains information read from whippet output files
#'
#' @name whippetDataSet-class
#' @rdname whippetDataSet-class
#' @exportClass whippetDataSet
#' @imports methods
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
#' @imports methods
#' @param whippetDataSet whippetDataSet generated from \code{readWhippetDataSet()}
setGeneric("diffSplicingResults",
           def=function(whippetDataSet)
           {
               standardGeneric("diffSplicingResults")
           }
)

#' @rdname diffSplicingResults-methods
#' @return differential splicing results data.frame
#' (originally from a whippet .diff file)
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
#' @imports methods
#' @param whippetDataSet whippetDataSet generated from \code{readWhippetDataSet()}
setGeneric("readCounts",
           def=function(whippetDataSet)
           {
               standardGeneric("readCounts")
           }
)
#' @rdname readCounts-methods
#' @return whippet read count data.frame
#' (originally from a whippet .psi file)
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
#' @imports methods
#' @param whippetDataSet whippetDataSet generated from \code{readWhippetDataSet()}
setGeneric("junctions",
           def=function(whippetDataSet)
           {
               standardGeneric("junctions")
           }
)
#' @rdname junctions-methods
#' @return junctions GRanges object
#' (originally from a whippet .jnc file)
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
#' @return whippet event coordinates GRanges
setMethod("coordinates", signature="whippetDataSet",
          definition=function(whippetDataSet)
          {
              return(whippetDataSet@coordinates)
          }
)
