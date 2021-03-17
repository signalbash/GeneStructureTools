#' Class whippetDataSet
#'
#' Class \code{whippetDataSet} contains information read from whippet output files
#'
#' @name whippetDataSet-class
#' @rdname whippetDataSet-class
#' @exportClass whippetDataSet
#' @import methods
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
#' @import methods
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
#' @family whippet data processing
#' @examples
#' whippetFiles <- system.file("extdata","whippet/",
#' package = "GeneStructureTools")
#' wds <- readWhippetDataSet(whippetFiles)
#'
#' diffSplicingResults <- diffSplicingResults(wds)
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
#' @import methods
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
#' @family whippet data processing
#' @examples
#' whippetFiles <- system.file("extdata","whippet/",
#' package = "GeneStructureTools")
#' wds <- readWhippetDataSet(whippetFiles)
#'
#' readCounts <- readCounts(wds)
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
#' @import methods
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
#' @family whippet data processing
#' @examples
#' whippetFiles <- system.file("extdata","whippet/",
#' package = "GeneStructureTools")
#' wds <- readWhippetDataSet(whippetFiles)
#'
#' junctions <- junctions(wds)
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
#' @return whippet splicing event coordinates as a GRanges object
#' @family whippet data processing
#' @examples
#' whippetFiles <- system.file("extdata","whippet/",
#' package = "GeneStructureTools")
#' wds <- readWhippetDataSet(whippetFiles)
#'
#' coordinates <- coordinates(wds)
setMethod("coordinates", signature="whippetDataSet",
          definition=function(whippetDataSet)
          {
              return(whippetDataSet@coordinates)
          }
)



##' Class rmatsDataSet
#'
#' Class \code{rmatsDataSet} contains information read from rmats output files
#'
#' @name rmatsDataSet-class
#' @rdname rmatsDataSet-class
#' @exportClass rmatsDataSet
#' @import methods
setClass("rmatsDataSet", slots=list(SE="data.frame",
                                    MXE="data.frame",
                                    RI="data.frame",
                                    A3SS="data.frame",
                                    A5SS="data.frame",
                                    filePath="character"))

#' Method extractEvent
#' @name extractEvent
#' @rdname extractEvent-methods
#' @exportMethod extractEvent
#' @import methods
#' @param rmatsDataSet rmatsDataSet generated from \code{readRmatsDataSet()}
#' @param eventType specific event type to extract results for. Must be SE/MXE/RI/A5SS/A3SS.
setGeneric("extractEvent",
           def=function(rmatsDataSet, eventType)
           {
               standardGeneric("extractEvent")
           }
)

#' @rdname extractEvent-methods
#' @return differential splicing results data.frame
#' (originally from a whippet .diff file)
#' @family rmats data processing
setMethod("extractEvent", signature="rmatsDataSet",
          definition=function(rmatsDataSet, eventType)
          {
              return(slot(rmatsDataSet, eventType))
          }
)

#' Class irfDataSet
#'
#' Class \code{irfDataSet} contains information read from irf output files
#'
#' @name irfDataSet-class
#' @rdname irfDataSet-class
#' @exportClass irfDataSet
#' @import methods
setClass("irfDataSet", slots=list(coordinates="GRanges",
                                  IRFresults="data.frame",
                                  filePath="character"))

#' Method irfResults
#' @name irfResults
#' @rdname irfResults-methods
#' @exportMethod irfResults
#' @import methods
#' @param irfDataSet irfDataSet generated from \code{readIRFDataSet()}
setGeneric("irfResults",
           def=function(irfDataSet)
           {
               standardGeneric("irfResults")
           }
)

#' @rdname irfResults-methods
#' @return differential splicing results data.frame
#' (originally from a whippet .diff file)
#' @family irf data processing
setMethod("irfResults", signature="irfDataSet",
          definition=function(irfDataSet)
          {
              return(slot(irfDataSet, "IRFresults"))
          }
)
