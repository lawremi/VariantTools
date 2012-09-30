### =========================================================================
### Various apply utilities
### -------------------------------------------------------------------------

safe_mclapply <- function(X, FUN, ...) {
  ans <- IRanges:::castList(mclapply(X, function(...) {
    try(FUN(...), silent = !getOption("verbose"))
  }, ...))
  cond <- batchCondition(ans)
  if (!is.null(cond))
    stop(cond)
  ans
}

safe_mcmapply <- function(FUN, ...) {
  ans <- IRanges:::castList(mcmapply(function(...) {
    try(FUN(...), silent = !getOption("verbose"))
  }, ..., SIMPLIFY = FALSE))
  cond <- batchCondition(ans)
  if (!is.null(cond)) {
    stop(cond)
  }
  ans
}

batchCondition <- function(x) {
  errors <- sapply(x, is, "try-error")
  if (any(errors)) {
    mesg <- paste(sum(errors), "job(s) encountered an error.")
    error_mesg <- unlist(x[errors], use.names = FALSE)
    unique_mesg <- unique(error_mesg)
    if (length(unique_mesg) == 1L)
      mesg <- paste(mesg, "The error was the same for all nodes:", unique_mesg)
    mesg <- paste(mesg,
                  "Call 'summary' or 'detail' on this condition object",
                  "for a report.")
    cond <- simpleError(mesg)
    cond$detail <- DataFrame(job = which(errors), error = error_mesg)
    class(cond) <- c("batchCondition", class(cond))
    cond
  } else NULL
}

setOldClass("batchCondition")

setMethod("summary", "batchCondition", function(object) {
  rename(as(xtabs(~ error, object$detail), "DataFrame"), Freq = "count")
})

setMethod("detail", "batchCondition", function(x) {
  x$detail
})

setMethod("show", "batchCondition", function(object) {
  cat(conditionMessage(object))
})

## x <- 1:10
## testFun <- function(x) {
##   if (x < 3)
##     stop("x must be >= 3")
##   if (x > 7)
##     stop("x must be <= 7")
##   x
## }
## safe_mclapply(x, testFun)
## tryCatch(safe_mclapply(x, testFun),
##          batchCondition = summary)
## tryCatch(safe_mcmapply(testFun, x),
##          batchCondition = summary)

applyByChromosome <- function(X, FUN, ...) {
  if (!is(X, "Seqinfo"))
    X <- seqinfo(X)
  gr <- as(X, "GenomicRanges")
  ind <- seq_len(length(gr))
  safe_mclapply(ind, function(i, ...) {
    FUN(gr[i], ...)
  }, mc.preschedule = FALSE, ...)
}
