## TODO: use dispatch on special file classes
file_ext_sans_gz <- function(x) {
  x <- sub("[.]b?gz$", "", x)
  file_ext(x)
}

loadVariants <- function(x, ...) {
  if (file_ext_sans_gz(x) == "vcf") {
    vcf <- readVcf(x, ...)
    geno(vcf)$DP <- NULL
    as(vcf, "VRanges")
  }
  else get(load(x))
}

calculateConcordanceMatrix <- function(variantFiles, ...) {
  if (length(variantFiles) < 2L) {
    stop("There must at least two variant sets.")
  }
  if (is.character(variantFiles)) {
      if (!all(file.exists(variantFiles))) {
          stop("Some of the files could not be found.")
      }
      nms <- variantFiles
  } else {
      if (!all(vapply(variantFiles, is, class2="VRanges", logical(1L))))
          stop("'variantFiles' must be a vector of file names or ",
               "a list of VRanges objects")
      nms <- names(variantFiles)
  }

  ## initialize vcmat, pairwise variance concordant score between dirs
  n <- length(variantFiles)
  vcmat <- matrix(NA, nrow=n, ncol=n)
  rownames(vcmat) <- colnames(vcmat) <- nms

  ## compute vcmat
  for (i in 1:(n - 1)) {
    avar <- if (is(variantFiles[[i]], "VRanges"))
                x
            else try(loadVariants(variantFiles[i], ...), silent=TRUE)
    if (class(avar) == "try-error") {
      stop("error: cannot load '", nms[i], "': ", avar)
    }
  
    for (j in (i+1):n) {
      bvar <- if (is(variantFiles[[j]], "VRanges"))
                  x
              else try(loadVariants(variantFiles[j], ...), silent=TRUE)
      if (class(bvar) == "try-error") {
        stop("error: cannot load '", nms[j], "': ", bvar)
      }
    
      ## compute variant concordance
      vc <- calculateVariantConcordance(avar, bvar)

      ## output results
      vcmat[i, i] <- 1
      vcmat[j, j] <- 1
      vcmat[i, j] <- vc[1]$fraction
      vcmat[j, i] <- vc[1]$fraction

      if (getOption("verbose"))
        message(paste("variantConcordance:",
                      nms[i],
                      nms[j],
                      vc$fraction, vc$denominator, sep="\t"))
    }
  }
  return(vcmat)  
}

callVariantConcordance <- function(concordanceMatrix,
                                   threshold) {
  if (!installed("RBGL"))
    stop("The 'RBGL' package is required for finding the concordant cliques")
  
  concordantCliques <- .getConcordantCliques(concordanceMatrix,
                                             threshold)
  ## test if connected components are cliques of
  ##at least 2 vertices

  n <- sum(elementNROWS(concordantCliques))
  vstate <- rep("non-concordant", n)
  names(vstate) <- unlist(concordantCliques)

  ##if there are multiple cliques of size 2 or greater, the entire set
  ##is undecidable
  numMultinodeCliques <- sum(elementNROWS(concordantCliques) >= 2)
  if (numMultinodeCliques >= 2L) {
    vstate <- rep("undecidable", n)
  } else if (numMultinodeCliques == 1L) {
    ##if in the clique with at least two elements, is concordant
    multinodeCliqueElements <-
      concordantCliques[[which(elementNROWS(concordantCliques) >= 2)]]
    vstate[multinodeCliqueElements] <- "concordant"
  }
  return(vstate)
}

##################
##HELPER FUNCTIONS
##################

.getConcordantCliques <- function(concordanceMatrix, threshold) {
  adj <- concordanceMatrix > threshold
  if (is.null(colnames(adj)))
    colnames(adj) <- paste0("id", seq_len(ncol(adj)))
  V <- colnames(adj)
  edL <- vector("list", length=length(V))
  names(edL) <- V
  for (i in seq_along(V)) {
    indicesOfConnectedNodes <- which(adj[, i] == TRUE)
    edL[[i]] <- list(edges=indicesOfConnectedNodes)
  }
  gR <- graph::graphNEL(nodes=V, edgeL=edL)
  RBGL::maxClique(gR)$maxCliques
}
