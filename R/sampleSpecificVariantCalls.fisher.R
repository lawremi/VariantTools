##' sampleSpecificVariantCalls.fisher estimates sample specific variant calls by comparing with the control sample
##' case: a GRanges object keeps the variant calss of the case sample
##' control: a GRanges object keeps the variant calss of the control sample
##' control.cov: the coverage of the control sample (rle.list)
##' min.cov.case: minimum coverage of the case sample
##' min.cov.ctrl: minimum coverage of the control sample
##' max.freq.ctrl: maximum variant frequency of the control sample corresponding to the variant of the case sample
##' max.fisher.pval: maximum pvalue of fisher.test between the case and control samples
##' fdr.adjust: whether to do fdr adjustment of the fisher.test pvalues
## Returned value: a GRanges object (filtered "case") with additional control sample information and fisher.test p.values
sampleSpecificVariantCalls.fisher <- function(case, control, control.cov, min.cov.case=10, min.cov.ctrl=10, max.freq.ctrl=0.03, max.fisher.pval=0.05, fdr.adjust=FALSE) {
	
	## match the variant and add missing control counts
	case <- annotateWithControlCounts(case, control, control.cov)
	
	## check whether there is high.quality count information, if not just use regular read counts
	if (is.null(case$high.quality) || is.null(control$high.quality)) {
		warning('high.quality count information is not available, regular count information will be used!\n')
		count.total <- case$count.total
		count <- case$count
		count.ref <- case$count.ref
	} else {
		count.total <- case$high.quality.total
		count <- case$high.quality
		count.ref <- case$high.quality.ref

		## Because the high.quality counts of indels are NAs, regular counts are used insteadh
		ind.indel <- which(nchar(case$alt) != 1 | nchar(case$ref) != 1)
		count.total[ind.indel] <- case$count.total[ind.indel]
		count[ind.indel] <- case$count[ind.indel]
		count.ref[ind.indel] <- case$count.ref[ind.indel]
	}

  ## filter on number of control and case coverage and frequency in control 
  freq.ctrl <- case$control.count / (case$control.count.total + 0.01)  ## adding 0.01 to avoid NAs
	filter.status <- (count.total >= min.cov.case) & 
									(case$control.count.total >= min.cov.ctrl) & 
									(freq.ctrl < max.freq.ctrl)
  selInd <- which(filter.status) ## remove NAs
	case <- case[selInd]

 	## filter based on fisher p-value
  mm <- cbind(values(case)$control.count,  count[selInd],  
				values(case)$control.count.total - values(case)$control.count,  count.ref[selInd])
  fisher.pval <- apply(mm, 1, function(x) fisher.test(matrix(x, nrow=2))$p.value)
	case$pvalue <- fisher.pval
	case$pvalue.adjusted <- p.adjust(fisher.pval, method='fdr')
	if (fdr.adjust) {
	  case <- case[which(case$pvalue.adjusted <= max.fisher.pval)]
	} else {
	  case <- case[which(fisher.pval <= max.fisher.pval)]
	}

  ## return filtered case
  return (case)
}
