#' Check that a list of tables is in the correct format as a metabarlist
#'
#' Test that a list of tables containing information on MOTU abundances, MOTU characteristics, pcr characteristics, and sample characteristics form a congruent \code{\link{metabarlist}}
#'
#' @param metabarlist a \code{metabarlist} object
#' @param groups a vector containing the identifier of replicate. The vector must has the same length of the table `PCRs` from a \code{\link{metabarlist}} object. Default = metabarlist$pcrs$sample_id
#' @param FUN a function returning a distance matrix. The distance matrix should be a object `dist` which has the same length of table `reads`.
#'
#' @name aggregate_reads
#'
#' @return TRUE or throws a stop
#'
#' @details
#'#'
#' @examples
#'
#' data(soil_euk)
#'
#' aggregate_reads(soil_euk)
#'
#' @author Frédéric Boyer
#' @export aggregate_reads


aggregate_reads <- function(metabarlist, groups=metabarlist$pcrs$sample_id, FUN=sum) {
  #TODO check that rows order in reads are in the same order as PCRs
  stopifnot(all(rownames(metabarlist$reads) == rownames(metabarlist$pcrs)))

  gu <- unique(metabarlist$pcrs$sample_id)
  newReads <- do.call(rbind,
    lapply(gu, function(group) {
      apply(metabarlist$reads[groups==group, , drop=F], MARGIN=2, FUN=FUN)
    })
  )
  rownames(newReads) <- gu
  newPCRs <- aggregate(metabarlist$pcrs, by=list(groups), FUN=function(x) {
    v<-unique(x);
    ifelse(length(v)==1, v, NA)
  })
  print(newPCRs)
  v <- apply(newPCRs, MARGIN=2, FUN=function(v) any(is.na(v)))
  newPCRs <- newPCRs[, !v, drop=F]
  rownames(newPCRs) <- newPCRs[,1]
  newPCRs <- newPCRs[,-1]
  newPCRs$control_type <- NA
  newReads <- newReads[rownames(newPCRs),]
  nm <- list(reads = newReads, pcrs = newPCRs, motus = metabarlist$motus, samples = metabarlist$samples)
  class(nm) <- "metabarlist"
  check_metabarlist(nm)
  return(nm)
}
