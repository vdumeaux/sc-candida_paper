if(R.version$major >= 3){
  require(parallel)
} else {
  require(multicore)
}

###
### file: bresat.R
### author: Ali Tofigh
###
### Contains functions related to Robert Lesurf's BreSAT ideas and
### signature database sigdb.


### sig.ranksum()
###
### Robert Lesurf's main BreSAT function that ranks the samples given a
### signature. In this implementation, patients with the same ranksum
### are assigned their "average" rank, i.e., they all receive the same
### rank value, which might not be an integer. One benefit of this
### strategy is that exchanging the up- and down-genes reverses the
### ranks, and ordering the patients in decreasing order (with a stable
### sort algorithm) will be the same as the order of the patients in
### increasing order before the exchange. This consistent behavior in
### turn makes sig.correlation behave consistently when computing
### correlations between patient ranks with respect to different gene
### signatures.
###
### Note that 'up' and 'dn' are used to index into exprdata. Hence, these can
### be either numerical indices or row names. NAs in exprdata are not
### supported.
###
### Arguments:
###     exprdata        the gene expression data, must be a matrix
###     up              indices of up-regulated genes
###     dn              indices of down-regulated genes
###     ns              indices of genes whose directions have not been specified.
###                     These genes will be clustered into two parts based on
###                     their correlation to each other and up and dn indices
###                     will be assigned by the function. It is currently not
###                     possible to mix 'up'/'dn' with 'ns'.
###     full.return     returns a struct when true, otherwise only the rank
###
### Returns:
###     if 'full.return' is FALSE, then the computed rank is
###     returned. Otherwise, a list is returned with the following members.
###
###     rank            rank of each patient, higher = better
###     dat             exprdata for the signature genes only and ordered both by
###                     patient and gene.
###     up.dn           vector with nrow(dat) values. -1, indicates down-gene
###                     and 1 indicates up-gene
###     pat.order       patient ordering, with best patient first.
###                     exprdata[, pat.order] orders the columns according to the
###                     strength of the signature.
###     gene.order      exprdata[gene.order, ] sorts rows of exprdata so that all
###                     up-genes are first and all down genes last.  Also, the
###                     most correlated up-gene is at top and most correlated
###                     down-gene is at the bottom.
sig.ranksum <- function(exprdata, up=NULL, dn=NULL, ns=NULL, full.return = FALSE, ranks.matrix=FALSE)
{
  if(ranks.matrix && full.return)
    stop("Cannot return full results of bresat ranksum with ranks as input")
  if(ranks.matrix && length(ns) > 0)
    stop("Cannot compute bresat ranksum with ranks as input on a signature without identified up and dn genes")
  
  if (length(dim(exprdata)) != 2)
        stop("'exprdata' must be a matrix")
    if (length(up) == 0 && length(dn) == 0 && length(ns) == 0)
        stop("no indices were specified")
    if (length(ns) > 0 && (length(up) > 0 || length(dn) > 0))
        stop("both directional and non-directional indices were specified")
    if (is.logical(up))
        up <- which(up)
    if (is.logical(dn))
        dn <- which(dn)
    if (is.logical(ns))
        ns <- which(ns)

    if (ncol(exprdata) < 2) {
      ret <- list()
      up <- c(up,ns)
      ret$rank <- 1
      ret$pat.order <- 1
      ret$gene.order <- c(up,dn)
      ret$dat <- exprdata[ret$gene.order,, drop=FALSE]
      ret$up.dn <- c(rep(1, length(up)), rep(-1, length(dn)))
      return(ret)
    }

    if (length(ns) > 0)
    {
        ## identify rows in exprdata that have zero variance
        library(matrixStats)
        tmp <- rowSds(exprdata[ns,,drop=F]) == 0 ## faster
        #tmp <- sapply(ns, function(idx) {sd(exprdata[idx, ], na.rm=TRUE)}) == 0
        zero.sd.idx <- ns[tmp]
        ns <- ns[!tmp]

        if (length(ns) == 1)
        {
            up <- ns
        }
        else if (length(ns) == 2)
        {
            if (cor(exprdata[ns[1], ], exprdata[ns[2], ], use="pairwise") < 0)
            {
                up <- ns[1]
                dn <- ns[2]
            }
            else
            {
                up <- ns
            }
        }
        else if (length(ns) > 2)
        {
            library(cluster)
            diss <- 1-cor(t(exprdata[ns, , drop=FALSE]), use="pairwise")
            diss[which(is.na(diss))] <- 1
            diss[diss >= 1] <- diss[diss >= 1] + 1

            clustering <- pam(diss, k=2, diss=TRUE, cluster.only=TRUE)
            up.cluster <- which.max(table(clustering))
            up <- ns[which(clustering == up.cluster)]
            dn <- ns[which(clustering != up.cluster)]
            up.dn.cor <- cor(t(exprdata[up, , drop=FALSE]), t(exprdata[dn, , drop=FALSE]), use="pairwise")
            if (sum(up.dn.cor < 0,na.rm=T) < length(up) * length(dn) / 2)
            {
                up <- ns
                dn <- NULL
            }
        }

        if (length(zero.sd.idx) > 0)
        {
            up <- c(up, zero.sd.idx)
        }
    }

    ranksum <- double(ncol(exprdata))
    col.counts <- rep(0, ncol(exprdata))
    if (length(up) != 0)
    {
        dat <- exprdata[up, , drop = FALSE]
        if(ranks.matrix)
          ranksum <- colSums(dat, na.rm=TRUE)
        else
        ranksum <- rowSums(apply(dat, 1, function(x) {.Internal(rank(x, length(x), "average"))}))
        col.counts <- colSums(!is.na(dat))
    }
    if (length(dn) != 0)
    {
        dat <- exprdata[dn, , drop = FALSE]
        if(ranks.matrix)
          ranksum <- ranksum + colSums(ncol(exprdata) - dat + 1, na.rm=TRUE)
        else
        ranksum <- ranksum + rowSums(ncol(exprdata) - apply(dat, 1, function(x) {.Internal(rank(x, length(x), "average"))}) + 1)
        col.counts <- col.counts + colSums(!is.na(dat))
    }

  ranksum <- ranksum / col.counts  
  rank <- .Internal(rank(ranksum, length(ranksum),"average"))

    if (full.return == FALSE)
        return(rank)

    if (length(up) == 0)
        up <- NULL
    if (length(dn) == 0)
        dn <- NULL

    ## computes correlation of a single row of the expression data with the
    ## patient ordering. Used to obtain ret$gene.order.
    gene.cor <- function(gene.idx, is.up, exprdata, pat.order)
    {
        gene.expr <- exprdata[gene.idx, pat.order]
        if (is.up == TRUE)
            gene.expr <- -gene.expr
        cor(rank(gene.expr), 1:ncol(exprdata))
    }

    ret <- list()
    ret$rank <- rank
    ## best first, patient nr pat.order[1] is the best patient for this sig so
    ## that dat[ , pat.order] sorts the patients correctly with nr 1 = best.
    ret$pat.order <- order(-rank)
    ## Also, gene.order orders the genes so that the up-genes come first,
    ## followed by the down-genes. The "best" up-gene first and the "best"
    ## down-gene last.
    up.cor <- sapply(up, gene.cor, TRUE, exprdata, ret$pat.order)
    dn.cor <- sapply(dn, gene.cor, FALSE, exprdata, ret$pat.order)
    ret$gene.order <- c(up[rev(order(up.cor))], dn[order(dn.cor)])
    ret$up <- up
    ret$dn <- dn
    ret$up.cor<-up.cor
    ret$dn.cor<-dn.cor
    ret$dat <- exprdata[ret$gene.order, ret$pat.order, drop=FALSE]
    ret$up.dn <- c(rep(1, length(up)), rep(-1, length(dn)))
    ret$ranksum <- ranksum

    return(ret)
}