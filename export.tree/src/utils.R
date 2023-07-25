prep.pseudobulk <- function(all.metadata=all.metadata, cluster.names = "leiden_scVI_2"){
  set.seed(321)

  
  # Named vector of cluster names
  kids <- purrr::set_names(levels(as.factor(all.metadata[[cluster.names]])))
  kids
  
  # Total number of clusters
  nk <- length(kids)
  
  # Named vector of sample names
  sids <- purrr::set_names(levels(all.metadata$orig.ident))
  
  # Total number of samples 
  ns <- length(sids)
 
  # Generate sample level metadata
  all.metadata$others <- ifelse(all.metadata$leiden_scVI_2=="others", TRUE, FALSE)
  
  
  ## Turn named vector into a numeric vector of number of cells per sample
  n_samples <- colSums(table(all.metadata$batch, all.metadata$drug_day) > 0)
  names(n_samples) <- unique(all.metadata$drug_day)
  
  single <- names(n_samples)[n_samples==1]
  
  all.metadata$n_samples <- n_samples[all.metadata$drug_day]
  
  all.metadata$orig.ident.split <- as.character(all.metadata$orig.ident)

  idx <- idx.1 <- idx.2 <- list()
  for (i in seq_along(single)){
    
    idx[[i]] <- which(all.metadata$drug_day==single[i])
    idx.1[[i]] <- sample(idx[[i]], length(idx[[i]])/2)
    idx.2[[i]] <- idx[[i]][!idx[[i]] %in% idx.1[[i]]]
    
    
    all.metadata$orig.ident.split[idx.1[[i]]] <- paste(as.character(all.metadata$orig.ident)[idx.1[[i]]], "1", sep = "_")
    all.metadata$orig.ident.split[idx.2[[i]]] <- paste(all.metadata$orig.ident[idx.2[[i]]], "2", sep = "_")
  }
  
  all.metadata$orig.ident.split <- as.factor(all.metadata$orig.ident.split)
  # Named vector of sample names
  sids <- purrr::set_names(levels(all.metadata$orig.ident.split))
  
  
  ## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
  m <- match(sids, all.metadata$orig.ident.split)
  n_cells <- as.numeric(table(all.metadata$orig.ident.split))
  
  
  ## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
  ei <- data.frame(all.metadata[m, ], 
                   n_cells, row.names = NULL) %>% 
    dplyr::select(-all_of(cluster.names))

  # Aggregate the counts per sample_id and cluster_id
  
  # Subset metadata to only include the cluster and sample IDs to aggregate across
  
  all.metadata$orig.ident2 <- gsub("_", ".", as.character(all.metadata$orig.ident.split))
  groups <- all.metadata[, c(cluster.names, "orig.ident2")]
  
  # Aggregate across cluster-sample groups
  pb <- Matrix.utils::aggregate.Matrix(round(raw.counts), 
                                       groupings = groups, fun = "sum") 
  
  # Not every cluster is present in all samples; create a vector that represents how to split samples
  splitf <- sapply(stringr::str_split(rownames(pb), 
                                      pattern = "_",  
                                      n = 2), 
                   `[`, 1)



  # Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
  pb <- split.data.frame(pb, factor(splitf)) %>%
    lapply(function(u) 
      set_colnames(t(u),  sapply(stringr::str_split(rownames(u), pattern = "_", n = 2), `[`, 2)))
  
  
  # prep. data.frame for plotting
  get_sample_ids <- function(x){
    pb[[x]] %>%
      colnames()
  }
  
  de_samples <- map(1:length(kids), get_sample_ids) %>%
    unlist()

  # Get cluster IDs for each of the samples
  
  samples_list <- map(1:length(kids), get_sample_ids)
  
  get_cluster_ids <- function(x){
    rep(names(pb)[x], 
        each = length(samples_list[[x]]))
  }
  
  de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
    unlist()
  
  # Create a data frame with the sample IDs, cluster IDs and condition
  
  gg_df <- data.frame(cluster_id = de_cluster_ids,
                      sample_id = de_samples)
  
  ei$sample_id <- gsub("_", ".",  ei$orig.ident.split)
  ei$group_id <- ei$drug_day
  
  gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id", "batch", "drug", "others")])
  metadata <- gg_df %>%
    dplyr::select(cluster_id, sample_id, group_id, batch, drug, others) 
  

  # Generate vector of cluster IDs
  metadata$cluster_id <- as.factor(metadata$cluster_id)

  res <- list(pb=pb, metadata = metadata)
  return(res)
}



myreadGAF = function(filename,evidence=NULL, aspect=c("P", "F", "C")){
  
  ## the column IDs of interest according to GAF 1.0 and 2.0
  gene.id.col = 2
  symbol.col = 3
  go.id.col = 5
  evidence.col = 7
  aspect.col = 9
  name.col = 10
  taxid.col = 13
  
  ## validity of parameters
  # if( !( is.null(evidence) | is.character(evidence) ) )
  #   stop("evidence must be NULL or a character vector.")
  # 
  if(!all(aspect %in% c("P", "F", "C")))
    stop("aspect must be a character vector with all entries in c(\"P\", \"F\", \"C\")")
  
  ## reading the file
  goa = read.delim(gzfile(filename), na.strings = "", header=FALSE, comment.char = "!", sep="\t")
  
  ## remove associations with NOT qualifiers
  goa <- goa[!is.na(goa$V4) > 0,]
  excl.not <- grep("NOT", goa$V4)
  goa <- goa[-excl.not,]
  
  goa = 
    data.frame ( 
      taxid = gsub("taxon:", "", goa[, taxid.col]),
      GOID = goa[,go.id.col],
      gene_symbol = goa[, gene.id.col],
      evidence = goa[, evidence.col],
      gene_id = goa[, symbol.col],
      name = goa[, name.col],
      aspect.code = goa[,aspect.col]
    )
  

  
  goa <- goa[goa$aspect.code %in% aspect, ]
  return(goa)}


rep.row <- function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

rep.col <- function(x,n){
  matrix(rep(x, n),ncol=n, byrow = FALSE)
}

grep.size <- function(x){
  frac <- unlist(lapply(strsplit(x, "\\("), "[", 2))
  num.frac <- unlist(lapply(strsplit(frac, "/"), "[", 1))
  num.frac <- as.numeric(num.frac)
  return(num.frac)
}


go.compare.dotplot <- function(sel){
  plot.dat <- list()
  for (group in c("day3", "both", "day6")){
    plot.dat1 <- tb.list[[group]] %>%
      dplyr::select(GO.cluster, IC, GO.ID, term, 
                    size.day3, paste(sel, "day3..log10_pvalue", sep="_"))
    colnames(plot.dat1)[5:6] <- c("size", "-log10(pvalue)")
    plot.dat1$day <- "day3"
    
    plot.dat2 <- tb.list[[group]] %>%
      dplyr::select(GO.cluster, IC, GO.ID, term,
                    size.day6, paste(sel, "day6..log10_pvalue", sep="_"))
    plot.dat2$day <- "day6"
    colnames(plot.dat2)[5:6] <- c("size", "-log10(pvalue)")
    
    plot.dat[[group]] <- rbind(plot.dat1, plot.dat2)
    plot.dat[[group]]$sigday <- group
  }
  
  plot.dat.all <- rbind(plot.dat$day3, plot.dat$both, plot.dat$day6)
  
  plot.dat.all$sigday <- factor(plot.dat.all$sigday, levels = c("day3", "both", "day6"), ordered = TRUE)
  level.order <- unique(plot.dat.all$term[order(plot.dat.all$GO.cluster)])
  plot.dat.all$term <- factor(plot.dat.all$term, levels = level.order, ordered = TRUE)
  plot.dat.all <- plot.dat.all[order(plot.dat.all$GO.cluster),]
  
  plot.dat.all$col <- as.character(lacroix_palette("PassionFruit", type = "continuous", n=16))[plot.dat.all$GO.cluster]
  
  return(plot.dat.all)
}