make_splici_txome <- function(gtf_path,
                              genome_path,
                              read_length,
                              flank_trim_length = 5,
                              output_dir,
                              extra_spliced=NULL,
                              extra_unspliced=NULL,
                              dedup_seqs=FALSE) {
  # if you get some error from .get_cds_IDX, please try to rerun the code again
  # read length is the scRNA-seq read length
  # flank trim length is used to avoid marginal case when dealing with junctional reads
  # assumes the following packages have been imported
  # eisaR, Biostrings, BSgenome, dplyr, stringr
  
  ########################################################################################################
  # Preprocessing
  ########################################################################################################
  
  suppressPackageStartupMessages({
    library(eisaR)
    library(Biostrings)
    library(BSgenome)
    library(stringr)
    library(GenomicFeatures)
  })
  
  if (!dir.exists(output_dir)) {
    dir.create(file.path(output_dir),recursive = TRUE, showWarnings = FALSE)
  }
  # make sure flank_length makes sense
  flank_length = read_length - flank_trim_length
  if (flank_length < 0 ){
    stop("flank trim length is larger than read length!")
  }
  # make sure gtf file exists
  if (!file.exists(gtf_path)) {
    stop("The following file does not exist: \n", gtf_path)
  }
  
  # make sure fasta file exists
  if (!file.exists(genome_path)) {
    stop("The following file does not exist: \n", genome_path)
  }
  
  # output file names
  file_name_prefix = paste0("transcriptome_splici_fl", flank_length)
  out_fa <- file.path(output_dir, paste0(file_name_prefix, ".fa"))
  out_t2g <- file.path(output_dir, paste0(file_name_prefix, "_t2g.tsv"))
  out_t2g3col <- file.path(output_dir, paste0(file_name_prefix, "_t2g_3col.tsv"))
  
  
  
  #########################################################################################################
  # Process gtf to get spliced and introns
  #########################################################################################################
  message("============processing gtf to get spliced and introns============")
  # fl is the flank length, here we set it to
  # the read length - 5
  grl <- suppressWarnings(getFeatureRanges(
    gtf = file.path(gtf_path),
    featureType = c("spliced", "intron"),
    intronType = "separate",
    flankLength = flank_length,
    joinOverlappingIntrons = TRUE,
    verbose = TRUE
  ))
  
  #########################################################################################################
  # Get spliced related stuffs
  #########################################################################################################
  
  # spliced ranges has no dash in it
  spliced_grl = grl[str_detect(names(grl), "-I", negate = TRUE)]
  
  #########################################################################################################
  # Get reduced introns
  #########################################################################################################
  
  # identify all introns and convert to GRanges
  intron_gr = unlist(grl[str_detect(names(grl), "-I")])
  # group introns by gene, then collapse ovelaping ranges!
  intron_grl = reduce(split(intron_gr, intron_gr$gene_id))
  
  # clean txp names and gene names
  intron_gr <- BiocGenerics::unlist(intron_grl)
  intron_gr$exon_rank <- 1L
  intron_gr$transcript_id <- word(names(intron_gr), 1, sep = '-I')
  intron_gr$gene_id <- intron_gr$transcript_id
  intron_gr$type <- "exon"
  intron_gr$transcript_id <- make.unique(paste0(intron_gr$transcript_id, "-I"), sep = '')
  intron_gr$gene_id <- paste0(intron_gr$gene_id, "-I")
  intron_gr$exon_id <- intron_gr$transcript_id
  names(intron_gr) <- NULL
  mcols(intron_gr) <-
    S4Vectors::mcols(intron_gr)[, c("exon_id", "exon_rank",
                                    "transcript_id", "gene_id", "type")]
  # remake intron GRangesList
  intron_grl <- BiocGenerics::relist(intron_gr, lapply(
    structure(seq_along(intron_gr),
              names = intron_gr$transcript_id), function(i) i))
  
  
  #########################################################################################################
  # extract sequences from genome
  #########################################################################################################
  
  message("============extracting spliced and intron sequences from genome============")
  
  # load the genome sequence
  x <- Biostrings::readDNAStringSet(file.path(genome_path))
  # get the first word as the name
  names(x) <- word(names(x), 1)
  
  
  grl = c(spliced_grl, intron_grl)
  
  # make sure introns don't out of boundary
  seqlevels(grl) <- seqlevels(x)
  seqlengths(grl) <- suppressWarnings(seqlengths(x))
  grl <- GenomicRanges::trim(grl)
  
  seqs <- GenomicFeatures::extractTranscriptSeqs(
    x = x,
    transcripts = grl
  )
  
  # If having duplicated sequences, only keep one
  if(dedup_seqs) {
    seqs = unique(seqs)
    grl = grl[names(seqs)]
  }
  
  
  # save some space
  rm(x)
  #########################################################################################################
  # process final outputs
  #########################################################################################################
  message("Writing outputs...")
  
  df <- getTx2Gene(grl)
  write.table(df, out_t2g, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
  df <- df %>%
    dplyr::mutate(gene_id = word(gene_id, 2, sep = '_'),
                  status = ifelse(str_detect(transcript_id, '-I'), 'U', 'S'))
  
  writeXStringSet(seqs, out_fa, format = "fasta")
  write.table(df, file.path(output_dir, paste0(file_name_prefix, "_t2g_3col.tsv")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
  
  # optional: adding extra spliced and unspliced sequences from an fasta file
  if (!is.null(extra_spliced)) {
    if (!file.exists(extra_spliced)) {
      warning("provided extra_sequences file does not exist, will ignore it")
    } else {
      fa = file(extra_spliced, open="r")
      lns = readLines(fa)
      close(fa)
      for (ln in lns) {
        if (startsWith(ln, ">")) {
          # it is a header, write to t2g files and fasta file
          txp_name = gsub(">", "", ln)
          write.table(matrix(c(txp_name, txp_name), nrow = 1), file = out_t2g, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
          write.table(matrix(c(txp_name, txp_name, "S"), nrow = 1), file = out_t2g3col, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
          write.table(ln, file = out_fa, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
        } else {
          # if not a header, just write to fasta file
          write.table(ln, file = out_fa, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
        }
      }
    }
  }
  
  if (!is.null(extra_unspliced)) {
    if (!file.exists(extra_unspliced)) {
      warning("provided extra_sequences file does not exist, will ignore it")
    } else {
      fa = file(extra_unspliced, open="r")
      lns = readLines(fa)
      close(fa)
      for (ln in lns) {
        if (startsWith(ln, ">")) {
          # it is a header, write to t2g file and fasta file
          txp_name = gsub(">", "", ln)
          write.table(matrix(c(txp_name, paste0(txp_name, "-U")), nrow = 1), file = out_t2g, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
          write.table(matrix(c(txp_name, txp_name, "U"), nrow = 1), file = out_t2g3col, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
          write.table(ln, file = out_fa, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
        } else {
          # if not a header, just write to fasta file
          write.table(ln, file = out_fa, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
        }
      }
    }
  }
  
  
  
  
  message("Done.")
}

load_fry <- function(frydir, which_counts = c('S', 'A'), verbose = FALSE) {
  suppressPackageStartupMessages({
    library(rjson)
    library(Matrix)
    library(SingleCellExperiment)
  })
  
  # read in metadata
  qfile <- file.path(frydir, "quant.json")
  if (!file.exists(qfile)) {
    qfile <- file.path(frydir, "meta_info.json")
  }
  
  meta_info <- fromJSON(file = qfile)
  ng <- meta_info$num_genes
  usa_mode <- meta_info$usa_mode
  
  if (usa_mode) {
    if (length(which_counts) == 0) {
      stop("Please at least provide one status in 'U' 'S' 'A' ")
    }
    if (verbose) {
      message("processing input in USA mode, will return ", paste(which_counts, collapse = '+'))
    }
  } else if (verbose) {
    message("processing input in standard mode, will return spliced count")
  }
  
  # read in count matrix
  af_raw <- readMM(file = file.path(frydir, "alevin", "quants_mat.mtx"))
  # if usa mode, each gene gets 3 rows, so the actual number of genes is ng/3
  if (usa_mode) {
    if (ng %% 3 != 0) {
      stop("The number of quantified targets is not a multiple of 3")
    }
    ng <- as.integer(ng/3)
  }
  
  # read in gene name file and cell barcode file
  afg <- read.csv(file.path(frydir, "alevin", "quants_mat_cols.txt"), 
                  strip.white = TRUE, header = FALSE, nrows = ng, 
                  col.names = c("gene_ids"), row.names = 1)
  afc <- read.csv(file.path(frydir, "alevin", "quants_mat_rows.txt"), 
                  strip.white = TRUE, header = FALSE,
                  col.names = c("barcodes"), row.names = 1)
  
  # if in usa_mode, sum up counts in different status according to which_counts
  if (usa_mode) {
    rd <- list("S" = seq(1, ng), "U" =  seq(ng + 1, 2 * ng),
               "A" =  seq(2 * ng + 1, 3 * ng))
    o <- af_raw[, rd[[which_counts[1]]], drop = FALSE]
    for (wc in which_counts[-1]) {
      o <- o + af_raw[, rd[[wc]], drop = FALSE]
    }
  } else {
    o <- af_raw
  }
  
  # create SingleCellExperiment object
  sce <- SingleCellExperiment(list(counts = t(o)),
                              colData = afc,
                              rowData = afg
  )
  sce
}