blast.clipped = function(df, indelScore = 80, pctAlign = 90, blastRef){
  library(rBLAST)
  
  # Align clipped sequences to reference
  bl = blast(db = blastRef)
  
  # Remove NA sequences before creating DNAStringSet object
  df = df[!is.na(df$clipped.consensus), ]
  
  # Create named DNAStringSet object
  seq = df$clipped.consensus
  names(seq) = df$rname_clippedPos_Orientation_ReadSide
  seq = DNAStringSet(seq, use.names = TRUE)
  
  # Perform BLAST
  results = predict(bl, seq, BLAST_args = "-dust no")
  names(results) = c("query", "subject", "pct_identity", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  results$qlen = results$qend - results$qstart + 1
  
  # Filter BLAST results by bitscore
  results = results[results$bitscore > indelScore, ]
  
  ##################################
  # Compare to mapping of anchored portion of read
  reference.clusters = data.frame(
    chr = df$RNAME,
    reference.start = df$clipped_pos,
    reference.stop = df$clipped_pos,
    clipped_seq = df$clipped.consensus,
    rname_clippedPos_Orientation_ReadSide = df$rname_clippedPos_Orientation_ReadSide,
    counts = df$counts,
    strand = "+",
    stringsAsFactors = FALSE
  )
  
  reference.clusters = merge(reference.clusters, results, by.x = "rname_clippedPos_Orientation_ReadSide", by.y = "query", all.x = TRUE)
  
  # Filter by percent length of alignment
  reference.clusters = reference.clusters[(nchar(reference.clusters$clipped_seq) / reference.clusters$length) > pctAlign, ]
  
  return(reference.clusters)
}

best.hits = function(reference.clusters){
  # Filter out rows with NA values in critical columns
  reference.clusters = reference.clusters[!is.na(reference.clusters$subject) & reference.clusters$chr == reference.clusters$subject, ]
  
  if(nrow(reference.clusters) > 0){
    reference.clusters$hit.len = reference.clusters$send - reference.clusters$sstart
    reference.clusters$pct_aligned = 100 * (reference.clusters$qlen) / reference.clusters$length
    reference.clusters$aligned.strand = sign(reference.clusters$send - reference.clusters$sstart)
    reference.clusters$dist.to.alignment.start = reference.clusters$sstart - reference.clusters$reference.start
    reference.clusters$dist.to.alignment.end = reference.clusters$send - reference.clusters$reference.start
    reference.clusters$sign.same = sign(reference.clusters$dist.to.alignment.end) == sign(reference.clusters$dist.to.alignment.start)
    
    # Remove alignments that overlapped anchored portion of read
    reference.clusters = reference.clusters[reference.clusters$sign.same, ]
    
    if(nrow(reference.clusters) == 0){
      return(reference.clusters)
    }
    
    reference.clusters$start.dist = abs(reference.clusters$sstart - reference.clusters$reference.start)
    
    ### Keep hits with highest score, closest to clipped read
    reference.clusters$min.dist = sapply(1:nrow(reference.clusters), function(i) min(c(abs(reference.clusters$dist.to.alignment.end[i]), abs(reference.clusters$dist.to.alignment.start[i]))))
    reference.clusters = reference.clusters[order(reference.clusters$rname_clippedPos_Orientation_ReadSide, reference.clusters$evalue, reference.clusters$min.dist), ]
    reference.clusters = reference.clusters[!duplicated(reference.clusters$rname_clippedPos_Orientation_ReadSide), ]
  }
  
  return(reference.clusters)
}
