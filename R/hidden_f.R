.check_extras = function(extras, paired, total.n){

  if(!('distr' %in% names(extras))){
    extras$distr = 'normal'
  }else{
    extras$distr = match.arg(extras$distr,
                             c('normal', 'empirical', 'custom'))
    if(extras$distr == 'custom' & !('custdens' %in% names(extras))){
      stop(.makepretty('to use custom fragment distribution, provide
                "custdens", a logspline object representing the distribution.'))
    }
  }

  # I don't love this--fraglen and fragsd aren't needed unless distr is normal.
  # but we store them anyway. should code better?
  if (!('fraglen' %in% names(extras))) {
    extras$fraglen = rep(250, total.n)
  } else {
    if (length(extras$fraglen) == 1) {
      extras$fraglen = rep(extras$fraglen, total.n)
    } else {
      stopifnot(length(extras$fraglen) == total.n)
    }
  }
  if (!('fragsd' %in% names(extras))) {
    extras$fragsd = rep(25, total.n)
  } else {
    if (length(extras$fragsd) == 1) {
      extras$fragsd = rep(extras$fragsd, total.n)
    } else {
      stopifnot(length(extras$fragsd) == total.n)
    }
  }

  if(!('readlen' %in% names(extras))){
    extras$readlen = 100
  }

  if(!('bias' %in% names(extras))){
    extras$bias = 'none'
  }else{
    extras$bias = match.arg(extras$bias, c('none', 'rnaf', 'cdnaf'))
  }

  if(!('error_model' %in% names(extras))){
    extras$error_model = 'uniform'
  }
  .check_error_model(extras, paired)

  if(!('error_rate' %in% names(extras))){
    extras$error_rate = 0.005
  }
  if(extras$error_model == 'custom'){
    extras$path = paste0(extras$model_path, '/', extras$model_prefix)
  }#this should work beause we already checked stuff.

  if(!('bias' %in% names(extras))){
    extras$bias = 'none'
  }else{
    extras$bias = match.arg(extras$bias, c('none', 'rnaf', 'cdnaf'))
  }

  if(!('lib_sizes' %in% names(extras))){
    extras$lib_sizes = rep(1, total.n)
  }else{
    stopifnot(is.numeric(extras$lib_sizes))
    stopifnot(length(extras$lib_sizes) == total.n)
  }

  if (!('frag_GC_bias' %in% names(extras))) {
    extras$frag_GC_bias <- 'none'
  } else {
    stopifnot(is.matrix(extras$frag_GC_bias))
    stopifnot(nrow(extras$frag_GC_bias) == 101)
    stopifnot(ncol(extras$frag_GC_bias) == total.n)
    stopifnot(all(extras$frag_GC_bias >= 0 & extras$frag_GC_bias <= 1))
  }

  if (!('strand_specific' %in% names(extras))) {
    extras$strand_specific <- FALSE
  }

  return(extras)

}

.check_error_model = function(extras, paired){

  # make sure it's an available model
  error_model = match.arg(extras$error_model,
                          c('uniform', 'illumina4', 'illumina5', 'custom'))

  # check uniform model --> error rate
  if(error_model == 'uniform'){
    if('error_rate' %in% names(extras)){
      error_rate = extras$error_rate
      stopifnot(is.numeric(error_rate))
      stopifnot(error_rate >= 0 & error_rate <= 1)
    }
  }

  # check paths and such for custom model
  if(error_model == 'custom'){
    if(!('model_path' %in% names(extras)) |
       !('model_prefix' %in% names(extras))){
      stop(.makepretty('with custom error models, you must provide both
                the path to the folder that holds your error model
                (model_path) and the prefix of your error model (model_prefix),
                where the prefix is whatever comes before _mate1 and _mate2
                (for paired reads) or _single (for single-end reads). You
                provided prefix when running build_error_models.py.'))
    }
    model_path = extras$model_path
    model_prefix = extras$model_prefix
    if(paired){
      if(!file.exists(paste0(model_path, '/', model_prefix, '_mate1')) |
         !file.exists(paste0(model_path, '/', model_prefix, '_mate2'))){
        stop('could not find error model')
      }
    }else{
      if(!file.exists(paste0(model_path, '/', model_prefix, '_single'))){
        stop('could not find error model')
      }
    }
  }
}

.check_fold_changes = function(fold_changes, num_reps, transcripts){

  # make sure fold change matrix is compatible with experiment size
  if(length(num_reps) == 1 | length(num_reps) == 2){
    stopifnot(is.numeric(fold_changes))
  }else{
    stopifnot(is.matrix(fold_changes))
    if(ncol(fold_changes) != length(num_reps)){
      stop(.makepretty('wrong number of columns in fold change matrix:
                need same number of columns as number of groups.'))
    }
    if(nrow(fold_changes) != length(transcripts)){
      stop(.makepretty('wrong number of rows in fold change matrix: need
                same number of rows as number of simulated transcripts. see
                count_transcripts to find out that number.'))
    }
  }
}

## internal sequencing function

sgseq = function(readmat, transcripts, paired, outdir, extras, reportCoverage=FALSE){
  #report theoretically perfect coverage if reportCoverage=TRUE, will write a file
  if(reportCoverage==TRUE){
    templates = unique(transcripts)
    coverage_matrices = list()
    for(i in 1:length(templates)){coverage_matrices = c(coverage_matrices, list(matrix(0, ncol=dim(readmat)[2], width(templates)[i])))}
    names(coverage_matrices) = names(templates)
  }

  for(i in seq_len(ncol(readmat))) {
    ##$ begin small chunk regarding fragment GC bias or not
    if (is.matrix(extras$frag_GC_bias)) {
      frag_GC_bias <- extras$frag_GC_bias[,i]
    } else {
      frag_GC_bias <- 'none'
    }
    ### end
    tObj = rep(transcripts, times=readmat[,i])
    iterations = ceiling(length(tObj) / 1e6L)
    offset = 1L
    for(iteration in seq_len(iterations)) {
      tSubset = tObj[offset:min(offset+999999L, length(tObj))] ## corrected value of integer added to offset to avoid duplicating reads
      tFrags = generate_fragments(tSubset, extras$fraglen[i], extras$fragsd[i],
                                  extras$readlen, extras$distr, extras$custdens,
                                  extras$bias, frag_GC_bias)

      if (!extras$strand_specific) {
        #reverse_complement some of those fragments
        tFrags = reverse_complement(tFrags)
      }

      #get reads from fragments
      reads = get_reads(tFrags, extras$readlen, paired)

      if(reportCoverage==TRUE){
        read_info = unique(names(reads))
        read_info_split = strsplit(read_info, ";mate1:|;mate2:")
        read_info_matrix = matrix(unlist(read_info_split), ncol=3, byrow=T)
        for(j in 1:dim(read_info_matrix)[1]){
          read = read_info_matrix[j,]
          target = which(names(coverage_matrices)==read[1])
          # ML: changing these to strsplit (str_split requires stringr depends or imports)
          coords1 = unlist(strsplit(read[2], "-"))
          coords2 = unlist(strsplit(read[3], "-"))
          coverage_matrices[[target]][coords1[1]:coords1[2],i]=coverage_matrices[[target]][coords1[1]:coords1[2],i]+1
          coverage_matrices[[target]][coords2[1]:coords2[2],i]=coverage_matrices[[target]][coords2[1]:coords2[2],i]+1
          save(coverage_matrices, file=file.path(outdir, 'sample_coverages.rda') )
        }
      }

      #add sequencing error
      if(extras$error_model == 'uniform'){
        errReads = add_error(reads, extras$error_rate)
      }else if(extras$error_model == 'custom'){
        errReads = add_platform_error(reads, 'custom', paired, extras$path)
      }else{
        errReads = add_platform_error(reads, extras$error_model, paired)
      }

      #write read pairs
      write_reads(errReads, readlen=extras$readlen,
                  fname=paste0(outdir, '/sample_', sprintf('%02d', i)), paired=paired,
                  gzip=extras$gzip, offset=offset)
      offset = offset + 1e6L
    }
  }

}

.write_info = function(extras, transcripts, num_reps, fold_changes, outdir,
                       group_ids, counts_matrix){

  if(!('transcriptid' %in% names(extras))){
    extras$transcriptid = names(transcripts)
  }

  if(is.numeric(fold_changes)){
    sim_info = data.frame(transcriptid=extras$transcriptid,
                          foldchange=fold_changes, DEstatus=fold_changes!=1)
  }else{
    fcv = apply(fold_changes, 1, function(x){
      paste(x, collapse=';')
    })
    DEstatus = rowSums(fold_changes) != ncol(fold_changes)
    sim_info = data.frame(transcriptid=extras$transcriptid,
                          foldchange=fcv, DEstatus=DEstatus)
  }

  write.table(sim_info, row.names=FALSE, quote=FALSE, sep="\t",
              file=paste0(outdir, '/sim_tx_info.txt'))

  rep_info = data.frame(
    rep_id=paste0('sample_', sprintf('%02d', 1:sum(num_reps))),
    group=group_ids, lib_sizes=extras$lib_sizes)

  write.table(rep_info, row.names=FALSE, quote=FALSE, sep='\t',
              file=paste0(outdir, '/sim_rep_info.txt'))

  rownames(counts_matrix) <- names(transcripts)
  colnames(counts_matrix) <- rep_info$rep_id
  save(counts_matrix, file=paste0(outdir, '/sim_counts_matrix.rda'))
}
