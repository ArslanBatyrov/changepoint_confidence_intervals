PELT = function(sumstat,pen=0, cost_func = "norm.mean", shape = 1, minseglen = 1, conf.set = FALSE){
  # function that uses the PELT method to calculate changes in mean where the segments in the data are assumed to be Normal
  n = length(sumstat[,1]) - 1
  if(n<2){stop('Data must have atleast 2 observations to fit a changepoint model.')}

  # conf.set is strictly TRUE/FALSE for now (scaffolding): FALSE = normal PELT
  # only, TRUE = also run the confidence-set execution in C. C sees 0 or 1.
  # (The confidence level itself stays on the R side for now.)
  # I must not forget to make a proposal to move the confidence-level part to C:
  # with high numbers of cpts and near-miss candidates it might overwhelm R and
  # could be faster within C
  if(!(isTRUE(conf.set) || isFALSE(conf.set))){stop('conf.set must be TRUE or FALSE.')}
  conf.flag = as.integer(isTRUE(conf.set))
  
  storage.mode(sumstat) = 'double'
  error=0
  
  lastchangelike = array(0,dim = n+1)
  lastchangecpts = array(0,dim = n+1)
  numchangecpts = array(0,dim = n+1)
  
  cptsout=rep(0,n) # sets up null vector for changepoint answer
  storage.mode(cptsout)='integer'
  
  answer=list()
  answer[[6]]=1
  on.exit(.C("FreePELT",answer[[6]]))
  
  storage.mode(lastchangelike) = 'double'
  storage.mode(lastchangecpts) = 'integer'
  storage.mode(numchangecpts) = 'integer'

  # conf_set output buffers for C
  # new design instead of flat vectors with lost grouping
  # we hand C a pre shaped matrix, row i = candidates of optimal cpt i, cols = candidate slots (Idea of Rebecca :)
  # so the grouping survives the C -> R trip for free
  # dims below are only a first guess; C writes the truly needed dims into cs_dim (arg 15)
  # and if the guess undershot in either direction we re run once with the exact size
  # (same logic as the v1 flat-vector version, just sending a matrix instead of a flat vector)
  nrow_guess = 10L
  ncol_guess = 30L

  # one C call to PELTC, building a helper function to make runs
  # simpler as we might call it twice in case our first guess was wrong
  run_peltc = function(nrows, ncols){
    .C('PELTC', cost_func, sumstat, as.integer(n), as.double(pen), cptsout,
       as.integer(error), as.double(shape), as.integer(minseglen),
       lastchangelike, lastchangecpts, numchangecpts, conf.flag,
       matrix(0L, nrows, ncols), matrix(0, nrows, ncols),
       as.integer(c(nrows, ncols)))
  }

  # answer=.C('PELT',cost_func, y3, y2,y,as.integer(n),as.double(pen),cptsout,as.integer(error),as.double(shape))
  answer = run_peltc(nrow_guess, ncol_guess)

  if(answer[[6]]>0){
    stop("C code error:",answer[[6]],call.=F)
  }
  if(conf.flag == 1){
    needed = answer[[15]]  # c(rows, cols) C reported it truly needed
    rows_have = nrow_guess # dims of the pass we are keeping
    cols_have = ncol_guess
    if(needed[1] > rows_have || needed[2] > cols_have){
      # If our first guess undershoot in either dimension: reallocate exact and re-run once.
      # Wrap the possible too large allocation in tryCatch
      # I did this as it was the way to make a universal upside cap
      # that is simple, as if we were to try to check for % of total RAM available
      # or what R has access to, we would make it ultra duper complicated and each device might need
      # its own approach, while using tryCatch will give the same outcome
      # asked whether this is a good idea in a convo with Rebecca and Owen and the idea was approved.
      answer2 = tryCatch(run_peltc(needed[1], needed[2]), error = function(e){
        warning("confidence-set buffer too large to allocate (", needed[1], " x ", needed[2],
                "); confidence set may be incomplete", call.=FALSE)
        NULL
      })
      if(!is.null(answer2)){
        if(answer2[[6]]>0){ stop("C code error:",answer2[[6]],call.=F) }
        answer = answer2
        rows_have = needed[1]
        cols_have = needed[2]
      }
    }
    # .C drops the dim attribute on return, so we re wrap the flat vector back
    # into the matrix shape and cut off the padding rows/cols we did not need
    pos = matrix(answer[[13]], nrow = rows_have, ncol = cols_have)
    lik = matrix(answer[[14]], nrow = rows_have, ncol = cols_have)
    keep_rows = seq_len(min(needed[1], rows_have))
    keep_cols = seq_len(min(needed[2], cols_have))
    return(list(lastchangecpts=answer[[10]],cpts=sort(answer[[5]][answer[[5]]>0]), lastchangelike=answer[[9]], ncpts=answer[[11]], checklist_positions=pos[keep_rows, keep_cols, drop=FALSE], checklist_likes=lik[keep_rows, keep_cols, drop=FALSE]))
  }
  return(list(lastchangecpts=answer[[10]],cpts=sort(answer[[5]][answer[[5]]>0]), lastchangelike=answer[[9]], ncpts=answer[[11]]))

}
