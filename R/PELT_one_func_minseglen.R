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

  # conf_set output buffers for C will be allocated here
  # the previous logic of minseglen^2 or n * ncpts both were inefficient
  # So, currently the next method we are going to use is two fold
  #  1: max(We allocate minseglen^2, 30) first
  #  2: if it is enough -> Operation is done, if not:
  #  3: C rewrites the cs_len and hands it to R, so that R knows exactly
  #  the memory that is required to be allocated. Second run occurs
  buflen = max(as.integer(minseglen) * as.integer(minseglen), 30L)

  # one C call to PELTC, building a helper function to make runs
  # simpler as we might call it twice in case our first guess memory allocation was wrong
  # meaning cs_len < needed
  run_peltc = function(buflen){
    .C('PELTC', cost_func, sumstat, as.integer(n), as.double(pen), cptsout,
       as.integer(error), as.double(shape), as.integer(minseglen),
       lastchangelike, lastchangecpts, numchangecpts, conf.flag,
       rep(0L, buflen), rep(0, buflen), as.integer(buflen))
  }

  # answer=.C('PELT',cost_func, y3, y2,y,as.integer(n),as.double(pen),cptsout,as.integer(error),as.double(shape))
  answer = run_peltc(buflen)

  if(answer[[6]]>0){
    stop("C code error:",answer[[6]],call.=F)
  }
  if(conf.flag == 1){
    needed = answer[[15]]          # count C reported it needed
    captured = min(needed, buflen) # how many the first pass wrote
    if(needed > buflen){
      # If our first guess undershoot: reallocate to the exact size and re-run once.
      # Wrap the possible too large allocation in tryCatch
      # I did this as it was the way to make a universal upside cap
      # that is simple, as if we were to try to check for % of total RAM available
      # or what R has access to, we would make it excessively complicated and each device might need
      # its own approach, while using tryCatch will give the same outcome
      # ask whether this is a good idea in a convo with Rebecca and Owen
      answer2 = tryCatch(run_peltc(needed), error = function(e){
        warning("confidence-set buffer too large to allocate (", needed,
                " slots); confidence set may be incomplete", call.=FALSE)
        NULL
      })
      if(!is.null(answer2)){
        if(answer2[[6]]>0){ stop("C code error:",answer2[[6]],call.=F) }
        answer = answer2
        captured = needed
      }
    }
    # trim off the unused padding so R only sees the real captured entries
    return(list(lastchangecpts=answer[[10]],cpts=sort(answer[[5]][answer[[5]]>0]), lastchangelike=answer[[9]], ncpts=answer[[11]], checklist_positions=answer[[13]][seq_len(captured)], checklist_likes=answer[[14]][seq_len(captured)]))
  }
  return(list(lastchangecpts=answer[[10]],cpts=sort(answer[[5]][answer[[5]]>0]), lastchangelike=answer[[9]], ncpts=answer[[11]]))

}
