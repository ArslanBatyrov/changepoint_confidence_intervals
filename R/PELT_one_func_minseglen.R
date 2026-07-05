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

  # empty space for the confidence-set execution in C to fill (only used when
  # conf.set is TRUE, but .C always needs the same number of arguments).
  # minseglen*minseglen, NOT n_cpts*n on a large dataset (I will need to chat with
  # Owen and Rebecca as this way of allocating memory has too many crash edge cases)
  # I will propose adding one more variable that both inputs to AND outputs from C:
  # the first one is a guess, in case it is not enough we reallocate space with the
  # exact number C requires
  buflen = as.integer(minseglen) * as.integer(minseglen)
  checklist_positions = rep(0L, buflen)
  checklist_likes = rep(0, buflen)
  storage.mode(checklist_positions) = 'integer'
  storage.mode(checklist_likes) = 'double'

  # answer=.C('PELT',cost_func, y3, y2,y,as.integer(n),as.double(pen),cptsout,as.integer(error),as.double(shape))
  answer=.C('PELTC',cost_func, sumstat,as.integer(n),as.double(pen),cptsout,as.integer(error),as.double(shape), as.integer(minseglen), lastchangelike, lastchangecpts,numchangecpts, conf.flag, checklist_positions, checklist_likes)

  if(answer[[6]]>0){
    stop("C code error:",answer[[6]],call.=F)
  }
  if(conf.flag == 1){
    # confidence set requested: also hand back the raw material C captured
    return(list(lastchangecpts=answer[[10]],cpts=sort(answer[[5]][answer[[5]]>0]), lastchangelike=answer[[9]], ncpts=answer[[11]], checklist_positions=answer[[13]], checklist_likes=answer[[14]]))
  }
  return(list(lastchangecpts=answer[[10]],cpts=sort(answer[[5]][answer[[5]]>0]), lastchangelike=answer[[9]], ncpts=answer[[11]]))

}
