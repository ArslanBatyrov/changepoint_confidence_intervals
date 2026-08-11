# Test-only mock for the confidence-set ingredients.
# A pure-R PELT (ported from the Owen's prototype PELT_lv) that KEEPS the
# per-time checklist of near-miss candidates and their pure costs, exactly
# the ingredients the real C capture will eventually dump.
# Differences from the prototype, both deliberate:
#  1 - cost is the package's mll_mean kernel (RSS around the segment mean),
#    NOT the prototype's version with the n*log(2*pi*sigma^2) term. At any
#    fixed time point that term is a constant shift across all candidates,
#    so which candidates enter the confidence set is unchanged, and this
#    matches what the C code will hand us later.
#  2 - only the checklists AT the optimal changepoints are returned (that is
#    all cs_core needs), not at every time point. The rest would be deadweight if returned as only the optimal cpts and points around them matter.
# minseglen is 1 (the prototype has no minseglen notion).
mock_pelt_checklist = function(data, penalty){
  n = length(data)
  lastchangeF = -penalty  # penalised optimal cost up to time t (entry t+1) ALSO, naming of F I took from the paper :)
  lastchangelike = 0      # PURE optimal cost up to time t, tracked for the CS
  lastchangecpts = 0      # backpointer table (entry t+1)
  checklist = NULL
  checklist.t = vector("list", n + 1)

  for(tstar in 1:n){
    tmpt = c(checklist, tstar - 1)
    tmplike = tmplike_CS = numeric(length(tmpt))
    for(i in seq_along(tmpt)){
      seg = data[(tmpt[i]+1):tstar]
      cost = sum((seg - mean(seg))^2) # the mll_mean kernel: x2 - x1^2/len
      tmplike[i]    = lastchangeF[tmpt[i]+1] + cost + penalty
      tmplike_CS[i] = lastchangelike[tmpt[i]+1] + cost
    }
    lastchangeF[tstar+1] = min(tmplike)
    lastchangecpts[tstar+1] = tmpt[which.min(tmplike)]
    keep = tmplike <= lastchangeF[tstar+1] + penalty # PELT pruning rule
    checklist = tmpt[keep]
    lastchangelike[tstar+1] = tmplike_CS[which.min(tmplike)]
    checklist.t[[tstar+1]] = list(positions=checklist, likes=tmplike_CS[keep])
  }

  # backtrack the optimal changepoints (n included), as PELT does
  fcpt = NULL
  last = n
  while(last != 0){
    fcpt = c(fcpt, last)
    last = lastchangecpts[last+1]
  }
  op_cpts = sort(fcpt)

  list(op_cpts = op_cpts,
       lastchangecpts = lastchangecpts,
       no_op_cpts = changepoint:::cs_no_op_cpts(lastchangecpts),
       checklists = lapply(op_cpts, function(t) checklist.t[[t+1]]),
       penalty = penalty,
       n = n)
}
