# Deleted the PELT_CS as the design changed and we need to work within PELTC.
# conf_set, as agreed, should be just a function that takes all of the values
# required to calculate the confidence sets and then calculates it.
# Also, a little later on, I will make confint call this.
# We will have two functions doing the same thing to make it more comfortable for users.
# In the docs, we will describe conf_set only.

confidence_set = function(object, level = 0.95) {
  if(!is(object, "cpt")){stop("object must be a cpt or cpt.range, cpt.reg is not supported")}
  if(object@method != "PELT"){stop("confidence sets are currently only implemented for method = 'PELT'")}

  # Turning the human-readable test.stat into a C-readable cost-function key
  if(test.stat(object) != "Normal"){stop("confidence set supports test.stat = 'Normal' only")}
  costfunc = switch(cpttype(object),
                    "mean" = "mean.norm",
                    "variance" = "var.norm",
                    "mean and variance" = "meanvar.norm",
                    stop("Unsupported changepoint type: ", cpttype(object)))
  if(pen.type(object) == "MBIC"){costfunc = paste0(costfunc, ".mbic")}

  # I copied the sumstat formula from data_input.R for use here
  data = data.set(object)
  mu = mean(data)
  sumstat = cbind(c(0, cumsum(coredata(data))), c(0, cumsum(coredata(data)^2)), cumsum(c(0, (coredata(data) - mu)^2)))

  # basic validation of the level here, C does not see it for now. In C (PELTC) we only have an off/on button
  if(!is.numeric(level) || length(level) != 1 || is.na(level) || level <= 0 || level >= 1){
    stop("level must be a single number strictly between 0 and 1") # this is not used by C yet
  }

  # Re-run PELT with the conf_set on
  pelt = PELT(sumstat, pen = pen.value(object), cost_func = costfunc, minseglen = object@minseglen, conf.set = TRUE)

  # just a sanity check that the re-run agrees with the changepoints we already stored
  if(!all(cpts(object) %in% pelt$cpts)){
    warning("re-run PELT did not reproduce the object's changepoints; confidence set may not correspond to the fit")
  }
  return(list(
    level = level,
    cpts = pelt$cpts,
    lastchangecpts = pelt$lastchangecpts,
    checklist_positions = pelt$checklist_positions, # fake pattern for now, will be real later
    checklist_likes = pelt$checklist_likes # fake pattern as well for now
  ))
}

# not a fully complete code, must be reviewed
# confidence set core, ported from Owens code (PELT_CS_fns.R).
# all internal pure R, costs arive alredy calcualted and we never recompute them.
# indexing: positon t is stored at slot t+1 (R has no slot 0)

# runmin from Owens code, renamed with cs_ to aviod clashes
cs_runmin = function(v, k){
  out = NULL
  for(i in 1:(length(v)-k+1)){
    out = c(out, min(v[i:(i+k-1)], na.rm=TRUE))
  }
  return(out)
}

# lookup table: how many optimal cpts stand before each positon,
# 0 is seeded with -1 so candiate 0 pays no penatly
cs_no_op_cpts = function(lastchangecpts){
  n = length(lastchangecpts) - 1
  no_op = rep(NA_integer_, n + 1)
  no_op[1] = -1
  for(t in 1:n){
    no_op[t+1] = no_op[lastchangecpts[t+1] + 1] + 1
  }
  return(no_op)
}

# confidnece level -> cost bugdet, fromula from Owens code as it is.
# smaller alpha -> bigger set (CS1)
cs_penalty = function(op_cpts, alpha){
  cpts = c(0, op_cpts)
  n = op_cpts[length(op_cpts)]
  if(length(cpts) == 2){ # no changepoints detected (op_cpts holds just n)
    nj.maxi = n/20
  } else {
    nj.maxi = min(max(cs_runmin(diff(cpts), 2)), n/4 - n*5/24)
  }
  return(3/2 + 2*log(2*nj.maxi) - 2*log(alpha))
}

# keeps every near miss candiate wihtin the budget from the optimal,
# row i of CS = candidates for the cpt before op_cpts[i], NA padded
cs_core = function(op_cpts, checklists, no_op_cpts, penalty, level=0.95, cs.penalty=NULL){
  alpha = 1 - level
  cpts = c(0, op_cpts)
  n = op_cpts[length(op_cpts)]
  if(is.null(cs.penalty)){ cs.penalty = cs_penalty(op_cpts, alpha) }

  CS = CS.likes = matrix(NA_real_, nrow=length(cpts)-1, ncol=n)
  for(i in length(cpts):2){
    x = checklists[[i-1]]$positions # near miss candidates at time cpts[i]
    y = checklists[[i-1]]$likes # thier pure costs
    # pure cost + one penalty per cpt on the candidates path
    semipurelikes = y + (no_op_cpts[x+1] + 1) * penalty
    op.index = which(x == cpts[i-1])
    if(length(op.index) != 1){
      stop("optimal changepoint ", cpts[i-1], " not found in the checklist at time ", cpts[i])
    }
    op.like = semipurelikes[op.index]
    CSvals.index = which(semipurelikes >= op.like & semipurelikes <= op.like + cs.penalty)
    CS.likes[i-1, seq_along(CSvals.index)] = semipurelikes[CSvals.index]
    CS[i-1, seq_along(CSvals.index)] = x[CSvals.index]
  }
  return(list(CS=CS, CS.likes=CS.likes, cs.penalty=cs.penalty))
}

# the "bunny hop", walks the backpionters back untill 0
cs_op_seg = function(fcpt=NULL, last, lastchangecpts){
  while(last != 0){
    fcpt = c(fcpt, last)
    last = lastchangecpts[last + 1]
  }
  return(sort(fcpt))
}

# CS matrix -> whole alternitive segmentaions, out[[k]] = the ones with k cpts
# (n counts as one aswell)
cs_backtrack = function(lastchangecpts, CS_mat, est_cpts){
  out = list()
  out[[length(est_cpts)]] = est_cpts
  fcpt = NULL
  for(cpt.i in nrow(CS_mat):1){
    fcpt = c(fcpt, est_cpts[cpt.i]) # fixing the optimal cpts later then this row
    CS.vec = CS_mat[cpt.i, !is.na(CS_mat[cpt.i,])]
    for(CS.vec.i in seq_along(CS.vec)){
      if(CS.vec[CS.vec.i] == 0){
        # candidate 0 = no earlier cpt, segmentation is just fcpt
        out[[length(fcpt)]] = rbind(out[[length(fcpt)]], sort(fcpt))
      } else {
        seg = cs_op_seg(fcpt=fcpt, last=CS.vec[CS.vec.i], lastchangecpts=lastchangecpts)
        if(length(seg) > length(out)){
          out[[length(seg)]] = seg
        } else {
          out[[length(seg)]] = rbind(out[[length(seg)]], seg)
        }
      }
    }
  }
  return(out)
}
