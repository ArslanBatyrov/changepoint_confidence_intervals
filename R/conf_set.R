# Deleted the PELT_CS as the design changed and we need to work within PELTC.
# conf_set, as agreed, should be just a function that takes all of the values
# required to calculate the confidence sets and then calculates it.
# Also, a little later on, I will make confint call this.
# We will have two functions doing the same thing to make it more comfortable for users.
# In the docs, we will describe conf_set only.

confidence_set = function(object, level = 0.95) {
  # check the object BEFORE touching the slot, otherwise conf.set() below blows up first. Added after I found out the bug and fixed it. 
  # with an ugly S4 dispatch error instead of our own message
  cs_check_object(object)
  cs_check_level(level)

  # use the slot as a cache. the pricey bit (PELT + capture) is in cs_ingredients and
  # doesnt care about the level, so do it once and stash it under "ingredients", reuse after.
  # if the fit changed (cpts/penalty dont match) just redo it, dont serve old leftovers
  existing = conf.set(object)
  ing = existing[["ingredients"]]
  if(is.null(ing) || !all(cpts(object) %in% ing$cpts) || !identical(ing$penalty, pen.value(object))){
    ing = cs_ingredients(object)
    existing[["ingredients"]] = ing
  }

  # the actual confidence set calculation, all validated against Owens code
  # CS / CS.likes / cs.penalty stay internal (backtrack needs them), the user does not see them
  # (MIGHT NEED TO ADD CS.likes to the output in the future though!)
  core = cs_core(ing$cpts, ing$checklists, cs_no_op_cpts(ing$lcc), ing$penalty, level = level)
  segs = cs_backtrack(ing$lcc, core$CS, ing$cpts)

  # tidy the segmentations for the user.
  # cs_backtrack gives segs[[k]] = the segmentations with k cpts, but with NULL holes at the ks that
  # have none and leftover rbind rownames, so: drop the holes, force every entry into a one row per
  # segmentation matrix, clear the junk rownames and name each entry by its cpt count (n counts as one)
  segs = segs[!vapply(segs, is.null, logical(1))]
  segs = lapply(segs, function(s){
    if(is.null(dim(s))){ s = matrix(s, nrow = 1) }
    dimnames(s) = NULL # drop the junk rbind rownames entirely, no col names either
    # "unique" keeps first occurence order and the backtrack seeds the optimal first,
    # thus the first row is ALWAYS a row of the optimal cpts
    unique(s)
  })
  names(segs) = vapply(segs, ncol, integer(1))

  out = list(
    level = level,
    cpts = ing$cpts, # optimal cpts, n included as the last one
    segmentations = segs # [["k"]] = all plausible segmentations with k cpts
  )
  # slap a tiny S3 class on top just so it prints nice, underneath it is still a plain list
  # and segmentations still has the whole conf set in it, optimal included, nothing dropped
  class(out) = "cpt.confset"

  # save under the level key and give back the whole object, not the set. thats what makes
  # confint just a twin of this. R wont edit in place so it only sticks if they reassign,
  # fit = confidence_set(fit, 0.95), then read it back with conf.set(fit)[["0.95"]]
  existing[[as.character(level)]] = out
  conf.set(object) = existing
  return(object)
}

# the half that doesnt need the level: checks, re-run PELT with capture on, shape the C output.
# level never shows up here so one run covers every level, thats the bit worth caching
cs_ingredients = function(object){
  cs_check_object(object)

  # Turning the human-readable test.stat into a C-readable cost-function key
  # opened up past Normal on mentors direction: we only wire the strings here (dist agnostic),
  # Rebecca will go over the statistical side of the non Normal dists herself later
  costfunc = switch(test.stat(object),
                    "Normal" = switch(cpttype(object),
                                      "mean" = "mean.norm",
                                      "variance" = "var.norm",
                                      "mean and variance" = "meanvar.norm",
                                      stop("Unsupported changepoint type: ", cpttype(object))),
                    "Exponential" = "meanvar.exp",
                    "Gamma" = "meanvar.gamma",
                    "Poisson" = "meanvar.poisson",
                    stop("Unsupported test statistic: ", test.stat(object)))
  if(pen.type(object) == "MBIC"){costfunc = paste0(costfunc, ".mbic")}

  # I copied the sumstat formula from data_input.R for use here
  data = data.set(object)
  mu = mean(data)
  sumstat = cbind(c(0, cumsum(coredata(data))), c(0, cumsum(coredata(data)^2)), cumsum(c(0, (coredata(data) - mu)^2)))

  # gamma is the only cost function in C that actually reads shape, the fit stores it
  # inside param.est (see the param method in cpt.class.R) so we can recover it from the object
  shape = 1
  if(test.stat(object) == "Gamma"){
    shape = param.est(object)$shape
    if(is.null(shape)){stop("cannot recover the gamma shape parameter from the object, refit with param.estimates = TRUE")}
  }

  # Re-run PELT with the conf_set on
  pelt = PELT(sumstat, pen = pen.value(object), cost_func = costfunc, shape = shape, minseglen = object@minseglen, conf.set = TRUE)

  # just a sanity check that the re-run agrees with the changepoints we already stored
  if(!all(cpts(object) %in% pelt$cpts)){
    warning("re-run PELT did not reproduce the object's changepoints; confidence set may not correspond to the fit")
  }

  # matrix rows from C -> cs_core style checklists, NA padding stripped
  # row i = near miss candidates of optimal cpt i and thier pure costs
  checklists = lapply(seq_len(nrow(pelt$checklist_positions)), function(i){
    keep = !is.na(pelt$checklist_positions[i, ])
    list(positions = as.numeric(pelt$checklist_positions[i, keep]),
         likes = as.numeric(pelt$checklist_likes[i, keep]))
  })

  return(list(
    cpts = pelt$cpts,
    checklists = checklists,
    lcc = as.numeric(pelt$lastchangecpts),
    penalty = pen.value(object)
  ))
}

# the object guards (safety checks), made own function so confidence_set can run them before it reads the slot efficiently. 
cs_check_object = function(object){
  if(!is(object, "cpt")){stop("object must be a cpt object, cpt.reg is not supported")}
  # cpt.range slips through the method check (CROPS stores method = PELT) but carries a
  # penalty RANGE not a single value, the maths would run on garbage. MUST add a cpt_range soon
  if(is(object, "cpt.range")){stop("confidence sets are not yet supported for BinSeg, SegNeigh or CROPS fits (cpt.range objects)")}
  if(object@method != "PELT"){stop("confidence sets are currently only implemented for method = 'PELT'")}
}

# quick level check, C never sees it (PELTC just has an on/off)
cs_check_level = function(level){
  if(!is.numeric(level) || length(level) != 1 || is.na(level) || level <= 0 || level >= 1){
    stop("level must be a single number strictly between 0 and 1") # this is not used by C yet
  }
}

# printing the WHOLE slot, ie what you get from conf.set(fit) with no key. hides the
# ingredients cache (its our internal plumbing, a huge blob of numbers, no use to anyone)
# and prints each stored level through the single set method below
print.cpt.confsets = function(x, n = getOption("changepoint.confset.n", 20), ...){
  lv = setdiff(names(x), "ingredients")
  if(length(lv) == 0){
    cat("No confidence sets stored yet, use confidence_set(fit, level) to add one.\n")
    return(invisible(x))
  }
  lv = lv[order(as.numeric(lv))] # nicer to read low to high, not the order they were added
  cat("Confidence sets stored at ", length(lv), " level(s): ", paste(lv, collapse = ", "), "\n", sep = "")
  for(l in lv){
    cat("\n")
    print(x[[l]], n = n) # hand the cap down so print(conf.set(fit), n = Inf) works too
  }
  if("ingredients" %in% names(x)){
    cat("\n(raw PELT output is cached in here too, hidden as its not for the user)\n") # This just a transparency note, so not necessary but just to let know.
  }
  invisible(x)
}

# printing the confidence set: optimal cpts first, then everything else as alternatives.
# the optimal is still inside segmentations, we just pull it out for the display
# n is how many alternatives to show. In the v1 I had a hard cap of 20, was pretty dumb as we need to give users flexibility, so I created a soft cap. It's just made so that a big set cannot flood the
# console, print(x, n = 50) or n = Inf shows more. the default can also be moved for good as 20 is in no way derived from a special calculation, just what seemed reasonable at the moment.
# with options(changepoint.confset.n = 50) so you dont have to keep typing it
print.cpt.confset = function(x, n = getOption("changepoint.confset.n", 20), ...){
  cat("Confidence set (level = ", x$level, ")\n\n", sep = "")
  cat("Optimal changepoints:\n  ", paste(x$cpts, collapse = " "), "\n", sep = "")

  # walk every segmentation group and keep the rows that are not the optimal one
  alts = list()
  for(g in x$segmentations){
    for(r in seq_len(nrow(g))){
      if(!identical(as.numeric(g[r, ]), as.numeric(x$cpts))){
        alts[[length(alts) + 1]] = g[r, ]
      }
    }
  }

  if(length(alts) == 0){
    cat("\nNo alternative segmentations at this level.\n")
  } else {
    cat("\nAlternative segmentations (", length(alts), "):\n", sep = "")
    # only a display cap, the data is all still there in $segmentations
    show_n = min(length(alts), n)
    for(i in seq_len(show_n)){
      cat("  ", paste(alts[[i]], collapse = " "), "\n", sep = "")
    }
    if(length(alts) > show_n){
      # tell them exactly what to type, otherwise they have no way of knowing they can
      cat("  ... and ", length(alts) - show_n, " more, use print(x, n = Inf) to see all\n", sep = "")
    }
  }
  invisible(x)
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
