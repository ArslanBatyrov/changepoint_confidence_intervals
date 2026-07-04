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
