# for now I am creating the empty function that will be able to receive
# data from the C execution; for now it is empty as the C execution is not implemented.
# I will test the code below with a simple C function that will return the value and
# make sure that it does what it is supposed to.

PELT_CS = function(sumstat, pen = 0, cost_func = "mean.norm", shape = 1, minseglen = 1){

  # As I do not know how many changepoints we have, I need to first run PELT to get the number
  # of changepoints, so that I can allocate the required space for each.
  pelt = PELT(sumstat, pen = pen, cost_func = cost_func, shape = shape, minseglen = minseglen)

  n = length(sumstat[, 1]) - 1
  n_cpts = length(pelt$cpts)

  # my rationale (as we discussed with Rebecca) is that per changepoint I will allocate
  # the max memory that could potentially be required, which is equal to the number of
  # data points overall. As such, I find the number of changepoints and multiply them by
  # the number of data points. Just in case, I pad unused columns with 0.
  # So the data space would be n_cpts * n columns. (n_cpts is always >= 1.)

  checklist_positions = matrix(0L, nrow = n_cpts, ncol = n) # candidate changepoint locations
  checklist_likes = matrix(0, nrow = n_cpts, ncol = n) # their (pure) likelihood values

  # Label each bucket with the exact C type it will hold.
  storage.mode(checklist_positions) = "integer"
  storage.mode(checklist_likes) = "double"

  # here I will add a C call later on,
  # so that afterwards we can return the now-filled buckets back out of answer.

  return(list(
    cpts = pelt$cpts,
    lastchangecpts = pelt$lastchangecpts, # backpointer table
    checklist_positions = checklist_positions, # empty space for C to fill
    checklist_likes = checklist_likes # empty space for C to fill
  ))
}
