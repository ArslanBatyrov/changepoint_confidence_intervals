# core tests: the maths only, nothing user facing.
# small helpers are checked against expected values computed from the formulas by hand, so
# the checks are independent of the code. the bigger blocks run on the R mock in
# helper-confset.R. one block runs the mock and the real C engine on the same data and
# they must agree, that is the whole point of having a second independent implementation!!
# I might have missed some tests that were needed and potentially add new tests in the future if such need occurs.

test_that("cs_runmin computes running minima", {
  expect_equal(changepoint:::cs_runmin(c(3,1,4,1,5), 2), c(1,1,1,1))
  expect_equal(changepoint:::cs_runmin(c(5,4,3,2), 3), c(3,2))
})

test_that("cs_no_op_cpts reconstructs the prototype's counts", {
  # one changepoint at 3: positions 1-3 have 0 cpts before them, 4-6 have 1
  lcc = c(0,0,0,0,3,3,3) # backpointer table (index t+1)
  expect_equal(changepoint:::cs_no_op_cpts(lcc), c(-1, 0,0,0, 1,1,1))
})

test_that("cs_penalty matches the prototype formula", {
  # write the formula out by hand so a "simplification" of cs_penalty would fail here.
  # no changepoints uses nj.maxi = n/20
  expect_equal(changepoint:::cs_penalty(op_cpts=c(100), alpha=0.05),
               3/2 + 2*log(2*100/20) - 2*log(0.05))
  # one changepoint uses nj.maxi = min(max(runmin(diff,2)), n/4 - 5n/24)
  expect_equal(changepoint:::cs_penalty(op_cpts=c(50,100), alpha=0.05),
               3/2 + 2*log(2*(100/4 - 100*5/24)) - 2*log(0.05))
})

test_that("cs_core keeps exactly the candidates within the budget", {
  # candidates 0, 4, 6 with costs 20, 18, 30 and penalty 5 -> semipure 20, 23, 35.
  # optimal (0) costs 20, budget 10, so it keeps 20 and 23 and drops 35.
  no_op = c(-1, rep(0, 10))
  res = changepoint:::cs_core(op_cpts = c(10),
                              checklists = list(list(positions=c(0,4,6), likes=c(20,18,30))),
                              no_op_cpts = no_op, penalty = 5,
                              cs.penalty = 10)
  expect_equal(res$CS[1, 1:2], c(0, 4))
  expect_true(all(is.na(res$CS[1, 3:10])))
  expect_equal(res$CS.likes[1, 1:2], c(20, 23))
  expect_equal(res$cs.penalty, 10)
  # the gap is semipure minus the optimal (20): 20-20=0 for the optimal, 23-20=3 for the other
  expect_equal(res$CS.gaps[1, 1:2], c(0, 3))
})

test_that("cs_core gaps are non-negative, within budget, and 0 for the optimal", {
  set.seed(7)
  data = c(rnorm(50), rnorm(50)+2.5)
  mock = mock_pelt_checklist(data, 2*log(100))
  core = changepoint:::cs_core(mock$op_cpts, mock$checklists, mock$no_op_cpts,
                               mock$penalty, level=0.95)
  g = core$CS.gaps[!is.na(core$CS.gaps)]
  expect_true(all(g >= 0))
  expect_true(all(g <= core$cs.penalty + 1e-9))
  # the optimal candidate sits in every row and its gap is exactly 0
  expect_true(all(apply(core$CS.gaps, 1, function(r) any(abs(r[!is.na(r)]) < 1e-9))))
})

test_that("cs_core errors if the optimal changepoint is missing from a checklist", {
  expect_error(
    changepoint:::cs_core(op_cpts = c(10),
                          checklists = list(list(positions=c(4,6), likes=c(18,30))),
                          no_op_cpts = c(-1, rep(0,10)), penalty = 5, cs.penalty = 10),
    "not found in the checklist")
})

test_that("cs_op_seg walks the backpointer chain", {
  # chain: 6 -> 3 -> 0
  lcc = c(0, 0,0,0, 0,0, 3)
  expect_equal(changepoint:::cs_op_seg(fcpt=c(10), last=6, lastchangecpts=lcc),
               c(3, 6, 10))
})

test_that("cs_backtrack rebuilds whole segmentations from the CS matrix", {
  # optimal cpts c(6,10); backpointer chain 4 -> 2 -> 0, everything else -> 0
  lcc = c(0,0,0,0,2,0,0,0,0,0,0)
  CS_mat = matrix(NA_real_, nrow=2, ncol=10)
  CS_mat[2, 1:2] = c(6, 4) # near-miss candidates for the cpt before 10
  CS_mat[1, 1:2] = c(0, 2) # near-miss candidates for the cpt before 6
  out = changepoint:::cs_backtrack(lcc, CS_mat, est_cpts=c(6,10))
  # the only 2-cpt segmentation is the optimal itself
  expect_equal(unname(unique(out[[2]])), matrix(c(6,10), nrow=1))
  # 3-cpt alternatives: candidate 4 drags in its chain (2,4,10), candidate 2 gives (2,6,10)
  expect_equal(as.numeric(out[[3]][1,]), c(2,4,10))
  expect_equal(as.numeric(out[[3]][2,]), c(2,6,10))
})

test_that("mock PELT agrees with the C engine, capture included", {
  set.seed(42)
  data = c(rnorm(40), rnorm(40)+3, rnorm(40)-2)
  pen = 2*log(120)
  mock = mock_pelt_checklist(data, pen)

  mu = mean(data)
  sumstat = cbind(c(0, cumsum(data)), c(0, cumsum(data^2)), cumsum(c(0, (data-mu)^2)))
  pelt = PELT(sumstat, pen=pen, cost_func="mean.norm", minseglen=1, conf.set=TRUE)

  # first the ordinary PELT answer
  expect_equal(mock$op_cpts, pelt$cpts)
  expect_equal(mock$lastchangecpts, as.numeric(pelt$lastchangecpts))

  # then the bit the mock actually exists for: the near-miss candidates C captured.
  # this is the only place the C capture is checked against a separate implementation
  expect_equal(nrow(pelt$checklist_positions), length(mock$checklists))
  for(i in seq_along(mock$checklists)){
    keep = !is.na(pelt$checklist_positions[i, ])
    expect_equal(as.numeric(pelt$checklist_positions[i, keep]),
                 as.numeric(mock$checklists[[i]]$positions))
    expect_equal(as.numeric(pelt$checklist_likes[i, keep]),
                 as.numeric(mock$checklists[[i]]$likes))
  }
})

test_that("optimal segmentation is always inside its own confidence set", {
  set.seed(7)
  data = c(rnorm(50), rnorm(50)+2.5)
  mock = mock_pelt_checklist(data, 2*log(100))
  core = changepoint:::cs_core(mock$op_cpts, mock$checklists, mock$no_op_cpts,
                               mock$penalty, level=0.95)
  # each optimal changepoint must appear in its own row of the set
  cpts0 = c(0, mock$op_cpts)
  for(i in seq_len(nrow(core$CS))){
    expect_true(cpts0[i] %in% core$CS[i, ])
  }
  # and the whole optimal segmentation must be one of the rebuilt ones
  segs = changepoint:::cs_backtrack(mock$lastchangecpts, core$CS, mock$op_cpts)
  k = length(mock$op_cpts)
  found = any(apply(rbind(segs[[k]]), 1, function(r) identical(as.numeric(r), as.numeric(mock$op_cpts))))
  expect_true(found)
  # sanity: every segmentation ends at n and is strictly increasing
  for(m in segs){
    if(is.null(m)) next
    m = rbind(m)
    for(r in seq_len(nrow(m))){
      row = as.numeric(m[r, ])
      expect_equal(row[length(row)], mock$n)
      expect_true(all(diff(row) > 0))
    }
  }
})

test_that("CS1: a higher confidence level gives a superset", {
  # the penalty matters here. with a small one PELT prunes so little that every level
  # returns the same set and this test would pass without proving anything, so use a big
  # penalty where the budget actually decides what gets in
  set.seed(1)
  data = c(rnorm(50), rnorm(50)+5, rnorm(50))
  mock = mock_pelt_checklist(data, 4*log(150))

  sig = function(level){
    core = changepoint:::cs_core(mock$op_cpts, mock$checklists, mock$no_op_cpts,
                                 mock$penalty, level=level)
    segs = changepoint:::cs_backtrack(mock$lastchangecpts, core$CS, mock$op_cpts)
    unlist(lapply(segs, function(m){
      if(is.null(m)) return(NULL)
      apply(rbind(m), 1, paste, collapse=",")
    }))
  }
  s90 = sig(0.90); s95 = sig(0.95); s99 = sig(0.99)
  expect_true(all(s90 %in% s95))
  expect_true(all(s95 %in% s99))
  # and the sets must really grow, otherwise the two checks above are just
  # "a set contains itself" and would pass even if the level were ignored
  expect_gt(length(unique(s99)), length(unique(s95)))
})

test_that("zero-changepoint case runs end to end", {
  set.seed(3)
  data = rnorm(50)
  mock = mock_pelt_checklist(data, 50) # huge penalty so no changes are found
  expect_equal(mock$op_cpts, 50)
  core = changepoint:::cs_core(mock$op_cpts, mock$checklists, mock$no_op_cpts,
                               mock$penalty, level=0.95)
  segs = changepoint:::cs_backtrack(mock$lastchangecpts, core$CS, mock$op_cpts)
  expect_equal(as.numeric(unique(rbind(segs[[1]]))[1,]), 50)
})
