# tests for the PUBLIC (user-facing) side: confidence_set() and confint().
# basically, if a user can type it, it should be tested somewhere in here.
# every bug I found by hand poking on 2026-08-07 is pinned down here so it cannot cuase any issues

fit_of = function(seed = 1){
  set.seed(seed); cpt.mean(c(rnorm(50,0,1), rnorm(50,5,1), rnorm(50,0,1)), method = "PELT")
}
set_of  = function(obj, lv = 0.95) conf.set(confidence_set(obj, lv))[[as.character(lv)]]
rows_of = function(cs) unlist(lapply(cs$segmentations, function(g) apply(g, 1, paste, collapse=" ")))

# spies on PELT and counts how many times it really ran. no other way to prove the cache is
# doing its job, the answers look the same either way, only the count gives it away.
# callers pass an env because <<- would not reach back into the test_that frame
count_pelt = function(expr){
  calls = 0
  suppressMessages(trace(PELT, tracer = function() calls <<- calls + 1, print = FALSE,
                         where = asNamespace("changepoint")))
  on.exit(suppressMessages(untrace(PELT, where = asNamespace("changepoint"))))
  force(expr); calls
}

test_that("confidence_set stores the set in the slot and hands the object back", {
  fit = fit_of()
  expect_identical(conf.set(fit), list()) # a fresh fit still looks like it always did

  out = confidence_set(fit, 0.95)
  expect_s4_class(out, "cpt")            # the fit comes back, not the set
  # R hands functions a copy, so the users own fit must come out untouched. this is the
  # whole reason they have to write fit = confidence_set(fit, ...) and not just call it
  expect_identical(conf.set(fit), list())

  cs = conf.set(out)[["0.95"]]
  expect_s3_class(cs, "cpt.confset")
  expect_equal(names(cs), c("level","cpts","segmentations","gaps","budget"))
  expect_equal(as.numeric(cs$cpts), as.numeric(c(cpts(fit), length(data.set(fit)))))
})

test_that("confint is an identical twin and levels accumulate from either name", {
  fit = fit_of()
  expect_identical(confidence_set(fit, level=0.95), confint(fit, level=0.95))

  fit = confidence_set(fit, 0.95); fit = confint(fit, level=0.99); fit = confidence_set(fit, 0.9)
  expect_setequal(setdiff(names(conf.set(fit)), "ingredients"), c("0.95","0.99","0.9"))
})

test_that("the slot cache runs PELT once, and a stale cache is not served", {
  e = new.env(); e$fit = fit_of()
  n = count_pelt({
    e$fit = confidence_set(e$fit, 0.9); e$fit = confidence_set(e$fit, 0.95)
    e$fit = confint(e$fit, level = 0.99)
  })
  expect_equal(n, 1)
  expect_true("ingredients" %in% names(conf.set(e$fit)))

  good = conf.set(e$fit)[["0.95"]]
  broken = conf.set(e$fit); broken$ingredients$cpts = c(7, 150) # vandalise the cache
  conf.set(e$fit) = broken
  expect_equal(count_pelt({ e$fit = confidence_set(e$fit, 0.95) }), 1) # it redid the work
  expect_identical(conf.set(e$fit)[["0.95"]], good)                    # and got it right
})

test_that("the optimal is always in the set and a higher level gives a superset (CS1)", {
  fit = fit_of(); prev = NULL
  for(lv in c(0.9, 0.95, 0.99)){
    cs = set_of(fit, lv)
    rows = rows_of(cs)
    expect_true(paste(cs$cpts, collapse=" ") %in% rows)
    # counting is not enough, it has to be the SAME segmentations plus more. anything
    # plausible at 90 must still be plausible at 99, you cannot lose one by asking for
    # more confidence
    if(!is.null(prev)) expect_true(all(prev %in% rows))
    prev = rows
  }
  # and it has to actually grow, else the check above is just "a set contains itself"
  expect_gt(length(rows_of(set_of(fit, 0.99))), length(rows_of(set_of(fit, 0.9))))
})

test_that("segmentations are tidy: no holes, no rownames, no duplicates, named by width", {
  cs = set_of(fit_of(), 0.99)
  expect_false(any(vapply(cs$segmentations, is.null, logical(1))))
  for(k in names(cs$segmentations)){
    g = cs$segmentations[[k]]
    expect_true(is.matrix(g))
    expect_null(dimnames(g))
    expect_equal(nrow(g), nrow(unique(g)))
    expect_equal(ncol(g), as.integer(k))
  }
})

test_that("each segmentation carries a gap: aligned, in [0, budget], optimal is 0", {
  cs = set_of(fit_of(), 0.99)
  for(k in names(cs$segmentations)){
    # one gap per segmentation row, and they line up
    expect_equal(length(cs$gaps[[k]]), nrow(cs$segmentations[[k]]))
  }
  g = unlist(cs$gaps)
  expect_true(all(g >= 0))
  expect_true(all(g <= cs$budget + 1e-9))
  # the optimal segmentation is the only gap-0 row (dedup keeps it once)
  opt_k = as.character(length(cs$cpts))
  expect_equal(cs$gaps[[opt_k]][1], 0)
  expect_equal(sum(g == 0), 1)
})

test_that("printing shows the gap column, sorted, with the budget in the header", {
  cs = set_of(fit_of(), 0.99)
  txt = capture.output(print(cs))
  expect_true(any(grepl("extra cost vs optimal \\(budget", txt)))
  expect_true(any(grepl("\\+[0-9]+\\.[0-9]{2}", txt))) # a +N.NN gap appears
  # alternatives are sorted smallest gap first: pull the printed gaps and check ascending
  printed = as.numeric(sub(".*\\+([0-9.]+)\\s*$", "\\1", grep("\\+[0-9]", txt, value=TRUE)))
  expect_false(is.unsorted(printed))
})

test_that("every PELT distribution and the MBIC keys run through the public path", {
  set.seed(2); xe = c(rexp(60,1), rexp(60,5))
  set.seed(3); xg = c(rgamma(60,shape=3,rate=1), rgamma(60,shape=3,rate=4))
  set.seed(4); xp = c(rpois(60,2), rpois(60,9))
  set.seed(6); xv = c(rnorm(60,0,1), rnorm(60,0,4))
  set.seed(1); xn = c(rnorm(50,0,1), rnorm(50,5,1), rnorm(50,0,1))
  fits = list(exp   = cpt.meanvar(xe, method="PELT", test.stat="Exponential"),
              gamma = cpt.meanvar(xg, method="PELT", test.stat="Gamma", shape=3),
              pois  = cpt.meanvar(xp, method="PELT", test.stat="Poisson"),
              var   = cpt.var(xv, method="PELT"),
              mbic  = cpt.mean(xn, method="PELT", penalty="MBIC"))
  for(nm in names(fits)){
    cs = set_of(fits[[nm]], 0.95)
    expect_s3_class(cs, "cpt.confset")
    expect_true(paste(cs$cpts, collapse=" ") %in% rows_of(cs)) # optimal is in every one
  }
})

test_that("edge fits work: no changepoints, and minseglen is respected", {
  set.seed(5); flat = cpt.mean(rnorm(100), method="PELT")
  expect_equal(ncpts(flat), 0)
  expect_equal(as.numeric(set_of(flat)$cpts), 100) # just n

  # 0.99 on purpose. at 0.95 the set is just the optimal, and PELT already guarantees THAT
  # respects minseglen, so we would be testing nothing. the alternatives are the risky ones
  set.seed(1); x = c(rnorm(50,0,1), rnorm(50,5,1), rnorm(50,0,1))
  cs = set_of(cpt.mean(x, method="PELT", minseglen=15), 0.99)
  expect_gt(length(rows_of(cs)), 1) # if this ever drops to 1 the loop below is pointless
  for(g in cs$segmentations){
    expect_true(all(apply(g, 1, function(r) all(diff(c(0, r)) >= 15))))
  }
})

test_that("the error guards fire, and cpt.range names all three methods", {
  set.seed(1); x = c(rnorm(50,0,1), rnorm(50,5,1), rnorm(50,0,1))
  for(bad in list(cpt.mean(x, method="BinSeg"),
                  suppressWarnings(cpt.mean(x, method="SegNeigh", Q=3, penalty="SIC")))){
    expect_error(confidence_set(bad), "BinSeg, SegNeigh or CROPS")
  }
  expect_error(confidence_set("not a cpt"), "must be a cpt object") # our message, not S4's

  fit = fit_of()
  for(lv in list(0, 1, 1.2, NA, c(0.9,0.95), "x")){
    expect_error(confidence_set(fit, lv), "strictly between 0 and 1")
  }
})

test_that("confint reads a positional 2nd argument as the level, not as parm", {
  fit = fit_of()
  # this used to quietly put 0.99 in parm and compute 0.95 instead
  expect_equal(setdiff(names(conf.set(confint(fit, 0.99))), "ingredients"), "0.99")
  # a character parm (lm style) means nothing here but must not error, level stays default
  expect_equal(setdiff(names(conf.set(confint(fit, parm="sigma"))), "ingredients"), "0.95")
})

test_that("the getter classes the list for printing and both setters strip it again", {
  fit = confidence_set(fit_of(), 0.95)
  got = conf.set(fit)
  expect_s3_class(got, "cpt.confsets")
  reg = new("cpt.reg")
  expect_silent(conf.set(fit) <- got) # cpt
  expect_silent(conf.set(reg) <- got) # cpt.reg has to behave the same
  expect_false(inherits(fit@conf.set, "cpt.confsets")) # the slot keeps a plain list
})

test_that("printing shows the optimal first, caps softly and hides the cache", {
  fit = confint(confidence_set(fit_of(), 0.99), level = 0.95)
  cs  = conf.set(fit)[["0.99"]]

  txt = capture.output(print(cs))
  expect_match(txt[1], "Confidence set \\(level = 0.99\\)")
  expect_true(any(grepl("Optimal changepoints", txt)))
  expect_true(any(grepl(paste(cs$cpts, collapse=" "), txt, fixed=TRUE)))

  full = length(capture.output(print(cs, n = Inf)))
  expect_lt(length(capture.output(print(cs, n = 2))), full)
  expect_true(any(grepl("print\\(x, n = Inf\\)", capture.output(print(cs, n = 2)))))
  old = getOption("changepoint.confset.n"); on.exit(options(changepoint.confset.n = old))
  options(changepoint.confset.n = 1)
  expect_lt(length(capture.output(print(cs))), full) # the option moves the default

  slot_txt = capture.output(print(conf.set(fit)))
  expect_match(slot_txt[1], "0.95, 0.99")                        # sorted, but 0.99 went in first
  expect_false(any(grepl("checklists|lastchangecpts", slot_txt))) # cache stays hidden
})
