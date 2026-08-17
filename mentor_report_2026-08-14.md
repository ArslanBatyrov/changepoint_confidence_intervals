---

---

# Confidence Sets for `changepoint`: status report

Google S of CODE! , 15 August 2026

Questions, explanations and demonstration for the meeting! :)

## Headline

The feature is built end to end and working. A PELT fit can now produce a confidence set: a collection of entire plausible segmentations, capturing uncertainty in both the number and the location of the changepoints. It is available through two functions, `confidence_set()` and `confint()`. All five planned milestones are done.

While testing it I found a real limitation of the method that needs your decision. It is in section 5, and Owen, it lines up with a line you commented out in `PELT_CS_fns.R`.

## Verified numbers

|  |  |
|------------------------------------|------------------------------------|
| Commits since 27 June | **36** |
| Lines added | **\~1,025** across 17 files |
| Full package suite | **10,049 passing, 0 failing** |
| Confidence set tests | **289 assertions, 25 blocks, 2 files** |
| Effect on the normal PELT path | **none**, `conf.set=FALSE` gives byte identical output |

## 1. The C capture

The idea: normal PELT throws away the near miss candidate changepoints as it runs. Those near misses are exactly what a confidence set is built from. So the C job is to catch them before they are discarded.

`PELTC` gained four arguments (the registration count went 11 to 15, and `changepoint_init.c` was updated to match):

-   `conf_set`, the on/off switch. When it is 0 the new code does not run at all, which is why the normal path is unchanged.
-   `checklist_positions`, an output buffer that R allocates and C fills with the candidate positions, one changepoint per row.
-   `checklist_likes`, the matching output buffer holding each candidate's cost, read as pairs with the positions.
-   `cs_dim`, a length 2 size handshake, covered in section 2.

When `conf.set=TRUE`, we run PELT normally, then **run the forward pass a second time** and, at each optimal changepoint, save the list of surviving candidates into that changepoint's row.

Why a second pass: we do not know where the optimal changepoints are until the first run has finished and we trace them back. Saving every candidate list at every time point during the first pass would use a lot of memory. A second pass is cheaper: PELT is linear, so running it twice is still linear. From `src/PELT_one_func_minseglen.c:221`:

``` c
/*real capture, runs only if conf_set = TRUE. we only know the optimal cpt times after the
  backtrack, so we just replay the forward pass again and this time at every optimal cpt we
  save the surviving candidates into that cpts matrix row (after the prune, like the mock).
  costs are saved PURE, no penatly: pure = tmplike - numchangecpts[x]*pen, so no new array.*/
```

### The pure cost trick, which I would like checked

The confidence set math needs each candidate's cost *without* the penalty added. But the value PELT already has, `lastchangelike`, has the penalty baked in. Instead of adding a whole new array to C to track the penalty free cost, I worked out that I can recover it by subtracting the penalties back off. The output I got was accurate yet I am unsure if this is a good way to go:

```         
lastchangelike[x]  =  pure(x) + numchangecpts[x]*pen - pen

tmplike[i]  =  lastchangelike[cl] + cost(cl,tstar) + pen
            =  pure(cl) + ncc(cl)*pen - pen + cost + pen
            =  [pure(cl) + cost]  +  ncc(cl)*pen
            =  purepath           +  ncc(cl)*pen

  =>  purepath  =  tmplike[i] - numchangecpts[cl]*pen
```

The starting cases check out: `lastchangelike[0] = -pen` with `numchangecpts = 0` gives `pure = 0`, and the `minseglen` prefill gives `pure = cost(0,j)`. So the capture is a single line and needs no extra memory. From `src/PELT_one_func_minseglen.c:275`:

``` c
checklist_positions[i*rowcap + r] = checklist[i];
checklist_likes[i*rowcap + r]     = tmplike[i] - numchangecpts[checklist[i]]*(*pen);
```

## 2. The memory contract (handing in the matrix from R to c)

The problem: when C hands back a flat list of candidates, R cannot tell where one changepoint's group of candidates ends and the next begins. We discussed this in out chat, Rebecca suggested to hand in the matrix to begin with to solve this issue.

R sets up a matrix in advance, one row per optimal changepoint, and C fills row `r` with the candidates for changepoint `r`. The grouping is then built in: R just reads it row by row and hands each row to `cs_core`.

C has no matrix type, so under the hood the matrix is a flat block of memory that R and C agree to read column by column. Cell (row `r`, column `i`) lives at `i*rowcap + r`.

The old single number `cs_len` became `cs_dim`, a length 2 vector holding rows and columns. It travels both ways: R sends in the size it allocated, and C writes back the size it actually needed. The argument count stays at 15, so no registration change was needed:

``` c
int rowcap = cs_dim[0];   /*what R allocated, this is what we index by*/
int colcap = cs_dim[1];
int rowneed = ncpts;      /*one row per optimal cpt, n counts as one*/
int colneed = 0;          /*widest checklist seen, reported back*/
...
cs_dim[0] = rowneed;      /*true dims back so R can re-run exact if it undershot*/
cs_dim[1] = colneed;
```

One care point: the index uses `rowcap`, the number of rows R allocated, not `rowneed`, the number it turned out to need. R fixed the memory layout when it allocated, so using the wrong row count would scatter the data diagonally across the matrix.

On the R side, if the first guess was too small in either direction, we allocate the exact size C asked for and run once more. From `R/PELT_one_func_minseglen.R:63`:

``` r
if(needed[1] > rows_have || needed[2] > cols_have){
  answer2 = tryCatch(run_peltc(needed[1], needed[2]), error = function(e){
    warning("confidence-set buffer too large to allocate (", needed[1], " x ", needed[2],
            "); confidence set may be incomplete", call.=FALSE)
    NULL
  })
```

The `tryCatch` is my own idea instead of using the 50% of the RAM warning. But we all discussed it and approved.

## 3. The R core

This is the math that turns the captured candidates into the actual confidence set. It is ported from your `PELT_CS_fns.R` and made general, so the costs arrive already computed from C and nothing is ever recalculated from a per distribution formula. That way one core serves all 12 cost functions.

| Function | Location | What it does |
|------------------------|------------------------|------------------------|
| `cs_runmin` | `R/conf_set.R:214` | running minimum, your `runmin` renamed |
| `cs_no_op_cpts` | `R/conf_set.R:224` | rebuilds `no.op.cpts` from `lastchangecpts` alone |
| `cs_penalty` | `R/conf_set.R:236` | your budget formula, unchanged |
| `cs_core` | `R/conf_set.R:249` | picks the candidates that fall within budget |
| `cs_op_seg` | `R/conf_set.R:277` | the bunny hop back through the backpointers |
| `cs_backtrack` | `R/conf_set.R:289` | rebuilds whole segmentations from the picks |

`cs_no_op_cpts` rebuilds the changepoint counts from `lastchangecpts` on purpose, rather than reading C's `numchangecpts`, because the two use a different convention that is off by one.

**Validation.** I checked our output against yours on 4 datasets by 3 confidence levels, running `PELT_lv` then `PELT_CS` then `conf.set.segment`. The optimal changepoints, the counts, the budgets, the per row picks, and the full set of segmentations all matched. The raw cost numbers differ by a constant (the `t*log(2*pi*sigma^2)` term), but that constant cancels when you ask "is this candidate within budget", so the sets that come out are identical.

## 4. The public API

`confidence_set()` and `confint()` do the same thing under two names, sharing one cached result. The expensive part (re-running PELT and shaping the candidates) does not depend on the level, so it lives in `cs_ingredients()` at `R/conf_set.R:71` and is stored once under the key `"ingredients"`. PELT then runs only once, and asking for a second level is nearly free. The cache is thrown out if the fit changed:

``` r
if(is.null(ing) || !all(cpts(object) %in% ing$cpts) || !identical(ing$penalty, pen.value(object))){
  ing = cs_ingredients(object)
  existing[["ingredients"]] = ing
}
```

C's matrix rows become the checklists `cs_core` expects, with the empty NA padding stripped, at `R/conf_set.R:112`:

``` r
checklists = lapply(seq_len(nrow(pelt$checklist_positions)), function(i){
  keep = !is.na(pelt$checklist_positions[i, ])
  list(positions = as.numeric(pelt$checklist_positions[i, keep]),
       likes     = as.numeric(pelt$checklist_likes[i, keep]))
})
```

It covers Normal (mean, variance, mean and variance), Exponential, Gamma and Poisson, with or without MBIC.

### Live output

```         
set.seed(1); x = c(rnorm(100), rnorm(100)+1.5, rnorm(100)-1)   # true changes at 100, 200

Confidence set (level = 0.95)
Optimal changepoints:
  100 203 300

Alternative segmentations (13), extra cost vs optimal (budget 13.93):
  100 200 300       +0.01     <- the TRUE changepoints
  101 203 300       +0.29
  103 203 300       +1.78
  100 202 300       +1.80
  ...
  100 200 203 300   +13.75    <- a DIFFERENT NUMBER of changepoints
```

Three things stand out. PELT's single best answer said 203, but the true segmentation is in the set at a gap of only +0.01, so statistically it is a dead heat (no difference literally, no rational justification of picking one over another). The set includes answers with a different *number* of changepoints, which is property CS6. And CS1 holds: 95% gives 14 members, 99% gives 26, and the 95% set sits inside the 99% set.

## 5. The finding: PELT's pruning cuts the confidence set short (MAIN ISSUE I AM THINKING ABOUT)

I tested whether the set keeps growing as the confidence level rises. It does not at the higher levels.

```         
pen = 2*log(300) = 11.41

level        budget    set size    max gap
0.95         13.93     14          13.83
0.99         17.15     26          17.00
0.999        21.75     26          17.00
0.99999      30.96     26          17.00
0.9999999    40.17     26          17.00
```

The budget quadruples, and the set does not move past 26.

### Why, in our C

We save the candidate list *after* the prune step, so we only ever see the candidates that survived pruning. From `src/PELT_one_func_minseglen.c:262`:

``` c
nchecktmp = 0;
for(i = 0; i < nchecklist; i++){
  if(tmplike[i] <= (lastchangelike[tstar]+*pen)){   /* the prune: the slack is exactly pen */
    checklist[nchecktmp] = checklist[i];
    tmplike[nchecktmp]   = tmplike[i];
    nchecktmp += 1;
  }
}
nchecklist = nchecktmp;                              /* the losers are overwritten and gone */
if(tstar == cptsout[ncpts-1-r]){
  /* our snapshot happens BELOW the prune, so it can only see survivors */
```

A candidate is kept by being copied to the front of the list. A candidate is dropped by simply not being copied: the list shrinks past it and the next write lands on top. There is no code that keeps the dropped ones, so nothing later can get them back. Once the budget is bigger than the pruning slack, the set is asking for candidates that no longer exist thus it stips expanding!

### The same thing happens in the prototype

This is not something the port introduced. In `PELT_CS_fns.R`, line 62 prunes and line 65 saves the CS costs using the *same* condition, so the prototype also keeps only survivors:

``` r
checklist = tmpt[tmplike <= lastchangeF[tstar+1]+penalty] # prunes the checklist.

lastchangelike[tstar+1] = tmplike_CS[which.min(tmplike)] # optimal pure likelihood
checklist.t[[tstar+1]] = list(checklist=checklist,
                              likes=tmplike_CS[tmplike <= lastchangeF[tstar+1]+penalty])
```

`tmplike <= lastchangeF[tstar+1]+penalty` is the same test in both codebases. That is why our output matches yours exactly: we are both saving the same survivors. Matching your output proves the port is faithful, but it does not prove either set is complete. So, the question arises, how best shall we solve this issue and shall I solve this issue if it is int he prototype itself or just port the prototype faithfully?

### Owen, this is your line 106

From `PELT_CS_fns.R:96` onwards, with the commented line at 106:

``` r
  if(is.null(CS.penalty)){
    if(length(cpts)==2){ # 0 cpts detected
      nj.maxi = n/20 # s.t. 2nj=n/10
    } else if(length(cpts)>2){
      nj.maxi = min(max(runmin(diff(cpts),2))
                    , n/4 - n*5/24)
    }
    CS.penalty = 3/2 + 2*log(2*nj.maxi)-2*log(alpha)
  }

  # CS.penalty[which(CS.penalty>PELT.lv.results$penalty)] = PELT.lv.results$penalty
  # if the CS threshold is bigger than beta, just set it equal to beta
```

That commented line caps the budget at the penalty, which is exactly the point where the cutting short begins. Turning it on would make the behaviour visible in the code (the set stops growing, and you can see why) instead of silent. If the intended rule is that the budget should never go above beta, then our saturation is correct and just needs documenting. If not, there is a real gap here that we misght need to work on.

### Why the freeze is at 17.00, not at 11.41

The prune is checked at *every* time step, against the best path so far, with the penalty included. The gap we print is measured *once*, at the changepoint we are conditioning on, with the penalty added back in. Because they are measured at different moments and on a slightly different basis, a candidate can stay within 11.41 at every step yet end up 17 behind by the time we measure it. There is no clean formula, so **the point where it freezes depends on the data.**

### A subtler worry: even 0.95 may be affected (This is the reason why I worry about issue the most, as it is so elusive, you just don't know its full impact)

As soon as the budget goes above the penalty, a gap opens. On this data it is `[11.41, 13.93]` at the 0.95 level. A candidate could drop 12 behind at some intermediate step, enough to get pruned, yet end up only 12 behind overall, which is *inside* the 0.95 budget. So even the 0.95 set may be missing a few candidates, not just the very high levels.

I cannot measure this yet, because everything I can validate against (`PELT_CS_fns.R`) prunes in the same way.

My proposed solution (but will we be able to make it withing the time we have is another question) :

### First step: measure how bad it is, with no C changes

Before touching any C code, I want to measure how severe this actually is. If pruning only cuts off a few candidates at extreme levels, we can document it and leave it. If it eats into the set at ordinary levels like 0.95, it needs a fix. Right now I do not know which, because everything I can validate against prunes the same way.

Segment Neighbourhood is the tool for this. It is the other exact method in the package, and it gives a complete confidence set because it never prunes. Running it next to PELT on the same data tells me exactly what pruning threw away. It works as a ruler for PELT for these reasons:

-   Same aim. Both minimise the same cost plus penalty, so they agree on the single best segmentation.
-   Both are exact. They find the true optimum, so the optimal cost is the same number in both.
-   Because the optimum is the same, any alternative segmentation's gap from it is also the same number in both. A gap is a property of the segmentation and the data, not of the algorithm that found it.
-   SegNeigh never prunes, so it holds every candidate, including the ones PELT deleted, and can put a gap on each.

So the difference between SegNeigh's set and PELT's set, at the same level on the same data, is exactly the amount PELT's pruning cut off. That turns the worry into a real number, and it needs no change to the C code, just running both and comparing.

### The fix, only if the measurement says we need it

If it turns out to be severe: pruning less is always safe for the answer. PELT's own theorem says a pruned candidate can never turn out to be optimal later, so keeping extra candidates cannot change the result, it only costs time. So when `conf.set=TRUE` we could widen the slack from `pen` to `max(pen, budget)`, giving the same changepoints and a complete set, with the extra cost paid only when someone asks for a confidence set.

The catch: the budget depends on the confidence level, and the level currently lives in R. C never sees it. So this would reverse the decision from 7 July. That is why I want the measurement first, before any C change.

## 6. Questions

**Q1, the pruning cut off (section 5).** My plan is to measure the severity first with Segment Neighbourhood, which needs no C change, and then decide: if it is small we document it, if it is large we widen the capture prune to `max(pen, budget)`, which means the level has to enter C after all. Does that order sound right? Owen: was line 106 removed on purpose?

**Q2, does it hold beyond Normal? ALREADY ANSWERED\< REBECCA WILL HAVE A LOOK AT IT AFTER THE PROJECT.** The budget formula was derived for the Normal case. Exponential, Gamma and Poisson are switched on as you asked, with the level flagged as provisional in NEWS. What would it take to make it solid, a simulation study, or is there exponential family theory to lean on? Same question for MBIC, where the extra length dependent term sits inside the cost rather than the penalty.

**Q3, the second pass.** Running PELT's forward pass twice, versus storing every candidate list during the first pass. Is running it twice the right trade?

**Q4, is answered**

**Q5, the Segment Neighbourhood code. ALREADY ANSWERED, I MUST EMAIL REBECCA AS SOON AS I AM DONE WITH PELT FULLY, yet to be done with PELT.** Owen, you mentioned that both the PELT and the Segment Neighbourhood code live either in the `changepoint` repo or in your own repo. I have looked, and I can only find `PELT_CS_fns.R`, which is PELT only, plus the `segneigh.*` detection functions in the package, which have no confidence set logic. **Could you send me the link to the Segment Neighbourhood confidence set code so I can port it?** With it, the next milestone is a port and validate job like the one I already did for PELT, roughly one week, and I would match your structure. Without it, it is new design work, roughly two weeks.

------------------------------------------------------------------------

## 7. What is next

### Ready now, written but not yet committed

Each segmentation now carries its **gap**, meaning how much more it costs than the best answer, and the set carries the **budget**. The printout sorts by gap, so the most plausible alternatives come first. Tests and documentation are written and the suite is green. Ready to commit.

### Depending on the answers:

### **1. Whatever Q1 decides**, either document the cut off or build the widened prune.

**2. Segment Neighbourhood.** It is worth building for its own sake, but it is also the tool for measuring Q1. SegNeigh is exact like PELT, so the two agree on the best answer and their gaps are on the same scale, but SegNeigh **never prunes**, so its confidence set is complete. The difference between the two sets is exactly the amount PELT's pruning cuts off, measured directly.

It is pure R. The candidate list I need is already computed and thrown away one line later, for example in `segneigh.mean.norm` at `R/multiple.norm.R:256`:

``` r
like = like.Q[q-1,v] + all.seg[v+1,j]    # the full candidate list, nothing pruned
like.Q[q,j] = max(like, na.rm=TRUE)      # only the winner is kept
cp[q,j] = which(like==max(like,na.rm=TRUE))[1]+(q-2)
```

The main differences from PELT:

-   The backpointers are 2D (`cp[q,j]`), so `cs_op_seg` and `cs_backtrack` need a 2D version. From `(q, j)` you hop to `(q-1, cp[q,j])`.
-   Uncertainty in the *number* of changepoints comes for free, because SegNeigh already computes the best segmentation for every number of segments and keeps the whole `criterion` vector.
-   The costs are on a maximised log likelihood scale, so they need converting to line up with the budget.
-   The six likelihood based `segneigh.*` functions are not quite uniform. The `v` range and the `cp[q,j]` offset differ between `mean.norm` and the rest, and three of them return `like.Q[,n]` while the others return `-2*like.Q[,n]`.

**3. BinSeg is not a candidate.** It is approximate, so its best answer is not the same best answer and the gaps would not line up. Of the three methods in the package, only SegNeigh can serve as ground truth for PELT.
