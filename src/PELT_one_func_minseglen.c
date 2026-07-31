#include <R.h> 
#include <Rmath.h>
#include <Rinternals.h> // RK addition
#include <R_ext/RS.h>  // RK addition
#include <R_ext/Lapack.h> // RK addition
#include <R_ext/BLAS.h> // RK addition
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
//#include "cost_general_functions.c"
#define SWAP(a,b)   { int t; t=a; a=b; b=t; }  // Macro for swapping

//static int *lastchangecpts;
//static double *lastchangelike;
static int *checklist;
static double *tmplike;
static int *tmpt;

void FreePELT(int* error)
{
  if(*error==0){
    //	free((void *)lastchangecpts);
    // free((void *)lastchangelike);
    free((void *)checklist);
    free((void *)tmplike);
    free((void *)tmpt);
  }
}

void PELTC(char** cost_func,
	  double* sumstat, // Summary statistic for the time serie
	  int* n, // Length of the time series
	  double* pen, // Penalty used to decide if a changepoint is significant
	  int* cptsout, // Vector of identified changepoint locations
	  int* error, // 0 by default, nonzero indicates error in code
	  double* shape, // only used when cost_func is the gamma likelihood 
	  int* minseglen, //minimum segment length
	  double* lastchangelike, // stores likelihood up to that time using
	                          //optimal changepoint locations up to that time
	  int* lastchangecpts, // stores last changepoint locations
	  int* numchangecpts,
	  int* conf_set, /*on/off switch*/
	  int* checklist_positions, // Output: near-miss candidate positions, an R matrix passed as
	                            // flat column-major memory (space allocated by R)
	  double* checklist_likes, // Output: paired costs, same matrix shape as positions
	  int* cs_dim //length 2, both input and output; arrives as c(rows, cols) capacity of the
	              //matrix R allocated (rows = optimal cpts, cols = candidate slots), and C
	              //writes back the dims it actually needed so R can re-run once with exact size
	  )
{
  // R code does know.mean and fills mu if necessary

  double (*costfunction)(double, double, double, int, double);
  double mll_var(double, double, double, int, double); 
  double mll_mean(double, double, double, int, double); 
  double mll_meanvar(double, double, double, int, double); 
  double mll_meanvar_exp(double, double, double, int, double); 
  double mll_meanvar_gamma(double, double, double, int, double); 
  double mll_meanvar_poisson(double, double, double, int, double); 
  double mbic_var(double, double, double, int, double); 
  double mbic_mean(double, double, double, int, double); 
  double mbic_meanvar(double, double, double, int, double); 
  double mbic_meanvar_exp(double, double, double, int, double); 
  double mbic_meanvar_gamma(double, double, double, int, double); 
  double mbic_meanvar_poisson(double, double, double, int, double); 

  if (strcmp(*cost_func,"var.norm")==0){
    costfunction = &mll_var;
  }
  else if (strcmp(*cost_func,"mean.norm")==0){
    costfunction = &mll_mean;
  }  
  else if (strcmp(*cost_func,"meanvar.norm")==0){
    costfunction = &mll_meanvar;
  }
  else if (strcmp(*cost_func,"meanvar.exp")==0){
    costfunction = &mll_meanvar_exp;
  }
  else if (strcmp(*cost_func,"meanvar.gamma")==0){
    costfunction = &mll_meanvar_gamma;
  }
  else if (strcmp(*cost_func,"meanvar.poisson")==0){
    costfunction = &mll_meanvar_poisson;
  }
  else if (strcmp(*cost_func,"mean.norm.mbic")==0){
    costfunction = &mbic_mean;
  }
  else if (strcmp(*cost_func,"var.norm.mbic")==0){
    costfunction = &mbic_var;
  }
  else if (strcmp(*cost_func,"meanvar.norm.mbic")==0){
    costfunction = &mbic_meanvar;
  }
  else if (strcmp(*cost_func,"meanvar.exp.mbic")==0){
    costfunction = &mbic_meanvar_exp;
  }
  else if (strcmp(*cost_func,"meanvar.gamma.mbic")==0){
    costfunction = &mbic_meanvar_gamma;
  }
  else if (strcmp(*cost_func,"meanvar.poisson.mbic")==0){
    costfunction = &mbic_meanvar_poisson;
  } 
  
  
  //  int *lastchangecpts;
  //  lastchangecpts = (int *)calloc(*n+1,sizeof(int));
  //  if (lastchangecpts==NULL)   {
  //    *error = 1;
  //    goto err1;
  //  }
  
  //double lastchangelike[*n]; /* stores likelihood up to that time using optimal changepoint locations up to that time */
  //    double *lastchangelike;
  //  lastchangelike = (double *)calloc(*n+1,sizeof(double));
  //  if (lastchangelike==NULL)   {
  //    *error = 2;
  //    goto err2;
  // }
  
  //int checklist[*n];
  int *checklist;
  checklist = (int *)calloc(*n+1,sizeof(int));
  if (checklist==NULL)   {
    *error = 1;
    goto err1;
  }
  

  
  int nchecklist;
  double minout;

  //double tmplike[*n];
  double *tmplike;
  tmplike = (double *)calloc(*n+1,sizeof(double));
  if (tmplike==NULL)   {
    *error = 2;
    goto err2;
  }
  
  
  //int tmpt[*n];
  int *tmpt;
  tmpt = (int *)calloc(*n+1,sizeof(int));
  if (tmpt==NULL)   {
    *error = 3;
    goto err3;
  }
  
  int tstar,i,whichout,nchecktmp;
  

  void min_which(double*, int, double*, int*);
  
  lastchangelike[0]= -*pen;
  lastchangecpts[0]=0; 
  numchangecpts[0]=0;
  // lastchangelike[1]=costfunction(*(sumstat+1),*(sumstat + *n + 1 + 1),*(sumstat + *n + *n + 2 + 1),1, *shape);
  // lastchangecpts[1]=0; lastchangecpts[*n+1]=1;
  
  int j; 
  
  for(j=*minseglen;j<(2*(*minseglen));j++){
    lastchangelike[j] = costfunction(*(sumstat+j),*(sumstat + *n + 1 + j),*(sumstat + *n + *n + 2 + j),j, *shape); 
    // lastchangelike[j] = mll_mean(n, sumstat, j, 0, j, *shape);
  }
  
  
  for(j=*minseglen;j<(2*(*minseglen));j++){ 
    lastchangecpts[j] = 0;
  }
  
  for(j=*minseglen;j<(2*(*minseglen));j++){ 
    numchangecpts[j] =1;
  }
  
  
  nchecklist=2;
  checklist[0]=0;
  checklist[1]=*minseglen;
  
  for(tstar=2*(*minseglen);tstar<(*n+1);tstar++){
    R_CheckUserInterrupt(); /* checks if user has interrupted the R session and quits if true */
    
    if ((lastchangelike[tstar]) == 0){ 
      for(i=0;i<(nchecklist);i++){
        tmplike[i]=lastchangelike[checklist[i]] + costfunction(*(sumstat+tstar)-*(sumstat+checklist[i]),*(sumstat + *n + 1 +tstar)-*(sumstat + *n + 1 +checklist[i]),*(sumstat + *n + *n + 2 +tstar)-*(sumstat + *n + *n + 2 +checklist[i]), tstar-checklist[i], *shape)+*pen;
      }
      min_which(tmplike,nchecklist,&minout,&whichout); /*updates minout and whichout with min and which element */
      lastchangelike[tstar]=minout;
      lastchangecpts[tstar]=checklist[whichout]; 
      numchangecpts[tstar]=numchangecpts[lastchangecpts[tstar]]+1;
      /* Update checklist for next iteration, first element is next tau */
      nchecktmp=0;
      for(i=0;i<nchecklist;i++){
	if(tmplike[i]<= (lastchangelike[tstar]+*pen)){
	  *(checklist+nchecktmp)=checklist[i];
	  nchecktmp+=1;
	}
      }
      nchecklist = nchecktmp;
    }
    
    
    *(checklist+nchecklist)=tstar-(*minseglen-1);// atleast 1 obs per seg
    nchecklist+=1;
    /*  nchecklist=nchecktmp;*/
  
  } // end taustar
  
  // put final set of changepoints together
  int ncpts=0;
  int last=*n;
  while(last!=0){
    *(cptsout + ncpts) = last; 
    last=lastchangecpts[last];
    ncpts+=1;
  }
  /*real capture, runs only if conf_set = TRUE. we only know the optimal cpt times after the
    backtrack, so we just replay the forward pass again and this time at every optimal cpt we
    save the surviving candidates into that cpts matrix row (after the prune, like the mock).
    costs are saved PURE, no penatly: pure = tmplike - numchangecpts[x]*pen, so no new array.
    matrix is column-major, cell (row r, col i) is at [i*rowcap + r], rowcap is the ALLOCATED
    rows not the needed ones. empty cells stay NA. row r = r-th smallest optimal cpt.*/
  if(*conf_set == 1){
    int rowcap = cs_dim[0]; /*matrix dims R allocated, our capacity*/
    int colcap = cs_dim[1];
    int rowneed = ncpts; /*one row per optimal cpt, n counts as one*/
    int colneed = 0; /*widest checklist seen, reported back so R can re-run exact*/
    int r = 0; /*next row = next optimal cpt (ascending) waiting for its capture*/
    int target;
    for(i = 0; i < rowcap*colcap; i++){
      checklist_positions[i] = NA_INTEGER;
      checklist_likes[i] = NA_REAL;
    }
    /*optimal cpts before the main loop starts (t < 2*minseglen) never had a checklist,
      the only possible candidate is 0 and its pure cost is lastchangelike itself
      (one segment, prefill added no penalty)*/
    while(r < ncpts && cptsout[ncpts-1-r] < 2*(*minseglen)){
      target = cptsout[ncpts-1-r];
      if(r < rowcap && colcap > 0){
        checklist_positions[r] = 0; /*cell (r,0): col 0 so index is just r*/
        checklist_likes[r] = lastchangelike[target];
      }
      if(colneed < 1){ colneed = 1; }
      r++;
    }
    /*replay: reset the checklist to its starting state and walk tstar forward again
      with the same cost calls and the same prune rule as the first pass*/
    nchecklist = 2;
    checklist[0] = 0;
    checklist[1] = *minseglen;
    for(tstar = 2*(*minseglen); tstar < (*n+1) && r < ncpts; tstar++){
      R_CheckUserInterrupt();
      for(i = 0; i < nchecklist; i++){
        tmplike[i] = lastchangelike[checklist[i]] + costfunction(*(sumstat+tstar)-*(sumstat+checklist[i]),*(sumstat + *n + 1 +tstar)-*(sumstat + *n + 1 +checklist[i]),*(sumstat + *n + *n + 2 +tstar)-*(sumstat + *n + *n + 2 +checklist[i]), tstar-checklist[i], *shape)+*pen;
      }
      /*prune exactly like the first pass, but compact tmplike alongside checklist
        so position/cost pairs stay aligned for the capture below*/
      nchecktmp = 0;
      for(i = 0; i < nchecklist; i++){
        if(tmplike[i] <= (lastchangelike[tstar]+*pen)){
          checklist[nchecktmp] = checklist[i];
          tmplike[nchecktmp] = tmplike[i];
          nchecktmp += 1;
        }
      }
      nchecklist = nchecktmp;
      if(tstar == cptsout[ncpts-1-r]){
        /*this tstar is an optimal cpt: dump the surviving candidates into row r*/
        for(i = 0; i < nchecklist; i++){
          if(r < rowcap && i < colcap){
            checklist_positions[i*rowcap + r] = checklist[i];
            checklist_likes[i*rowcap + r] = tmplike[i] - numchangecpts[checklist[i]]*(*pen);
          }
        }
        if(nchecklist > colneed){ colneed = nchecklist; }
        r++;
      }
      *(checklist+nchecklist) = tstar-(*minseglen-1);
      nchecklist += 1;
    }
    cs_dim[0] = rowneed; /*reporting the true dims back so R can compare with its guess*/
    cs_dim[1] = colneed;
  }

  free(tmpt);
 err3:  free(tmplike);
 err2:  free(checklist);
  // err3:  free(lastchangelike);
  // err2:  free(lastchangecpts);
 err1:  return;
}


