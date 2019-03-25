*! version 0.2 December27-2018
*! R. Andres Castaneda - acastaneda@worldbank.org
*! Paul Corral - pcorralrodas@worldbank.org  -- Corresponding author
*! Leonardo Lucchetti - llucchetti@worldbank.org
*! World Bank Group - Poverty and Equity Global Practice 
*! Equity policy lab

cap program drop lassopmm
program define lassopmm, eclass byable(recall)
	version 14

#delimit;
	syntax varlist(min=2 numeric fv) [if] [in] [aw],
		uniqid(varlist)
		add(integer)
	[
		SORTy
		knn(numlist max=1 int)
		lambda(real -1)
		numlambda(integer 100)
		numfolds(integer 10)
		psu(varlist max=1 numeric)
		seed(integer 12345)
		NOIsily
		EPSIlon(real 0.001)
	];
#delimit cr

set more off

qui{
*===============================================================================
// 1. HOUSE KEEPING
*===============================================================================

	//Mlong?
	local mlong `_dta[_mi_style]'
	
	//confirm mi registered
	if ("`_dta[_mi_marker]'"==""){
		dis as error "You must mi set your data before using this command, see help file"	
		error 119
		exit
	}
	
	//mark sample for use
	marksample touse
	
	// Mark for the other variables in the estimation
	foreach j of varlist `psu' `uniqid'{
		replace `touse' = 0 if missing(`j')
	}
	
	
	tokenize `varlist'
	// Local for dependent variable
	local depvar `1'
	
	// obtain the independent variables
	macro shift 
	local _my_x `*'
	
	if ("`_dta[_mi_ivars]'"!="`depvar'"){
		dis as error "You must mi register `depvar' as imputed before using this command"	
		error 198
		exit
	}	
	
	if ("`sorty'"!="" & "`knn'"!=""){
		//Only display warning!
		dis as error "Note that KNN option only works when the default option of PMM is used."
	}
	
	
	//Weights for estimation
	tempvar wi _psu _group
	if missing(`"`weight'"') generate byte `wi' = 1
	else generate `wi' `exp'
	
	//Check unique id
	isid `uniqid' 
	
	//PSU for bootstrapping
	if ("`psu'"==""){
		sort `uniqid'
		gen `_group' = 1
	}
	else{
		sort `psu' `uniqid'
		clonevar `_group' = `psu'		
	}
	
*===============================================================================
// 2. Get bootstrap permutation vectors and load data
*===============================================================================
	*---------- Preserve original data
	preserve 
	
	set seed `seed'	
	noi mata: r = lpmmset()         ;  /* 
         */ lpmm_getmyboot(r, "`_group'", "`touse'")
	
	
*===============================================================================
// 3. Run lasso and imputation on different bootstrapped data
*===============================================================================
	//Get numsim vectors of Y	
	forval i=1/`add' {		
		
		noi mata: lpmm_iteration(r, "`touse'", "`depvar' `_my_x' `wi'", `i')
		
		`noisily' lassoregress `depvar' `_my_x' if `touse'==1 [aw=`wi'],    /* 
			 */	                    numlambda(`numlambda') numfolds(`numfolds') /* 
			 */                     lambda(`lambda') epsilon(`epsilon')
				
		tempname _beta BB
		tempvar touse1 touse2
		
		mat `_beta' = e(b)
		
		gen `touse1' = e(sample)
		gen `touse2'= `touse1'==0 & missing(`depvar')
		
		
		local a=1
		local chosen
		foreach x of local _my_x {
			if (`_beta'[1,`a']!=0){
				local chosen `chosen' `x'
				mat `BB' = nullmat(`BB'),`_beta'[1,`a']
			}
			local ++a
		}
		mat `BB' = `BB',`_beta'[1,`a']
		
		foreach x of local chosen{
			replace `touse2' = 0 if missing(`x')
		}
			
			
		mata: lpmm_load(r, "`BB'", "`chosen'", "`touse1'", "`touse2'", "`depvar'", "`wi'")
		mata: y1 = lpmm_execute(r)
			
		restore, preserve 
	
	} //End of sim loop
	

*===============================================================================
// 3. Export simulations to mi set data
*===============================================================================

	restore   // restore original data.
	if ("`mlong'"=="mlong"){
		tempvar touse2 nn
		gen `touse2' = `touse'==0 & missing(`depvar')
		qui:count 
		local hu0 = r(N)
			preserve
				qui:keep if `touse2'==1
				tempfile _nnio
				qui:save `_nnio'
				qui:count
				local hu=r(N)
			restore
		local a=1
		while (`a'<=`add'){
			qui:append using `_nnio', gen(`nn')
			qui:replace _mi_m=`a' if `nn'==1
			local ++a
			drop `nn'			
		}
		qui: char _dta[_mi_M] `add'
		qui: char _dta[_mi_n] `hu'
		qui: char _dta[_mi_N] `hu0'		
		replace `touse2'=`touse2'==1 & _mi_m>0
		mata: st_store(.,st_varindex(tokens("`depvar'")),"`touse2'",y1)				
	}
	else{
		tempvar touse2 nn
		gen `touse2' = `touse'==0 & missing(`depvar')
		replace _mi_miss = 1 if `touse2'==1
		local a=1
		local toimp
		while (`a'<=`add'){
			gen double _`a'_`depvar' = `depvar' if `touse'==1
			
			local toimp `toimp' _`a'_`depvar'
			local ++a				
		}
		qui: char _dta[_mi_M] `add'
		mata: st_store(.,st_varindex(tokens("`toimp'")),"`touse2'",y1)				
	}	
}
end

*===============================================================================
// Annex: MATA functions
*===============================================================================
set matastrict on
cap mata drop lpmm*()

mata:

struct lpmm  {

	// -------- Input variables

	real matrix     x      // covariates for Complete data based on selected betas
	real colvector  w      // weights of Complete data based on selected betas
	real colvector  y      // dependent var with COMPLETE data	
	real matrix     x1     // covariates for Incomplete data
	real colvector  w1     // weights for Incomplete data
	real colvector  psu    // psu variable
	
	real matrix     b      // selected Betas from lasso (estimates)
	
	// input options
	
	real scalar sorty     // if sorty option selected
	real scalar mlong     // if mlong form selected
	real scalar knn       // No. of neighbors
	real scalar sim       // No. of similations
	
	//-----------derived
	real colvector yh1     // predicted y for COMPLETE data
	real colvector yh2     // predicted y for INconplete data
	real matrix    ind     // indicator of bootstrapingb
	real colvector myy     // order of Knn in Y
  real matrix    info    // panelsetup of psu
  real colvector index   // index of psu 
  real colvector indexm  // Subset index of psu 
  real scalar    dim     // dimension (rows) of indexm
	
	//------- Output
	real matrix    y1      // output matrix
	
}


struct lpmm scalar lpmmset() {
	struct lpmm scalar r
	
	 r.x   = J(0,0, .z)
	 r.x1  = J(0,0, .z)
	 r.y1  = J(0,0, .z)
	 r.b   = J(0,0, .z)
	 r.y   = J(0,1, .z)
	 r.w   = J(0,1, .z)
	 r.w1  = J(0,1, .z)
	
	if (st_local("sorty") != "") r.sorty = 1
	else                         r.sorty = 0
	
	if (st_local("mlong") == "mlong") r.mlong = 1
	else                              r.mlong = 0
	
	r.knn   = strtoreal(st_local("knn"))
  r.sim   = strtoreal(st_local("add"))
	
	return(r)
}



void lpmm_iteration(struct lpmm scalar r,
                    string scalar touse,
									  string scalar YXw,
									  real   scalar i) {
	
	real matrix myvar
	
	st_view(myvar=., .,YXw, touse)
	myvar[.,.] = myvar[r.ind[.,i],.]
	 
}


void lpmm_load(struct lpmm scalar r,
               string scalar BB, 
               string scalar chosen,
							 string scalar touse1,
							 string scalar touse2,
							 string scalar depvar,
							 string scalar wi) {
 
	
	r.b = st_matrix(BB)
	st_view(r.x   =., ., chosen , touse1)
	st_view(r.w   =., ., wi     , touse1)
	st_view(r.y   =., ., depvar , touse1)
	st_view(r.x1  =., ., chosen , touse2)
	st_view(r.w1  =., ., wi     , touse2)
	
	// Yhats
	
	r.yh1 = quadcross((r.x  , J(rows(r.x) , 1,1))', r.b')
	r.yh2 = quadcross((r.x1 , J(rows(r.x1), 1,1))', r.b')
	
}


//Function for bootstrap permutation vectors
void lpmm_getmyboot(struct lpmm scalar r, 
                    string scalar group,
										string scalar touse) {
		
	st_view(r.psu=., ., group, touse)
	real scalar rr

	r.info  = panelsetup(r.psu,1)
	r.index = runningsum(J(rows(r.psu),1,1))
	
	rr = rows(r.info)
	
	for(i=1; i<=rr; i++){
		m = r.info[i,1],1 \ r.info[i,2],1
		r.indexm = r.index[|m|]
		r.dim = rows(r.indexm)
		if (i==1) r.ind =         lpmm_sampleepsi(r)
		else      r.ind = r.ind \ lpmm_sampleepsi(r)
	}
}	

real matrix lpmm_sampleepsi( struct lpmm scalar r, 
														 | transmorphic f) {
	
	real scalar n, N
	real matrix sige2
	pointer(real matrix) scalar eps
	
	if (args() == 1) { // for lpmm_getmyboot
		n   = r.sim
		eps = &r.indexm
	}
	else {     // for lpmm_randomleo() and lpmm_Mpmm()
	  n   = f
		eps = &r.myy
	}
	
	
	sige2 = J(r.dim,n,0)
	
	N = rows(*eps)
	
	if (cols(*eps)==1)  {
		for(i=1; i<=n; i++) {	
			sige2[.,i] = (*eps)[ceil(N*runiform(r.dim,1)),1]
		}
	}
	
	else  {
		for(i=1; i<=n; i++) {
			sige2[.,i] = (*eps)[ceil(N*runiform(r.dim,1)),i]
		}
	}

	return(sige2)
		
}

//Function will return an index selection vector for PMM Y

 
real matrix lpmm_execute(struct lpmm scalar r)  {
	
	real colvector ytemp
	
	if (r.sorty == 1) ytemp = lpmm_randomleo(r)
	else              ytemp = lpmm_Mpmm(r)

	if (rows(r.y1) == 0) {
		r.y1 = ytemp
	}
	else {
		if (r.mlong == 1 ) {
			r.y1 = r.y1 \ ytemp
		}
		else {
			r.y1 = r.y1 , ytemp
		}
	
	}
	
	return(r.y1) // this could be eliminated if we call r.y1 instead of y1. 
	
}

real colvector lpmm_Mpmm (struct lpmm scalar r) {
	//Search distance to yh2
	
	real colvector mynn
	real scalar ry1, ry2
	
	ry2 = rows(r.yh2)
	ry1 = rows(r.yh1)  // what is this for?
	
	r.dim = 1
	
	if (r.knn > 1){
		for(i=1; i<=ry2; i++) {
			
			 r.myy = order(abs(r.yh1 :- r.yh2[i]),1)[|1\r.knn|]				
			 
			 if (i==1) mynn =        lpmm_sampleepsi(r, 1)
			 else      mynn = mynn \ lpmm_sampleepsi(r, 1) 
			
		}
	}
	else {
		for(i=1; i<=ry2; i++){
			if (i==1) mynn =        order(abs(r.yh1:-r.yh2[i]),1)[1]	
			else      mynn = mynn \ order(abs(r.yh1:-r.yh2[i]),1)[1]	
		}	
	}
	return(r.y[mynn])
}

//Function, selects a random set of observed y, of length rows unobserved. 
// sorts the set and assigns the value to the sorted xb from unobserved
real colvector lpmm_randomleo(struct lpmm scalar r) {
	
	real matrix tosort
	r.myy = r.y
	r.dim  = rows(r.yh2)
	//random sample of y
	tosort      = sort((runningsum(J(r.dim,1,1)),r.yh2),2)
	tosort[.,2] = sort(lpmm_sampleepsi(r, 1),1)
	_sort(tosort,1)
	
	return(tosort[.,2])
}	



end 




exit
/* End of do-file */

><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

Notes:
1.
2.
3.


Version Control: