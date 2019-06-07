set more off
set matastrict off
clear 
run "C:\Users\WB378870\OneDrive - WBG\000.my_ados\lassopmm\postlassopmm.ado"
set seed 5678
sysuse auto, clear

gen psu = 7 if _n<74
replace psu = 6 if _n<60
replace psu = 5 if _n<50
replace psu = 4 if _n<40
replace psu = 3 if _n<30
replace psu = 2 if _n<20
replace psu = 1 if _n<10

gen _numobs = _n
preserve
	sample 10
	sum price
	replace price = .
		
	tempfile uno
	save `uno'
restore
append using `uno', gen(samples)

gen _numobs11 = _n

//Local with all candidate variables
local _x mpg headroom trunk weight length turn displacement gear_ratio foreign
//Local with dependent variable
local _y price 

mi set wide

mi register imputed price 


lassopmm `_y' `_x' [aw=weight], knn(1) add(1) psu(psu) seed(12388) uniqid(_numobs11) postlasso

replace price = _1_price if price==. 

drop _*_price  _numobs11 _mi_miss _1_price

