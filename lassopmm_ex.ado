/*==================================================
project:       Examples in lassopmm.sthl
Author:        Andres Castaneda, Paul Corral, Leo Lucchetti
----------------------------------------------------
Creation Date:    14 Dec 2018 - 11:29:32
==================================================*/

/*==================================================
                        0: Program set up
==================================================*/
program define lassopmm_ex
version 14
set more off
	`0'
end

program Msg
	di as txt
	di as txt "-> " as res `"`0'"'
end

program Xeq
	di as txt
	di as txt `"-> "' as res _asis `"`0'"'
	`0'
end


program ex1
if ("`1'" == "trace") local trace "set trace on"
Msg preserve

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

`trace'
lassopmm `_y' `_x' [aw=weight], knn(5) add(5) psu(psu) seed(12388) uniqid(_numobs11)
mi estimate: mean price if samples==1 [aw=weight]

Msg restore 
end 

exit
/* End of do-file */

><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
