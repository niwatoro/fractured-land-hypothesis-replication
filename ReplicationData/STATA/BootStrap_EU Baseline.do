clear



import delimited "mainfolder2/herfEU.csv"

drop v1
xpose, clear



* Initializing Empty Storage Matrix
mat observe = J(500,8,.)

* Generating Preliminary Estimates of h-index
*local a=1
*	foreach v of varlist v1-v500{
*		sum `v', meanonly
*		mat observe[`a',1]=r(mean)
*		loc ++a
*		}

mat list observe

* Generating Bootstrap Estimates of h-index
local a=1
	foreach v of varlist v1-v500{
		quietly bootstrap r(mean), reps(10000) size(30): summarize `v', detail
		mat observe[`a',1]=e(b)
		mat observe[`a',2]=e(b_bs)
		mat observe[`a',3]=e(bias)
		mat observe[`a',4]=e(se)
		mat observe[`a',5]=r(table)[3,1]
		mat observe[`a',6]=r(table)[4,1]
		mat observe[`a',7]=e(ci_bc)[1,1]
		mat observe[`a',8]=e(ci_bc)[2,1]
		mat list observe
		loc ++a
		}

mat list observe

* Naming Rows and Columns
local rownames
forvalues j = 1/500 {
    local rownames `rownames' `=0+`j''
}
matrix rownames observe = `rownames'
matrix list observe

matrix colnames observe = initial bs_est bs_bias bs_se bs_z bs_pvalue bs_cilower bs_ciupper
matrix list observe

*Export Results into Excel

putexcel set "mainfolder2/ResultsEU_10000.xlsx", sheet("ResultsEU") replace
putexcel A1=matrix(observe), names
