*Importing the results from bootstrapping CN simulation estimates of h-index
clear
import excel "mainfolder2/ResultsCN_10000.xlsx", sheet("CN") firstrow


*To ensure that imported column names are distinct when datasets are merged
rename * *_CN

*Converting the index to a numeric variable, used as key to merge later
gen period = real(A_CN)
sort period

*Saving a copy of the resultsCN data into a dta file, to be merged later
save "mainfolder2/CNData.dta", replace

*Importing the results from bootstrapping EU simulation estimates of h-index
clear
import excel "mainfolder2/ResultsEU_10000.xlsx", sheet("ResultsEU") firstrow

*To ensure that imported column names are distinct when datasets are merged
rename * *_EU

*Converting the index to a numeric variable
gen period = real(A_EU)
save "mainfolder2/EUData.dta", replace

*Merging the datasets together, and saving it
merge 1:1 period using "mainfolder2/CNData.dta"
save "mainfolder2/MergedData.dta", replace

*Code to plot the superimposed graphs for CN bootstrap estimates and CI 
twoway line initial_CN period, color(red)   lwidth(thick) || rarea   bs_ciupper_CN    bs_cilower_CN period,  fcolor(red%50) lcolor(red%50)  lwidth(thick) ylabel(0(.2)1)  xtitle("Time Period", size(large)) ytitle("H-Index", size(large)) legend(off) scheme(s1color) aspectratio(0.5)
*twoway rarea   bs_ciupper_CN    bs_cilower_CN period,  fcolor(red%50) lcolor(red%50)|| line initial_CN period, color(red)   lwidth(thick) ylabel(0(.1)1)   xtitle("Time Periods") ytitle("h-index") legend(off) scheme(s1color)
*twoway rarea   bs_ciupper_CN    bs_cilower_CN period,  fcolor(red%50) lcolor(red%50)|| line initial_CN period, color(red)   lwidth(thick) ylabel(0(.1)1)  title("Plot of H-index over 500 periods - China (10000 iterations)") xtitle("Time Periods") ytitle("h-index") legend(order(2 "China" 1 "China 95% CI") pos(6) row(1)) scheme(s1color)

graph export "mainfolder2/hindex_graphs_china_10000.pdf", replace as(pdf)

*Code to plot the superimposed graphs for EU bootstrap estimates and CI 
twoway rarea   bs_ciupper_EU    bs_cilower_EU period,  fcolor(green%50) lcolor(green%50) ||  line initial_EU period, color(green)   lwidth(thick) ylabel(0(.1)1) xtitle("Time Period", size(large)) ytitle("H-Index", size(large)) legend(off) scheme(s1color) aspectratio(0.5)
*twoway rarea   bs_ciupper_EU    bs_cilower_EU period,  fcolor(green%50) lcolor(green%50) ||  line initial_EU period, color(green)   lwidth(thick) ylabel(0(.1)1) title("Plot of H-index over 500 periods - EU (10000 iterations)") xtitle("Time Periods") ytitle("h-index") legend(order(2 "EU" 1 "EU 95% CI") pos(6) row(1)) scheme(s1color)

graph export "mainfolder2/hindex_graphs_EU_10000.pdf", replace as(pdf)

*Code to plot the superimposed graphs for both CN and EU bootstrap estimates and CI 
twoway line initial_CN period, color(red)   lwidth(thick) || rarea   bs_ciupper_CN    bs_cilower_CN period,  fcolor(red%50) lcolor(red%50) || rarea   bs_ciupper_EU    bs_cilower_EU period,  fcolor(green%50) lcolor(green%50) || line initial_EU period, color(green)   lwidth(thick) ylabel(0(.2)1)  xtitle("Time Period", size(large)) ytitle("H-Index", size(large)) legend(off) scheme(s1color) aspectratio(0.5)
*twoway line initial_CN period, color(red)   lwidth(thick) || rarea   bs_ciupper_CN    bs_cilower_CN period,  fcolor(red%50) lcolor(red%50) || rarea   bs_ciupper_EU    bs_cilower_EU period,  fcolor(green%50) lcolor(green%50) || line initial_EU period, color(green)   lwidth(thick) ylabel(0(.1)1)  xtitle("Time Periods") ytitle("h-index") legend(off) scheme(s1color)
*twoway line initial_CN period, color(red)   lwidth(thick) || rarea   bs_ciupper_CN    bs_cilower_CN period,  fcolor(red%50) lcolor(red%50) || rarea   bs_ciupper_EU    bs_cilower_EU period,  fcolor(green%50) lcolor(green%50) || line initial_EU period, color(green)   lwidth(thick) ylabel(0(.1)1)  title("Plot of H-index over 500 periods (10000 iterations)") xtitle("Time Periods") ytitle("h-index") legend(order(1 "China" 2 "China 95% CI" 4 "EU" 3 "EU 95% CI") pos(6) row(2)) scheme(s1color)


graph export "mainfolder2/hindex_graphs_both_10000.pdf", replace as(pdf)
