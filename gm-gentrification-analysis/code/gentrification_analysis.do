/*==============================================================================
Gentrification and Housing Affordability Analysis
Author: Sumedh Brahmadevara
Description: Analysis of gentrification impacts on housing affordability using 
             MSOA-level census data (1991-2021)
==============================================================================*/

clear all
set more off

/*------------------------------------------------------------------------------
Step 1: Data Preparation - Link MSOA codes to shapefile IDs
------------------------------------------------------------------------------*/

* Set working directory (users should modify this path)
global data_path "[INSERT_YOUR_DATA_PATH]"
global shape_path "[INSERT_YOUR_SHAPE_PATH]"

* Load and prepare MSOA identification data
use "$shape_path/gm_msoa_data.dta", clear
keep msoa11cd _ID
duplicates drop
save "$shape_path/msoa_ids.dta", replace

/*------------------------------------------------------------------------------
Step 2: Merge census panel data with spatial identifiers
------------------------------------------------------------------------------*/

use "$data_path/census_panel.dta", clear
merge m:1 msoa11cd using "$shape_path/msoa_ids.dta"
keep if _merge == 3
drop _merge
save "$data_path/census_panel_spatial.dta", replace

/*------------------------------------------------------------------------------
Step 3: Panel data setup and cleaning
------------------------------------------------------------------------------*/

use "$data_path/census_panel_spatial.dta", clear
drop occ_other

* Set up panel structure
xtset _ID year
xtdescribe 

* Check for balanced panel
egen n_years = count(year), by(_ID)
list _ID if n_years != 4

* Remove incomplete observations
drop if _ID == 6109
xtset _ID year
xtdescribe

* Create alternative MSOA identifier
encode msoa11cd, gen(msoa_id)
xtset msoa_id year 

* Set up spatial structure
spset _ID
spset, modify shpfile(msoa_coords)

/*------------------------------------------------------------------------------
Step 4: Variable Construction
------------------------------------------------------------------------------*/

* Housing affordability measure
gen affordability_ratio = (weekly_income * 52) / houseprice
gen log_houseprice = log(houseprice)

* Gentrification components
gen pct_degree = (educ_lev3 + educ_lev4) / (educ_lev1 + educ_lev2 + educ_lev3 + educ_lev4)
gen pct_prof_occ = (occ_professional + occ_managers_seniors + occ_skilled) / occ_total
gen pct_25_44 = age25_44 / pop
gen pct_childless = (hh_married_nochildren + hh_1person_other) / hh_total

/*------------------------------------------------------------------------------
Step 5: Standardized Gentrification Index (Base Year: 1991)
------------------------------------------------------------------------------*/

* Initialize storage variables
foreach var in mean_deg sd_deg mean_occ sd_occ mean_ind sd_ind mean_age sd_age mean_child sd_child {
    gen `var' = .
}

* Calculate 1991 baseline statistics
quietly {
    summarize pct_degree if year == 1991
    replace mean_deg = r(mean)
    replace sd_deg = r(sd)

    summarize pct_prof_occ if year == 1991
    replace mean_occ = r(mean)
    replace sd_occ = r(sd)

    summarize quaternary_share if year == 1991
    replace mean_ind = r(mean)
    replace sd_ind = r(sd)

    summarize pct_25_44 if year == 1991
    replace mean_age = r(mean)
    replace sd_age = r(sd)

    summarize pct_childless if year == 1991
    replace mean_child = r(mean)
    replace sd_child = r(sd)
}

* Forward-fill baseline statistics
foreach var in mean_deg sd_deg mean_occ sd_occ mean_ind sd_ind mean_age sd_age mean_child sd_child {
    replace `var' = `var'[_n-1] if missing(`var')
}

* Create standardized scores
gen z_edu    = (pct_degree - mean_deg) / sd_deg
gen z_occ    = (pct_prof_occ - mean_occ) / sd_occ
gen z_ind    = (quaternary_share - mean_ind) / sd_ind
gen z_age    = (pct_25_44 - mean_age) / sd_age
gen z_child  = (pct_childless - mean_child) / sd_child

* Main gentrification index (3-component)
gen gentrification_index = (z_edu + z_age + z_child) / 3

/*------------------------------------------------------------------------------
Step 6: Change-based Gentrification Measures
------------------------------------------------------------------------------*/

* Store baseline values for each MSOA
foreach var in pct_degree pct_prof_occ quaternary_share pct_25_44 pct_childless {
    gen base_`var' = .
    replace base_`var' = `var' if year == 1991
    bysort _ID (year): replace base_`var' = base_`var'[_n-1] if missing(base_`var')
}

* Calculate changes from baseline
gen d_degree = pct_degree - base_pct_degree
gen d_occ    = pct_prof_occ - base_pct_prof_occ
gen d_ind    = quaternary_share - base_quaternary_share
gen d_age    = pct_25_44 - base_pct_25_44
gen d_child  = pct_childless - base_pct_childless

* Change-based gentrification index
gen gentrification_change_index = (d_degree + d_occ + d_ind + d_age + d_child) / 5

/*------------------------------------------------------------------------------
Step 7: Housing Market Variables
------------------------------------------------------------------------------*/

* Private rental share (L variable)
gen L = (private_rented_total / total_households) * 100

* Lagged variables
sort msoa11cd year
gen gentrification_index_lag = .
gen L_lag = .
bysort msoa11cd (year): replace gentrification_index_lag = gentrification_index[_n-1]
bysort msoa11cd (year): replace L_lag = L[_n-1]

* Private rental changes
bysort msoa11cd (year): gen private_rented_change = private_rented_total - private_rented_total[_n-1]
gen private_rented_change_pct = private_rented_change / total_households

* Housing supply measures (S variables)
gen net_add_dwellings = dwelling_total - dwelling_total[_n-1]
gen net_add_dwellings_pc = (net_add_dwellings / dwelling_total[_n-1]) * 100
gen vacancy_rate = (dwelling_total - total_emp) / dwelling_total
gen structural_supply = (type_detach + type_semi) / dwelling_total

/*------------------------------------------------------------------------------
Step 8: Interaction Terms
------------------------------------------------------------------------------*/

gen S_G = structural_supply * gentrification_index
gen L_G = L * gentrification_index

/*------------------------------------------------------------------------------
Step 9: Control Variables
------------------------------------------------------------------------------*/

gen pct_young = (age0_14 + age15_24 + age25_44 + age45_64) / pop
gen pct_families = (hh_married_withchildren + hh_married_nochildren) / hh_total
gen unemp_rate = 1 - (occ_total / pop)
gen log_pop_density = log(pop / area_km2)

* Standardized amenity variables
egen z_school_count = std(school_count_1km)
egen z_pct_greenspace = std(pct_greenspace)

/*------------------------------------------------------------------------------
Step 10: Descriptive Statistics
------------------------------------------------------------------------------*/

* Summary statistics
summarize affordability_ratio gentrification_index_lag private_rented_change_pct ///
    net_add_dwellings_pc pct_young pct_families log_pop_density ///
    z_school_count z_pct_greenspace dist_to_station

* Correlation matrix
pwcorr affordability_ratio gentrification_index_lag private_rented_change_pct ///
    net_add_dwellings_pc pct_young pct_families log_pop_density ///
    z_school_count z_pct_greenspace dist_to_station, sig

* Distribution of gentrification index
summarize gentrification_index, detail
histogram gentrification_index, width(0.5) frequency normal ///
    title("Distribution of Gentrification Index")

gen gentr_index_rounded = round(gentrification_index, 0.5)
tabulate gentr_index_rounded

/*------------------------------------------------------------------------------
Step 11: Spatial Weight Matrices
------------------------------------------------------------------------------*/

* Create spatial weight matrices
spmatrix create contiguity W if year == 1991, normalize(row)
spmatrix create idistance M if year == 1991, normalize(row)

/*------------------------------------------------------------------------------
Step 12: Regression Analysis
------------------------------------------------------------------------------*/

* Standard fixed effects model
xtreg affordability_ratio gentrification_index L net_add_dwellings_pc ///
    pct_young pct_families log_pop_density z_school_count z_pct_greenspace, fe

* Spatial autoregressive model with contiguity weights
spxtregress affordability_ratio gentrification_index L net_add_dwellings_pc ///
    pct_young pct_families log_pop_density z_school_count z_pct_greenspace ///
    dist_to_station, fe dvarlag(W) force
estimates store model_sar
etable, estimates(model_sar) showstars showstarsnote ///
    stars(.1 "*" .05 "**" .01 "***", attach(_r_b)) mstat(r2)

* Spatial Durbin model with lagged variables
spxtregress affordability_ratio gentrification_index_lag L_lag net_add_dwellings_pc ///
    pct_young pct_families log_pop_density z_school_count z_pct_greenspace ///
    dist_to_station, fe dvarlag(W) ///
    ivarlag(W:gentrification_index_lag L_lag net_add_dwellings_pc) force

/*------------------------------------------------------------------------------
Step 13: Alternative Dependent Variable (Log Specification)
------------------------------------------------------------------------------*/

* Create time index for dynamic models
gen time_index = .
replace time_index = 1 if year == 1991
replace time_index = 2 if year == 2001
replace time_index = 3 if year == 2011
replace time_index = 4 if year == 2021
xtset _ID time_index

* Log affordability measure
gen log_income = log(weekly_income)
gen afford = log_income - log_houseprice
gen L_afford = L.afford

/*------------------------------------------------------------------------------
Step 14: Model Comparison and Diagnostics
------------------------------------------------------------------------------*/

* Basic model progression
quietly xtreg afford gentrification_index dist_to_station unemp_rate i.year, fe
estimates store m1
estat ic

quietly xtreg afford gentrification_index structural_supply L dist_to_station ///
    unemp_rate i.year, fe
estimates store m2
estat ic

quietly xtreg afford gentrification_index structural_supply S_G L L_G ///
    dist_to_station unemp_rate i.year, fe
estimates store m3
estat ic

* Hausman test for FE vs RE
quietly xtreg afford gentrification_index structural_supply S_G L L_G ///
    dist_to_station unemp_rate i.year, re
estimates store re_model
quietly xtreg afford gentrification_index structural_supply S_G L L_G ///
    dist_to_station unemp_rate i.year, fe
estimates store fe_model
hausman fe_model re_model

* Dynamic panel model (LSDVC)
xtlsdvc afford L_afford gentrification_index structural_supply S_G L L_G ///
    dist_to_station unemp_rate, initial(ab) vcov(50)

/*------------------------------------------------------------------------------
Step 15: Diagnostic Tests
------------------------------------------------------------------------------*/

* Serial correlation test
xtserial afford gentrification_index dist_to_station unemp_rate

* Heteroskedasticity test
xttest3

* Model specification tests
linktest
estat ovtest

/*------------------------------------------------------------------------------
Step 16: Spatial Visualization
------------------------------------------------------------------------------*/

* Create maps for each year
levelsof year, local(years)
local maps_gent

foreach y of local years {
    preserve
        keep if year == `y'
        
        if `y' == 2011 {
            spmap gentrification_index using "$shape_path/msoa_coords.dta", ///
                id(_ID) fcolor(YlOrRd) clmethod(custom) ///
                clbreaks(-0.15 -0.04 0 0.04 0.08 0.18 0.30 0.40) ///
                osize(vthin) title("`y'", size(small)) ///
                legend(title("Gentrification Index", size(vsmall)) ///
                       size(vsmall) ring(0) pos(5)) ///
                name(map_gent_`y', replace)
        }
        else {
            spmap gentrification_index using "$shape_path/msoa_coords.dta", ///
                id(_ID) fcolor(YlOrRd) clmethod(custom) ///
                clbreaks(-0.15 -0.04 0 0.04 0.08 0.18 0.30 0.40) ///
                osize(vthin) title("`y'", size(small)) legend(off) ///
                name(map_gent_`y', replace)
        }
    restore
    
    local maps_gent `maps_gent' map_gent_`y'
}

* Combine maps into panel
graph combine `maps_gent', cols(2) ///
    title("Gentrification Index by MSOA (1991–2021)", size(medium)) ///
    iscale(1) imargin(0 0 0 0) xsize(12) ysize(10)

graph export "gentrification_index_panel.png", width(2500) replace

/*==============================================================================
End of Analysis
==============================================================================*/
