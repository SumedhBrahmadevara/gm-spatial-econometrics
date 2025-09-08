/*==============================================================================
MSOA Data Processing Pipeline
================================================================================
Purpose: Process and merge MSOA (Middle Layer Super Output Area) census data
         across multiple time periods (1991, 2001, 2011, 2021)
Author:  Sumedh Brahmadevara
Version: 2.0

Notes:   This script consolidates data processing for all census years and
         house price data into a unified pipeline
==============================================================================*/

clear all
set more off
version 16.0

/*------------------------------------------------------------------------------
Global Path Configuration
------------------------------------------------------------------------------*/
// Set relative paths - adjust these to match your project structure
global PROJECT_ROOT "."
global RAW_DATA     "${PROJECT_ROOT}/data/raw"
global EXT_DATA     "${PROJECT_ROOT}/data/external"
global PROCESSED    "${PROJECT_ROOT}/data/processed"
global MAPPING      "${PROJECT_ROOT}/data/mapping"

// Create output directories if they don't exist
foreach dir in "${PROCESSED}" "${PROCESSED}/1991" "${PROCESSED}/2001" "${PROCESSED}/2011" "${PROCESSED}/2021" {
    cap mkdir "`dir'"
}

/*------------------------------------------------------------------------------
Program: convert_to_dta
Purpose: Convert various file formats to Stata format with error handling
------------------------------------------------------------------------------*/
program define convert_to_dta
    args filepath filetype varname_row
    
    local filename = subinstr("`filepath'", ".xlsx", "", .)
    local filename = subinstr("`filename'", ".csv", "", .)
    
    cap confirm file "`filepath'"
    if _rc {
        di as error "File not found: `filepath'"
        exit 601
    }
    
    if "`filetype'" == "excel" {
        cap import excel using "`filepath'", firstrow clear
        if _rc {
            di as error "Failed to import Excel file: `filepath'"
            exit _rc
        }
        cap destring _all, replace
    }
    else if "`filetype'" == "csv" {
        if "`varname_row'" == "" local varname_row = 1
        cap import delimited using "`filepath'", varnames(`varname_row') clear
        if _rc {
            di as error "Failed to import CSV file: `filepath'"
            exit _rc
        }
        cap drop v1  // Remove potential row number column
        cap destring _all, replace
    }
    
    save "`filename'.dta", replace
    di as result "Converted: `filepath' -> `filename'.dta"
end

/*------------------------------------------------------------------------------
Program: merge_census_files
Purpose: Merge all census files for a given year with consistent ID variable
------------------------------------------------------------------------------*/
program define merge_census_files
    args year id_var file_list raw_path
    
    di as text _n "Processing `year' census data..."
    
    // Start with first file (age)
    use "`raw_path'/age.dta", clear
    
    // Merge remaining files
    foreach file in `file_list' {
        cap confirm file "`raw_path'/`file'.dta"
        if !_rc {
            merge 1:1 `id_var' using "`raw_path'/`file'.dta", nogen
            di as result "  Merged: `file'.dta"
        }
        else {
            di as error "  Warning: `file'.dta not found, skipping"
        }
    }
    
    gen year = `year'
    di as result "Added year variable: `year'"
end

/*------------------------------------------------------------------------------
Main Processing Pipeline
------------------------------------------------------------------------------*/

/*============================
2021 Census Data Processing
============================*/
di as text "{hline 80}"
di as text "PROCESSING 2021 CENSUS DATA"
di as text "{hline 80}"

// Convert source files
local file_types_2021 "age comp educ occ tenure htype"
foreach f of local file_types_2021 {
    convert_to_dta "${RAW_DATA}/2021/`f'.xlsx" "excel"
}
convert_to_dta "${RAW_DATA}/2021/ind.csv" "csv"

// Merge census files
merge_census_files 2021 "msoa21cd" "comp educ occ tenure htype ind" "${RAW_DATA}/2021"

// Map MSOA21 to MSOA11 and aggregate
cap merge 1:1 msoa21cd using "${MAPPING}/msoa21_to_msoa11_mapping.dta", nogen
replace msoa11cd = msoa21cd if missing(msoa11cd)

// Aggregate to MSOA11 level
ds, has(type numeric)
local numeric_vars `r(varlist)'
collapse (sum) `numeric_vars', by(msoa11cd)
gen year = 2021

// Add external data
local external_files_2021 "houseprice_2021 spatial_2021 income_2021"
foreach x of local external_files_2021 {
    cap merge 1:1 msoa11cd using "${EXT_DATA}/`x'.dta", nogen
    if _rc == 0 {
        di as result "  Added external data: `x'"
    }
}

save "${PROCESSED}/2021/msoa11_2021.dta", replace
di as result "Saved: ${PROCESSED}/2021/msoa11_2021.dta"

/*============================
2011 Census Data Processing  
============================*/
di as text "{hline 80}"
di as text "PROCESSING 2011 CENSUS DATA"
di as text "{hline 80}"

// Convert source files (all CSV format for 2011)
local file_types_2011 "age comp educ occ tenure htype ind"
foreach f of local file_types_2011 {
    convert_to_dta "${RAW_DATA}/2011/`f'.csv" "csv"
}

// Merge census files
merge_census_files 2011 "msoa11cd" "comp educ occ tenure htype ind" "${RAW_DATA}/2011"

// Add external data
local external_files_2011 "houseprice_2011 spatial_2011 income_2011"
foreach x of local external_files_2011 {
    cap merge 1:1 msoa11cd using "${EXT_DATA}/`x'.dta", nogen
    if _rc == 0 {
        di as result "  Added external data: `x'"
    }
}

save "${PROCESSED}/2011/msoa11_2011.dta", replace
di as result "Saved: ${PROCESSED}/2011/msoa11_2011.dta"

/*============================
2001 Census Data Processing
============================*/
di as text "{hline 80}"
di as text "PROCESSING 2001 CENSUS DATA"  
di as text "{hline 80}"

// Convert source files
local file_types_2001 "age comp educ occ tenure htype"
foreach f of local file_types_2001 {
    convert_to_dta "${RAW_DATA}/2001/`f'.xlsx" "excel"
}
convert_to_dta "${RAW_DATA}/2001/ind.csv" "csv"

// Merge census files
merge_census_files 2001 "msoa01cd" "comp educ occ tenure htype ind" "${RAW_DATA}/2001"

// Handle MSOA boundary changes (2001 -> 2011)
di as text "Handling MSOA boundary splits..."

// Define areas that need to be split
local split2_codes "E02001054 E02001058"
local split3_code "E02001060"

// Expand observations for split areas
expand 2 if inlist(msoa01cd, `split2_codes')
expand 3 if msoa01cd == "`split3_code'"

// Adjust values proportionally for splits
ds, has(type numeric)
local numeric_vars `r(varlist)'
foreach v of local numeric_vars {
    replace `v' = `v'/2 if inlist(msoa01cd, `split2_codes')
    replace `v' = `v'/3 if msoa01cd == "`split3_code'"
}

// Create unique identifiers for duplicated areas
duplicates tag msoa01cd, gen(dup_tag)
gen obs_id = _n
bysort msoa01cd (obs_id): gen dup_index = _n if dup_tag > 0
gen msoa11cdx = msoa01cd
replace msoa11cdx = msoa01cd + "_dup" + string(dup_index) if !missing(dup_index)

// Map specific duplicate codes to new MSOA11 codes
local mapping_pairs `" "E02001054_dup1" "E02006913" "E02001054_dup2" "E02006915" "E02001058_dup1" "E02006902" "E02001058_dup2" "E02006912" "E02001060_dup1" "E02006914" "E02001060_dup2" "E02006916" "E02001060_dup3" "E02006917" "'

local n_pairs = wordcount(`"`mapping_pairs'"')/2
forval i = 1/`n_pairs' {
    local old_code = word(`"`mapping_pairs'"', `i'*2-1)
    local new_code = word(`"`mapping_pairs'"', `i'*2)
    replace msoa11cdx = "`new_code'" if msoa11cdx == "`old_code'"
}

// Apply MSOA01 to MSOA11 mapping
merge m:1 msoa01cd using "${MAPPING}/msoa01merge_to_msoa11_mapping.dta", nogen
replace msoa11cd = msoa11cdx if !missing(msoa11cdx)
replace msoa11cd = msoa01cd if missing(msoa11cd)

// Collapse to MSOA11 level
collapse (sum) `numeric_vars', by(msoa11cd)
drop dup_tag obs_id dup_index
gen year = 2001

// Add external data
local external_files_2001 "houseprice_2001 spatial_2001 income_2001"
foreach x of local external_files_2001 {
    cap merge 1:1 msoa11cd using "${EXT_DATA}/`x'.dta", nogen
    if _rc == 0 {
        di as result "  Added external data: `x'"
    }
}

save "${PROCESSED}/2001/msoa11_2001.dta", replace
di as result "Saved: ${PROCESSED}/2001/msoa11_2001.dta"

/*============================
1991 Census Data Processing
============================*/
di as text "{hline 80}"
di as text "PROCESSING 1991 CENSUS DATA"
di as text "{hline 80}"

// Import and clean 1991 data (single file format)
import delimited using "${RAW_DATA}/1991/1991.csv", varnames(1) clear
rename geo_code msoa11cd

// Remove "_w" suffix from variable names (weighted variables)
ds, has(type numeric)
local numeric_vars `r(varlist)'
foreach v of local numeric_vars {
    if substr("`v'", -2, .) == "_w" {
        local new_name = substr("`v'", 1, length("`v'") - 2)
        rename `v' `new_name'
        di as result "  Renamed: `v' -> `new_name'"
    }
}

gen year = 1991

// Add external data
local external_files_1991 "houseprice_1991 spatial_1991 income_1991"
foreach x of local external_files_1991 {
    cap merge 1:1 msoa11cd using "${EXT_DATA}/`x'.dta", nogen
    if _rc == 0 {
        di as result "  Added external data: `x'"
    }
}

save "${PROCESSED}/1991/msoa11_1991.dta", replace
di as result "Saved: ${PROCESSED}/1991/msoa11_1991.dta"

/*==============================================================================
House Price Data Processing Module
==============================================================================*/
di as text "{hline 80}"
di as text "PROCESSING HOUSE PRICE DATA"
di as text "{hline 80}"

// Create house price output directory
cap mkdir "${PROCESSED}/houseprices"

// Process each sheet from the house price Excel file
local sheets "1a 2a 3a"
local sheet_names "Sheet1a Sheet2a Sheet3a"

foreach s of local sheets {
    di as text "Processing house price sheet: `s'"
    
    // Import selected sheet with specific range
    cap import excel using "${RAW_DATA}/houseprices/medianhousepricesbymiddlelayersuperoutputarea.xlsx", ///
        sheet("`s'") cellrange(A4341:DP4693) firstrow clear
    
    if _rc {
        di as error "Could not import sheet `s', skipping..."
        continue
    }
    
    // Keep only required columns (MSOA codes, names, and price columns)
    keep C D E I M Q U Y AC AG AK AO AS AW BA BE BI BM BQ BU BY CC CG CK CO CS CW DA DE DI DM
    
    // Rename ID and name columns
    rename C msoa21cd
    rename D msoa21nm
    
    // Rename year columns systematically
    local price_vars "E I M Q U Y AC AG AK AO AS AW BA BE BI BM BQ BU BY CC CG CK CO CS CW DA DE DI DM"
    local start_year = 1995
    local col_num = 0
    
    foreach var of local price_vars {
        cap confirm variable `var'
        if !_rc {
            local target_year = `start_year' + `col_num'
            rename `var' houseprice`target_year'
            local ++col_num
        }
    }
    
    // Clean price data (remove currency symbols and commas)
    ds houseprice*, has(type string)
    foreach v of varlist `r(varlist)' {
        replace `v' = subinstr(`v', ",", "", .)
        replace `v' = subinstr(`v', "£", "", .)
        replace `v' = trim(`v')
        replace `v' = "" if inlist(`v', "..", "N/A", "#N/A")
        destring `v', replace force
    }
    
    // Reshape to long format for trend extrapolation
    reshape long houseprice, i(msoa21cd) j(year)
    gen fit_sample = inrange(year, 1995, 2023)
    
    // Fit linear trends and extrapolate to 1991-1994
    tempfile trend_coefs
    statsby _b_cons=_b[_cons] _b_slope=_b[year], by(msoa21cd) saving("`trend_coefs'"): ///
        regress houseprice year if fit_sample
    
    // Generate missing early years using trend extrapolation
    preserve
        keep msoa21cd
        duplicates drop
        expand 4
        bysort msoa21cd: gen year = 1995 - _n
        merge 1:1 msoa21cd using "`trend_coefs'", nogen keep(match)
        gen houseprice = _b_cons + _b_slope * year
        keep msoa21cd year houseprice
        keep if inrange(year, 1991, 1994)
        tempfile early_years
        save "`early_years'", replace
    restore
    
    // Append extrapolated early years
    append using "`early_years'"
    
    // Reshape back to wide format
    reshape wide houseprice, i(msoa21cd) j(year)
    
    // Map MSOA21 to MSOA11
    merge 1:1 msoa21cd using "${MAPPING}/msoa21_to_msoa11_mapping.dta", nogen
    replace msoa11cd = msoa21cd if missing(msoa11cd)
    
    // Aggregate to MSOA11 level
    ds houseprice*, has(type numeric)
    local price_vars `r(varlist)'
    collapse (mean) `price_vars', by(msoa11cd)  // Use mean for price aggregation
    
    // Keep only target census years and export
    keep msoa11cd houseprice1991 houseprice2001 houseprice2011 houseprice2021
    
    // Save processed data
    save "${PROCESSED}/houseprices/houseprice_`s'.dta", replace
    export delimited using "${PROCESSED}/houseprices/houseprice_`s'.csv", replace
    
    di as result "Processed and saved house price data for sheet `s'"
}

/*==============================================================================
Final Summary and Validation
==============================================================================*/
di as text "{hline 80}"
di as text "PROCESSING COMPLETE - SUMMARY"
di as text "{hline 80}"

// Check output files exist
local expected_files "${PROCESSED}/1991/msoa11_1991.dta ${PROCESSED}/2001/msoa11_2001.dta ${PROCESSED}/2011/msoa11_2011.dta ${PROCESSED}/2021/msoa11_2021.dta"

foreach file of local expected_files {
    cap confirm file "`file'"
    if !_rc {
        di as result "✓ Created: `file'"
        use "`file'", clear
        di as text "  - Observations: " as result _N
        di as text "  - Variables: " as result c(k)
    }
    else {
        di as error "✗ Missing: `file'"
    }
}

di as text _n "House price files:"
foreach s in 1a 2a 3a {
    local file "${PROCESSED}/houseprices/houseprice_`s'.dta"
    cap confirm file "`file'"
    if !_rc {
        di as result "✓ Created: `file'"
    }
    else {
        di as error "✗ Missing: `file'"
    }
}

di as text _n "{hline 80}"
di as text "All processing tasks completed!"
di as text "Check the ${PROCESSED} directory for output files."
di as text "{hline 80}"

/*==============================================================================
End of Script
==============================================================================*/
