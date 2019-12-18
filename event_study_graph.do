// state is treated units
// year is time variable
// year_treated is date of treatment 
// treated is ever_treated x post_treatment


capture program drop create_event_study_graph
program create_event_study_graph, 
syntax, ///
	first_year(real) ///
	last_year(real) ///
	first_treated_year(real) ///
	last_treated_year(real) ///
	low_event_cap(real) ///
	high_event_cap(real) ///
	low_event_cap_graph(real) ///
	high_event_cap_graph(real) ///
	y_var(string) ///
	dd_control_list(string) ///
	control_list(string) ///
	fixed_effects(string) ///
	cluster_var(string) ///
	x_label(string) ///
	x_lab_step_size(real) ///
	y_label(string) ///
	graph_label(string) ///
	output_file_name(string) 

	////////////////////////////////////////////////////////////////////////////
	// Set up timing variables
	
	capture confirm variable year_treated
	if !_rc {
		di "year_treated exists"
	}
	else {
		*Make an adoption year variable.
		capture drop exp_year
		gen exp_year = year if ever_treated == 1
		egen year_treated = min(exp_year), by(state)
	}
	
	
	capture confirm variable post
	if !_rc {
		di in red "post"
	}
	else {
		di in red "post does not exist"
	}
	
	if "`dd_control_list'" == "none" {
		local dd_control_list ""
	}
	
	if "`control_list'" == "none" {
		local control_list ""
	}
	
	
	* gen event-time (adoption_date = treatment adoption date )
	gen event_time_bacon = year - year_treated

	* make sure untreated units are included, but also don't get dummies (by giving them "-1")
	recode event_time_bacon (.=-1) (-1000/`low_event_cap'=`low_event_cap') (`high_event_cap'/1000=`high_event_cap')

	*ensure that "xi" omits -1
	char event_time_bacon[omit] -1

	*mak dummies
	xi i.event_time_bacon, pref(_T)

	* Position of -2
	if `low_event_cap' > ( `first_year' - `last_treated_year') {
	 local pos_of_neg_2 = abs(2 - (-`low_event_cap' + 1))
	}	

	if `low_event_cap' <= ( `first_year' - `last_treated_year') {
	 local pos_of_neg_2 = abs(2 - (-( `first_year' - `last_treated_year') + 1))
	}

	* Position of 0
	local pos_of_zero = `pos_of_neg_2' + 2

	* Position of max
	if `high_event_cap' >= ( `last_year' - `first_treated_year') {
	 local pos_of_max = `last_year' - `first_treated_year' + `pos_of_zero'
	}	

	if `high_event_cap' < ( `last_year' - `first_treated_year') {
	 local pos_of_max = `high_event_cap' + `pos_of_zero'
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Run DD regression
	reghdfe `y_var' ///
		treated ///
		`dd_control_list' ///
		`control_list' ///
	, absorb(`fixed_effects') cluster(`cluster_var')
	
	// Store DD Coefs and 95% CI and p-value

	
	/////////////////////////////////////////////////////////////////////////////
	// Run event-study regression
	reghdfe `y_var' ///
		_T* ///
		`control_list' ///
		`dd_control_list' ///
	, absorb(`fixed_effects') cluster(`cluster_var')

	// Store estimates
	capture drop order
	capture drop b 
	capture drop high 
	capture drop low

	gen order = .
	gen b =. 
	gen high =. 
	gen low =.
	
	* Graph start position
	local graph_start  = `low_event_cap_graph' - `low_event_cap' + 1
	
	* Position of -2
	local pos_of_neg_2 = abs(2 - (-`low_event_cap' + 1))

	* Position of 0
	local pos_of_zero = `pos_of_neg_2' + 2

	* Position of max
	 local pos_of_max = `high_event_cap_graph' + `pos_of_zero'
	
	
	local i = 1
	
	forvalues j = `graph_start'(1)`pos_of_neg_2'{
		local event_time = `j' - 2 -  `pos_of_neg_2'
		
		replace order = `event_time' in `i'
		
		replace b    = _b[_Tevent_tim_`j'] in `i'
		replace high = _b[_Tevent_tim_`j'] + 1.96*_se[_Tevent_tim_`j'] in `i'
		replace low  = _b[_Tevent_tim_`j'] - 1.96*_se[_Tevent_tim_`j'] in `i'
			
		local i = `i' + 1
	}

	replace order = -1 in `i'

	replace b    = 0  in `i'
	replace high = 0  in `i'
	replace low  = 0  in `i'

	local i = `i' + 1

	forvalues j = `pos_of_zero'(1)`pos_of_max'{
		local event_time = `j' - 2 -  `pos_of_neg_2'

		replace order = `event_time' in `i'
		
		replace b    = _b[_Tevent_tim_`j'] in `i'
		replace high = _b[_Tevent_tim_`j'] + 1.96*_se[_Tevent_tim_`j'] in `i'
		replace low  = _b[_Tevent_tim_`j'] - 1.96*_se[_Tevent_tim_`j'] in `i'
			
			
		local i = `i' + 1
	}
	
	// Plot estimates	
	twoway rarea low high order if order <= `high_event_cap_graph' & order >= `low_event_cap_graph' , fcol(gs14) lcol(white) msize(3) /// estimates
		|| connected b order if order <= `high_event_cap_graph' & order >= `low_event_cap_graph', lw(1.1) col(white) msize(7) msymbol(s) lp(solid) /// highlighting
		|| connected b order if order <= `high_event_cap_graph' & order >= `low_event_cap_graph', lw(0.6) col("71 71 179") msize(5) msymbol(s) lp(solid) /// connect estimates
		|| scatteri 0 `low_event_cap_graph' 0 `high_event_cap_graph', recast(line) lcol(gs8) lp(longdash) lwidth(0.5) /// zero line 
			xlab(`low_event_cap_graph'(`x_lab_step_size')`high_event_cap_graph' ///
					, nogrid labsize(7) angle(0)) ///
			ylab(, nogrid labs(7)) ///
			legend(off) ///
			xtitle("`x_label'", size(5)) ///
			ytitle("`y_label'", size(5)) ///
			subtitle("`graph_label'", size(5) pos(11)) ///
			xline(-.5, lpattern(dash) lcolor(gs7) lwidth(1)) 	
			
	// Save graph
	graph export "`output_file_name'.png", replace

end
////////////////////////////////////////////////////////////////////////////////
// Change directory

cd "~/Documents/GitHub/event_study_graphs/"


////////////////////////////////////////////////////////////////////////////////
// Run with test data 

// Set seed
set seed 1234

// Create fake data
clear

// Set 51 observations
set obs 51

// Generate id 
gen state = _n

// Generate treated
gen random_draw = runiform()
gen ever_treated = 0
replace ever_treated = 1 if random_draw > .75
drop random_draw

// Generate year treated 
gen year_treated  = ceil(24 * uniform()) + 1979
replace year_treated = . if ever_treated == 0

// Generate state component 
gen state_portion = rnormal()

// Expand dataset
expand 24

// Generate year 
bysort state: gen year = _n + 1979
sort state year

// Generate treated
gen treated = 0 
replace treated = 1 if year >= year_treated & ever_treated == 1

// Gen year portion 
gen year_portion = runiform()
bysort year: replace year_portion = year_portion[1]

// Gen control variable 
gen x = rnormal()

// Gen y 
gen y = 1.2 + 2*treated + .2*x + 3*state_portion + 2*year_portion + rnormal()


create_event_study_graph, ///
	first_year(1980) ///
	last_year(2003) ///
	first_treated_year(1982) ///
	last_treated_year(2003) ///
	low_event_cap(-11) ///
	high_event_cap(11) ///
	low_event_cap_graph(-10) ///
	high_event_cap_graph(10) ///
	y_var(y) ///
	dd_control_list(none) ///
	control_list(x) ///
	fixed_effects("state year") ///
	cluster_var("state") ///
	x_label("Years since treatment") ///
	x_lab_step_size(2) ///
	y_label("Y-variable") ///
	graph_label("Event study estimates with imposed treatment effect of 2 at t = 0") ///
	output_file_name("event-study") 

////////////////////////////////////////////////////////////////////////////////
// Run test data 

// Open Wolfers data (http://users.nber.org/~jwolfers/data/DivorceData(QJE).zip)
clear
use "Til Death MASTER DATA streamlined.dta" 

// Keep only women
keep if sex == 2

// Generate treated
gen ever_treated = 0 
replace ever_treated = 1 if !missing(divyear)

// Generate post 
gen post = 0 
replace post = 1 if year >= divyear & !missing(divyear)

// Generate treated 
gen treated = 0 
replace treated = 1 if post == 1 & ever_treated == 1

// Rename divyear to year_treated
rename divyear year_treated


create_event_study_graph, ///
	first_year(1960) ///
	last_year(1996) ///
	first_treated_year(1950) ///
	last_treated_year(2000) ///
	low_event_cap(-11) ///
	high_event_cap(21) ///
	low_event_cap_graph(-10) ///
	high_event_cap_graph(20) ///
	y_var(suiciderate_elast_jag) ///
	dd_control_list(none) ///
	control_list("afdcmax afdcrolls roevwade femep stateur lnpersinc sa* sh*") ///
	fixed_effects("st year") ///
	cluster_var("st") ///
	x_label("Years since unilateral divorce law") ///
	x_lab_step_size(5) ///
	y_label("Aggregate suicide rate") ///
	graph_label("Change in suicide rate by years since law change") ///
	output_file_name("event-study-wolfers")
	
