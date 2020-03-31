
////////////////////////////////////////////////////////////////////////////////
// Change directory

cd "~/Documents/GitHub/event_study_graphs/"

///////////////////////////////////////////////////////////////////////////////
// Import program for event study
do "event_study_program.do"
do "unique_treated_units_by_time_program.do"


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

// Generate time treated 
gen time_treated  = ceil(24 * uniform()) + 1979
replace time_treated = . if ever_treated == 0

// Generate state component 
gen state_portion = rnormal()

// Expand dataset
expand 24

// Generate time 
bysort state: gen time = _n + 1979
sort state time

// Generate treated
gen treated = 0 
replace treated = 1 if time >= time_treated & ever_treated == 1

// Gen time portion 
gen time_portion = runiform()
bysort time: replace time_portion = time_portion[1]

// Gen control variable 
gen x = rnormal()

// Gen y 
gen y = 1.2 + 2*treated + .2*x + 3*state_portion + 2*time_portion + rnormal()


///////////////////////////////////////////////////////////////////////////////
// Create event study graph

create_event_study_graph, ///
	first_time(1990) ///
	last_time(2003) ///
	first_treated_time(1982) ///
	last_treated_time(2003) ///
	low_event_cap(-11) ///
	high_event_cap(11) ///
	low_event_cap_graph(-10) ///
	high_event_cap_graph(10) ///
	y_var(y) ///
	dd_control_list(none) ///
	control_list(x) ///
	fixed_effects("state time") ///
	cluster_var("state") ///
	x_label("Years since treatment") ///
	x_lab_step_size(2) ///
	y_label("Y-variable") ///
	graph_label("Event study estimates with imposed treatment effect of 2 at t = 0") ///
	output_file_name("event-study") ///
	time_variable(time) ///
	id_variable(state) ///
	weight_var(none)
	
	
create_unique_treat_by_yr_graph, ///
	low_event_cap(-11) ///
	high_event_cap(10) ///
	y_var(y) ///
	x_label("test") ///
	x_lab_step_size(2) ///
	y_label("test-y") ///
	graph_label("none") ///
	output_file_name("common-support") ///
	time_variable(time) ///
	id_variable(state) 

exit
///////////////////////////////////////////////////////////////////////////////
// Create event study graph after dropping event times

drop if event_time_bacon == -1 | event_time_bacon == 0 | event_time_bacon == 9
drop if event_time_bacon >= -8 &  event_time_bacon <= -3
keep state ever_treated time_treated state_portion time treated time_portion x y

create_event_study_graph, ///
	first_time(1990) ///
	last_time(2003) ///
	first_treated_time(1982) ///
	last_treated_time(2003) ///
	low_event_cap(-11) ///
	high_event_cap(11) ///
	low_event_cap_graph(-10) ///
	high_event_cap_graph(10) ///
	y_var(y) ///
	dd_control_list(none) ///
	control_list(x) ///
	fixed_effects("state time") ///
	cluster_var("state") ///
	x_label("Years since treatment") ///
	x_lab_step_size(2) ///
	y_label("Y-variable") ///
	graph_label("Event study estimates with imposed treatment effect of 2 at t = 0") ///
	output_file_name("event-study-missing-data") ///
	time_variable(time) ///
	id_variable(state) ///
	weight_var(none)
	

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
rename divyear time_treated	

create_event_study_graph, ///
	first_time(1960) ///
	last_time(1996) ///
	first_treated_time(1950) ///
	last_treated_time(2000) ///
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
	output_file_name("event-study-wolfers") ///
	time_variable(year) ///
	id_variable(st) ///
	weight_var(none)
	
