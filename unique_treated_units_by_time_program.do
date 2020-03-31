// state is treated units
// time is time variable
// time_treated is date of treatment 
// treated is ever_treated x post_treatment


capture program drop create_unique_treat_by_yr_graph
program create_unique_treat_by_yr_graph, 
syntax, ///
	low_event_cap(real) ///
	high_event_cap(real) ///
	y_var(string) ///
	x_label(string) ///
	x_lab_step_size(real) ///
	y_label(string) ///
	graph_label(string) ///
	output_file_name(string) ///
	time_variable(string) ///
	id_variable(string) 

	////////////////////////////////////////////////////////////////////////////
// Set up variables
	
	capture confirm variable time_treated
	if !_rc {
		di "time_treated exists"
	}
	else {
		*Make an adoption time variable.
		capture drop temp_time
		gen temp_time = `time_variable' if ever_treated == 1
		egen time_treated = min(temp_time), by(`id_variable')
		drop temp_time
	}
	
	
	* gen event-time (adoption_date = treatment adoption date )
	capture drop event_time_bacon
	gen event_time_bacon = `time_variable' - time_treated

	* make sure untreated units are included, but also don't get dummies (by giving them "-1")
	recode event_time_bacon (-1000/`low_event_cap'=`low_event_cap') (`high_event_cap'/1000=`high_event_cap')

	// Determine number of unique "event times"
	levelsof(event_time_bacon) if !missing(`y_var') , matrow(event_time_names)
	mat list event_time_names
				
	local unique_event_times = r(r)
	di `unique_event_times'
	
	// Find out which event time in the pre-period is closest to -1
	capture drop event_time_names1
	svmat event_time_names
	preserve 
		drop if event_time_names1 >= 0
		*drop if event_time_names1 < `low_event_cap'
		egen ref_event_time = max(event_time_names1)
		sum ref_event_time
		local ref_event_time = r(mean)
	restore
	di `ref_event_time'
	capture drop event_time_names1
	
	recode event_time_bacon (.=`ref_event_time')
	*ensure that "xi" omits -1
	char event_time_bacon[omit] `ref_event_time'
	
	
		xi i.event_time_bacon, pref(_T)
		
		
	// Determine number of unique "event times"		
	capture drop x_value 
	capture drop y_value 
	
	gen x_value = .
	gen y_value = . 
	
	local order_index = 1
	forvalues i = `low_event_cap'(1)`high_event_cap' {
		di "`i'"
		replace x_value = `i' in `order_index'
		unique `id_variable' if event_time_bacon == `i'
		replace y_value = `r(unique)' in  `order_index'
		
		local order_index = `order_index' + 1
	}
	
	graph bar y_value, over(x_value)
				
	// Save graph
	graph export "`output_file_name'-all.png", replace
	
	// Now for actual treated 
	drop if ever_treated == 0 
	
	capture drop x_value 
	capture drop y_value 
	
	gen x_value = .
	gen y_value = . 
	
	local order_index = 1
	forvalues i = `low_event_cap'(1)`high_event_cap' {
		di "`i'"
		replace x_value = `i' in `order_index'
		unique `id_variable' if event_time_bacon == `i'
		replace y_value = `r(unique)' in  `order_index'
		
		local order_index = `order_index' + 1
	}
	
	graph bar y_value, over(x_value)
			
	// Save graph
	graph export "`output_file_name'-only-treated.png", replace

end
