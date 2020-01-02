// state is treated units
// time is time variable
// time_treated is date of treatment 
// treated is ever_treated x post_treatment


capture program drop create_event_study_graph
program create_event_study_graph, 
syntax, ///
	first_time(real) ///
	last_time(real) ///
	first_treated_time(real) ///
	last_treated_time(real) ///
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
	output_file_name(string) ///
	time_variable(string) ///
	id_variable(string) ///
	weight_var(string)

	////////////////////////////////////////////////////////////////////////////
	// Set up variables
	capture drop weight
	if "`weight_var'" == "none" {
		gen weight = 1
	}
	if "`weight_var'" != "none" {
		gen weight = `weight_var'
	}
	
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

	
	if `unique_event_times' < (`high_event_cap' - `low_event_cap' + 1) {
		local gap = 1
		di in red "There are gaps in treatment event time."
		di in red "The event times are `r(levels)'"
		di in red "The reference time is `ref_event_time'"

	}
	if `unique_event_times' == (`high_event_cap' - `low_event_cap' + 1) {
		local gap = 0
		di in red "There are no gaps in treatment event time."

	}	
	
	*mak dummies
	xi i.event_time_bacon, pref(_T)
	/*
	* Position of -2
	if `low_event_cap' > ( `first_time' - `last_treated_time') {
	 local pos_of_neg_2 = abs(2 - (-`low_event_cap' + 1))
	}	

	if `low_event_cap' <= ( `first_time' - `last_treated_time') {
	 local pos_of_neg_2 = abs(2 - (-( `first_time' - `last_treated_time') + 1))
	}

	* Position of 0
	local pos_of_zero = `pos_of_neg_2' + 2

	* Position of max
	if `high_event_cap' >= ( `last_time' - `first_treated_time') {
	 local pos_of_max = `last_time' - `first_treated_time' + `pos_of_zero'
	}	

	if `high_event_cap' < ( `last_time' - `first_treated_time') {
	 local pos_of_max = `high_event_cap' + `pos_of_zero'
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Run DD regression
	reghdfe `y_var' ///
		treated ///
		`dd_control_list' ///
		`control_list' ///
		[aw = weight] ///
	, absorb(`fixed_effects') cluster(`cluster_var')
	*/
	// Store DD Coefs and 95% CI and p-value

	
	/////////////////////////////////////////////////////////////////////////////
	// Run event-study regression
	reghdfe `y_var' ///
		_T* ///
		`control_list' ///
		`dd_control_list' ///
		[aw = weight] ///
		, absorb(`fixed_effects') cluster(`cluster_var')
		
		quietly {		
			*if `gap' == 1 {
				// https://www.statalist.org/forums/forum/general-stata-discussion/general/1466327-trying-to-get-the-square-root-of-the-diagonal-of-e-v

				*mata st_matrix("srdvcovbt",sqrt(diagonal(st_matrix("e(V)"))))
				*matrix list srdvcovbt
				// List of coefficients, 95% CI, and event time
				mat b = e(b)'
				mat list b 
				mata st_matrix("b",select(st_matrix("b"),st_matrix("b")[.,1]:!=0))
				mat list b

				mata st_matrix("ll", (st_matrix("e(b)"))'- invt(st_numscalar("e(df_r)"), 0.975)*sqrt(diagonal(st_matrix("e(V)"))))
				mat list ll
				mata st_matrix("ll",select(st_matrix("ll"),st_matrix("ll")[.,1]:!=0))
				mat list ll
				
				mata st_matrix("ul", (st_matrix("e(b)"))'+ invt(st_numscalar("e(df_r)"), 0.975)*sqrt(diagonal(st_matrix("e(V)"))))
				mat list ul
				mata st_matrix("ul",select(st_matrix("ul"),st_matrix("ul")[.,1]:!=0))
				mat list ul		
				
				*local ref_event_time = -1
				
				// Make an event time list that is in the same order as the coefs. 
				levelsof(event_time_bacon) if event_time_bacon != `ref_event_time' & !missing(`y_var') , matrow(event_time_names_without_ref)
				mat list event_time_names_without_ref
				
				// Find out max order
				capture drop event_time_names_without_ref1
				svmat event_time_names_without_ref
				preserve 
					gen temp_order = _n if !missing(event_time_names_without_ref1)
					sum temp_order
					local max_pos = r(max)
				restore
				di `max_pos'
				capture drop event_time_names_without_ref1		
				
				
				// Combine event time, coef, lower ci, and upper ci into one matrix
				mat results =  event_time_names_without_ref, b[1..`max_pos',.], ll[1..`max_pos',.], ul[1..`max_pos',.]
				mat list results
			
			// Drop rows bigger or smaller than cap 
				*Fix hard coding of 10
				mata st_matrix("results_new",select(st_matrix("results"),st_matrix("results")[.,1]:>=`low_event_cap_graph'))
				mata st_matrix("results_new",select(st_matrix("results_new"),st_matrix("results_new")[.,1]:<=`high_event_cap_graph'))
				mat list results_new
				
				// Make an event time list without referecen and without caps. 

				
				levelsof(event_time_bacon) if event_time_bacon != `ref_event_time' & event_time_bacon >= `low_event_cap_graph' & event_time_bacon <= `high_event_cap_graph'  & !missing(`y_var'), matrow(event_time_names_wo_ref_capped)
				mat list event_time_names_wo_ref_capped		
				
				// Find out order of event time in the pre-period is right before reference event-time
				capture drop event_time_names_wo_ref_capped1
				svmat event_time_names_wo_ref_capped
				preserve 
					drop if event_time_names_wo_ref_capped1 >= 0
					*drop if event_time_names_without_ref1 < -10
					egen time_b4_ref_event_time = max(event_time_names_wo_ref_capped1)
					gen temp_order = _n if !missing(event_time_names_wo_ref_capped1)
					sum temp_order if time_b4_ref_event_time == event_time_names_wo_ref_capped1
					local time_b4_ref_event_time = r(mean)
				restore
				*di `time_b4_ref_event_time'
				capture drop event_time_names_wo_ref_capped1
				
				// Find out order of event time in the post-period is right after reference event-time
				capture drop event_time_names_wo_ref_capped1
				svmat event_time_names_wo_ref_capped
				preserve 
					*drop if event_time_names_without_ref1 < 0
					*drop if event_time_names_without_ref1 > 10
					egen time_after_ref_event_time = min(event_time_names_wo_ref_capped1) if event_time_names_wo_ref_capped1 >= 0
					gen temp_order = _n if !missing(event_time_names_wo_ref_capped1)
					sum temp_order if time_after_ref_event_time == event_time_names_wo_ref_capped1
					local time_after_ref_event_time = r(mean)
					
					sum temp_order
					local max_pos = r(max)
				restore
				*di `time_after_ref_event_time'
				*di `max_pos'
				capture drop event_time_names_wo_ref_capped1		
				

				// Add row where reference group should be
				mat results_new = results_new[1..`time_b4_ref_event_time',.] \ `ref_event_time', 0, 0, 0 \ results_new[`time_after_ref_event_time'..`max_pos',.]
				mat list results_new
				
				
				
				capture drop order
				capture drop b 
				capture drop high 
				capture drop low
					
				svmat results_new 
				
				rename results_new1 order
				rename results_new2 b
				rename results_new3 low
				rename results_new4 high
		*	}
		}
		
		/*
	quietly {		
		if `gap' == 0 {
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
		}
	}
	*/
	
	// Plot estimates	
	twoway rarea low high order if order <= `high_event_cap_graph' & order >= `low_event_cap_graph' , fcol(gs14) lcol(white) msize(3) /// estimates
		|| connected b order if order <= `high_event_cap_graph' & order >= `low_event_cap_graph', lw(1.1) col(white) msize(7) msymbol(s) lp(solid) /// highlighting
		|| connected b order if order <= `high_event_cap_graph' & order >= `low_event_cap_graph', lw(0.6) col("71 71 179") msize(5) msymbol(s) lp(solid) /// connect estimates
		|| scatteri 0 `low_event_cap_graph' 0 `high_event_cap_graph', recast(line) lcol(gs8) lp(longdash) lwidth(0.5) /// zero line 
			xlab(`low_event_cap_graph'(`x_lab_step_size')`high_event_cap_graph' ///
					, nogrid labsize(7) angle(0)) ///
			ylab(, nogrid labs(7)) ///
			legend(off) ///
			xtitle("`x_label'", size(4)) ///
			subtitle("`y_label'", size(4) pos(11)) ///
			title("`graph_label'", size(5) pos(11)) ///
			xline(-.5, lpattern(dash) lcolor(gs7) lwidth(1)) 	
			
	// Save graph
	graph export "`output_file_name'.png", replace

end
