## Clear memory
rm(list=ls())

## Set-up

### Load packages

if (!require("pacman")) {
    install.packages("pacman")
}
pacman::p_load(
    haven, tidyverse, data.table, lubridate, zoo, fixest, knitr, kableExtra, lazyeval
)

### Set output options

options("tidylog.display" = NULL)

### Read in data for analysis

divorce_data <- haven::read_dta("~/Documents/GitHub/event_study_graphs/Til Death MASTER DATA streamlined.dta")

### Manipulate data so it matches the stata program I made. 

divorce_data <- divorce_data %>%
    filter(sex == 2) %>% #  Keep only women
    mutate(ever_treated = if_else(is.na(divyear), 0, 1)) %>% # Generate treated
    mutate(post = if_else(year >= divyear & !is.na(divyear), 1, 0)) %>% # Generate post
    mutate(treated = ever_treated*post)  %>%
    rename(time_treated = divyear)


### Create event-time variables and dummies

event_time <- function(dataset, id, time, treated, high_cap, low_cap) {
    
    # Create a temporary dataset that creates "event-time" and caps this variable. 
    # Never treated units get set to -1
    
    temp_dataset <- dataset %>% 
        group_by({{ id }}) %>% 
        filter({{ treated }} == 1) %>%
        slice(which.min({{ time }})) %>%
        mutate(first_time = {{ time }}) %>%
        select({{ id }}, first_time) %>%
        right_join(dataset, by = c(.data[[id]])) %>%
        mutate(event_time = {{ time }} - first_time) %>%
        mutate(event_time_capped = case_when(
            event_time > high_cap ~ high_cap + 1, 
            event_time < low_cap ~ low_cap - 1, 
            is.na(event_time) ~ -1, 
            event_time >= low_cap & event_time <= high_cap ~ event_time
        )) %>%
        select(event_time, event_time_capped, {{ id }}, {{ time }})
    
    dataset <- dataset %>% 
        left_join(temp_dataset, by = c(.data[[id]], .data[[year]]))
    
    # Create a list of unique levels of the capped event-time variable
    unique_vals <- levels(as.factor(dataset$event_time_capped))
    
    # Remove -1 event-time (so it will be forced to be reference group)
    unique_vals <- unique_vals[!(unique_vals == -1)]
    
    # Force the dataset to be a data.table
    dataset = data.table(dataset)
    
    # Create a dummy variable matrix. 
        # Note this is strictly binary 
        # You could imagine wanting NA values 
            # (e.g., can you know if T-5 is really zero or one if the current year is 2020 and the data year is 2019?)
            # Technically this should be NA IMO. 
    
    for (i in unique_vals) {
        dataset[event_time_capped == i, paste0('_T', i) := 1]
    }
    
    for (i in seq(ncol(dataset)-length(unique_vals)+1,ncol(dataset))) set(dataset, i=which(is.na(dataset[[i]])), j=i, value=0)
    
    # Output is the dataset
    return(dataset)
}

# Create event time dataset
divorce_data_event_time <- event_time(dataset = divorce_data, 
           id = st, 
           time = year, 
           treated = treated, 
           low_cap = -10, 
           high_cap = 20)

# Function to give preview of structure in dataset. 
display_structure <- function(dataset, id, time, treated, rowstart = 1, rowend = 5) {
    dataset %>%
        group_by(.data[[id]], .data[[time]]) %>%
        slice(1) %>%
        ungroup() %>%
        select(c(.data[[id]], .data[[time]], .data[[treated]], event_time_capped, starts_with("_T"))) %>%
        slice(rowstart:rowend) %>%
        kbl() %>%
        kable_styling()  
}

# Run function to display structure
display_structure(
    dataset = divorce_data_event_time, 
    id = "st", 
    time = "year", 
    treated = "treated", 
    rowstart = 10, 
    rowend = 25)


## Run DID regression

### First. Run a hard coded version of the diff-in-diff

# Note: The wildcard * in stata seems waaaaay superior to this verbose format needed for R. 
# There must be an easier way of doing this. 

feols(
    as.formula(
        paste0(
            "suiciderate_elast_jag ~ treated + ",
            paste(grep("^sa", names(divorce_data), value=TRUE), collapse = " + "), 
            " + ",
            paste(grep("^sh", names(divorce_data), value=TRUE), collapse = " + "),
            " + afdcmax + afdcrolls + roevwade + femep + stateur + lnpersinc", 
            " | st + year")),  divorce_data, cluster = ~st)


### Function for making did table

did_regression <- function(dataset, id, time, treated, y_var, dd_control_list = NULL, control_list = NULL, fixed_effects, cluster_var){
    
    # Format any lists. 
    if (!is.null(dd_control_list)) {
        dd_control_list_formatted <- paste0(" + ", paste(dd_control_list, collapse = " + "))
    }
    if (is.null(dd_control_list)) {
        dd_control_list_formatted <- ""
    }
    
    if (!is.null(control_list)) {
        control_list_formatted <- paste0(" + ", paste(control_list, collapse = " + "))
    }
    if (is.null(control_list)) {
        control_list_formatted <- ""
    }
    

    # Estimate model 
    model <- eval(
                substitute(
                    feols(as.formula(paste0(y_var, 
                                " ~ ", 
                                treated, 
                                dd_control_list_formatted,
                                control_list_formatted,
                                     " | ", 
                                fixed_effects
                                )
                         ), 
                   data = dataset,
                   cluster = ~cluster_var
                 )
                 )
                )



    # Display estimates
    table <- etable(model, 
                    subtitles = c("1"), 
                    keep = c(treated), 
                    digits = 2)
    
    return(
        table %>% 
        kbl() %>%
        kable_styling()
    )
}

### Run DID table. Later I want to add a feature to extract the treatment effect estimate. 

did_regression(dataset = divorce_data_event_time,
                id = st,
                time = year,
                treated = "treated",
                y_var = "suiciderate_elast_jag",
                dd_control_list = NULL,
                control_list = c("afdcmax", "afdcrolls", "roevwade", "femep", "stateur", "lnpersinc"),
                fixed_effects = "st + year",
                cluster_var = st)


## Run event-study regression
event_study_regression <- function(dataset, id, time, treated, y_var, dd_control_list = NULL, control_list = NULL, fixed_effects, cluster_var, high_cap, low_cap){
    
    # Format any lists. 
    if (!is.null(dd_control_list)) {
        dd_control_list_formatted <- paste0(" + ", paste(dd_control_list, collapse = " + "))
    }
    if (is.null(dd_control_list)) {
        dd_control_list_formatted <- ""
    }
    
    if (!is.null(control_list)) {
        control_list_formatted <- paste0(" + ", paste(control_list, collapse = " + "))
    }
    if (is.null(control_list)) {
        control_list_formatted <- ""
    }
    
    # Estimate model 
    model <- eval(
        substitute(
            feols(as.formula(paste0(y_var, 
                                    " ~ `",
                                    paste(grep("_T", names(dataset), value=TRUE), 
                                          collapse = "` + `"), 
                                    "`", 
                                    dd_control_list_formatted,
                                    control_list_formatted,
                                    " | ", 
                                    fixed_effects
            )
            ), 
            data = dataset,
            cluster = ~cluster_var
            )
        )
    )
    ests <- 
        broom::tidy(model)
    
    ests <- ests %>%
        add_row(term = "`_T-1`", estimate = 0, std.error = 0, .before = abs(low_cap) + 1)
    
    ests <- ests %>%
        slice(1:(abs(low_cap) + high_cap + 3))
    
    ests$x <- c(seq(low_cap - 1, high_cap + 1))

    # Display estimates
    table <- etable(model, 
                    subtitles = c("1"), 
                    digits = 2)
    
    table <- table %>% 
        kbl() %>%
        kable_styling()
    
    output <- list(ests, table)
    names(output) <- c("ests", "table")
    return(output)
        
    
}




event_study_results <- event_study_regression(dataset = divorce_data_event_time,
               id = st,
               time = year,
               treated = "treated",
               y_var = "suiciderate_elast_jag",
               dd_control_list = NULL,
               control_list = c("afdcmax", "afdcrolls", "roevwade", "femep", "stateur", "lnpersinc"),
               fixed_effects = "st + year",
               cluster_var = st, 
               low_cap = -10,
               high_cap = 20)


event_study_results$ests
event_study_results$table

# Make a function that plots the event-study results
plot_event_study <- function(event_study_results, high_cap, low_cap, breaks, plot_title, plot_subtitle, x_title){
   
     # Drop the caps
    event_study_results$ests <- event_study_results$ests %>%
        filter(x != high_cap + 1 & x != low_cap - 1)
    
    event_study_plot <- ggplot(event_study_results$ests) +
        geom_hline(yintercept = 0, color = "grey50", size = 1, linetype = "dashed") +
        geom_vline(xintercept = -.5, color = "grey25", size = 2, linetype = "longdash") +
        geom_ribbon(aes(x = x, ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error), alpha = 0.2) +
        geom_line(aes(x = x, y = estimate, size = 1), color = "red", size = 1) + 
        geom_point(aes(x = x, y = estimate, size = 1), color = "white", size = 6, shape = 15) + 
        geom_point(aes(x = x, y = estimate, size = 1), color = "red", size = 5, shape = 15) + 
        theme_minimal() +
        theme(
            legend.position = "none",
            axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24),
            axis.title.x = element_text(size = 24), axis.title.y = element_text(size = 24),
            panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),
            axis.line = element_line(colour = "black")
        ) +
        scale_x_continuous(breaks = seq(from = low_cap, to = high_cap, by = breaks)) +
        labs(
            x = x_title,
            y = "",
            title = plot_title,
            subtitle = plot_subtitle
        )
    
    return(event_study_plot)
}

event_study_plot <- plot_event_study(event_study_results = event_study_results,
                                     low_cap = -10,
                                     high_cap = 20,
                                     breaks = 5, 
                                     plot_title = "Change in suicide rate by years since law change", 
                                     x_title = "Years since unilateral divorce law", 
                                     plot_subtitle = "Aggregate suicide rate")



ggsave(filename = "~/Documents/GitHub/event_study_graphs/event-study-wolfers-in-R.pdf", 
       plot = event_study_plot, 
       width = 8, 
       height = 6)


