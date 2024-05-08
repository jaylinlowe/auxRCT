#' Create interactive Shiny app to calculate sample sizes using reloop
#'
#' @details `run_app` uses other `reloop.samp` functions to display a Shiny app that walks users
#' through the process of defining subgroups and calculating the sample size that would be
#' necessary if the RCT population resembled each subgroup. Sample sizes for each subgroup
#' without using ReLOOP are also displayed. Additionally, users can explore diagnostics
#' of the initial random forest results and can compare each subgroup to the rest of the
#' remnant for any covariate.
#'
#' @import dplyr
#' @import ggplot2
#' @import shiny
#' @import shinydashboard
#' @import shinyalert
#' @import scales
#'
#' @param ...
#'
#' @export
#'
run_app <- function(...) {
  ui <- dashboardPage(skin = "purple",
                      dashboardHeader(title = "Sample Size Calculations"),

                      dashboardSidebar(
                        collapsed = TRUE,
                        sidebarMenu()
                      ),
                      dashboardBody(
                        navbarPage("", id = "myNavbar",
                                   tabPanel("Upload Dataset",
                                            fluidPage(
                                              fileInput("dataset", "Choose CSV File",
                                                        multiple = TRUE,
                                                        accept = c("text/csv",
                                                                   "text/comma-separated-values,text/plain",
                                                                   ".csv")),

                                              checkboxInput("header", "Header", TRUE),
                                              helpText("Note: Please do not replace missing values with other characters in the data frame."),
                                              selectInput("out", "Select output column:", c('Upload dataset to see column choices'), multiple = F),
                                              actionButton(inputId = "button1", label = "Next: Choose Random Forest Variables")
                                            ))
                        )
                      )

  )

  server <- function(input, output, session) {

    #PAGE 1 - UPLOAD DATASET
    data <- reactive({
      req(input$dataset)
      df <- read.csv(input$dataset$datapath, header = input$header)
      return(df)
    })

    observeEvent(data(), {
      updateSelectInput(session, "out", choices=colnames(data()))
    })

    #PAGE 2 - INITIAL RANDOM FOREST PARAMETERS
    # set up second tab
    button1_counter <- reactiveVal(0)
    observeEvent(input$button1, {
      button1_counter(button1_counter()+1)
      if (button1_counter() == 1 ) {
        insertTab(inputId = "myNavbar",
                  tabPanel("Initial Random Forest Parameters",
                           helpText("Select variables to include in random forest: "),
                           actionLink("selectall","Select All"),
                           checkboxGroupInput("variables", "", c("Upload dataset to see column choices")),
                           actionButton(inputId = "button2", label = "Next: Run Random Forest")),
                  position = 'after',
                  target = "Upload Dataset"
        )
      }
      updateTabItems(session, "myNavbar", selected = "Initial Random Forest Parameters")
    })

    # update variable choices based on selected data set
    observeEvent(input$button1, {
      updateCheckboxGroupInput(session, "variables", choices=setdiff(colnames(data()), input$out))
    })

    # create select all option for variable check boxes
    observe({
      if(is.null(input$selectall)) return(NULL)
      else if (input$selectall%%2 == 0)
      {
        updateCheckboxGroupInput(session,"variables",choices=setdiff(colnames(data()), input$out))
      }
      else
      {
        updateCheckboxGroupInput(session,"variables",choices=setdiff(colnames(data()), input$out),selected=setdiff(colnames(data()), input$out))
      }
    })



    #PAGE 3 - RANDOM FOREST

    #setup for tab 3
    button2_counter <- reactiveVal(0)
    observeEvent(input$button2, {
      button2_counter(button2_counter()+1)
      if (button2_counter() == 1 ) {
        insertTab(inputId = "myNavbar",
                  tabPanel("Random Forest Variable Investigation",
                           helpText("This tab displays information about the variables included in the random forest.
                         It may be useful if you wish to pair down the variables in the original random forest or identify variables
                         that define subgroups of particular interest. If you do not wish to do so,
                          you can proceed directly to calculating sample sizes by clicking the button below:"),
                           actionButton(inputId = "button3", label = "Next: Calculate Sample Sizes"),
                           fluidRow(
                             column(6, tableOutput("var_table1")),
                             column(6, tableOutput("var_table2"))),
                           fluidRow(plotOutput("corr_plot")),
                           fluidRow(
                             selectInput("var_invest", "Variable to Investigate Further:", c("Run random forest to see column choices")),
                             plotOutput("plot"),
                             textOutput("na_text"))),
                  position = 'after',
                  target = "Initial Random Forest Parameters"
        )
      }
      updateTabItems(session, "myNavbar", selected = "Random Forest Variable Investigation")
    })

    #get table of counts of NA values per column prior to mean imputation
    na_table <- reactive({
      X <- data()[, !names(data()) == input$out]
      return(X %>% summarise(across(everything(), ~ sum(is.na(.)))))
    })

    #run initial random forest
    rf_results <- reactive({
      out_name <- input$out
      Y <- data()[[out_name]]
      X <- mean_imputation(data()[, !names(data()) == out_name])
      variables <- get_var_names(input$variables, X )
      #print(variables)
      #print(colnames(cbind(Y,X)))
      f <- as.formula(paste("Y ~", paste(variables, collapse = " + ")))
      rf.model <- randomForest::randomForest(f, data = cbind(Y, X))
      preds <- predict(rf.model)
      df <- cbind(Y, preds, X)
      return(list(df, rf.model))
    }) %>%
      bindEvent(input$button2)

    #results from random forest
    data_preds <- reactive({rf_results()[[1]]})
    rf_model <- reactive({rf_results()[[2]]})
    Y <- reactive({data_preds()[,1]})
    X <- reactive({data_preds()[,c(-1, -2)]})
    preds <- reactive({data_preds()[,2]})

    # get variable importance values
    diagnostics <- reactive({
      scores <- caret::varImp(rf_model()) %>%
        as.data.frame() %>%
        arrange(desc(Overall)) %>%
        head(20)

      imp_table <- data.frame(Variable = rownames(scores), `Importance` = scores$Overall)
      return(imp_table)
    })

    var_table_split <- reactive({importance_split(diagnostics())})
    output$var_table1 <- renderTable(var_table_split()[[1]])
    output$var_table2 <- renderTable(var_table_split()[[2]])

    # Correlation plot
    # corrs <- reactive({
    #   variables <- diagnostics()$Variable
    #   my_cor <- data_preds() %>%
    #     select(all_of(variables)) %>%
    #     select_if(is.numeric) %>%
    #     cor() %>%
    #     reshape2::melt() %>%
    #     ggplot(aes(x = Var1, y = Var2)) +
    #     geom_tile(aes(fill = value)) +
    #     geom_text(aes(label = round(value, 2))) +
    #     labs(x = "", y = "", title = "Correlations for Numerical Variables") +
    #     scale_fill_gradient2(
    #       low = muted("orange"),
    #       mid = "white",
    #       high = muted("navyblue"),
    #       midpoint = 0,
    #       limits = c(-1, 1)) +
    #     scale_x_discrete(limits = rev)
    #   return(my_cor)
    #
    # })

    output$corr_plot <- renderPlot(correlation_plot(X(), diagnostics()$Variable))

    #Variable Specific Investigation
    observeEvent(input$button2, {
      updateSelectInput(session, "var_invest", choices=diagnostics()$Variable)
    })

    var_plot <- reactive({
      var <- input$var_invest

      if (class(data_preds()[[var]]) == "numeric") {
        p <- ggplot(data_preds()) + geom_histogram(aes(x = get(var))) +
          labs(x = aes_string(var))
      }
      else {
        p <- ggplot(data_preds()) + geom_bar(aes(x = get(var))) +
          labs(x = aes_string(var))
      }
      return(p)
    })

    output$plot <- renderPlot(var_plot())


    #Get NA count for chosen variable
    na_text <- reactive({
      var <- input$var_invest
      count_na <- na_table()[[var]]
      return(paste("Number of NA Observations: ", count_na))
    })


    #PAGE 4 - INITIAL RESULTS
    button3_counter <- reactiveVal(0)
    observeEvent(input$button3, {
      button3_counter(button3_counter()+1)
      if (button3_counter() == 1 ) {
        insertTab(inputId = "myNavbar",
                  tabPanel("Initial Results",
                           fluidRow(
                             column(3,
                                    wellPanel(
                                      fluidRow(radioButtons("effect_options", label = "Effect Size Units", choices = c("Standard Deviations of Outcome", "Raw Number"))),
                                      fluidRow(numericInput("effect_size", "Effect Size: ", value = 0.2)),
                                      fluidRow(numericInput("alpha", "Alpha: ", 0.05, min = 0.01, max = 0.1)),
                                      fluidRow(numericInput("beta", "Beta: ", value = 0.20, min = 0.01, max = 0.30)),
                                      fluidRow(actionButton("run_everybody", "Calculate General Sample Sizes", style="background-color: #5F5CA3; color: #FFFFFF")),
                                      fluidRow(actionLink("moreinfo_params", "Click here for more information about these inputs"))
                                    )
                             ),
                             column(9,
                                    valueBoxOutput("ss_without_reloop"),
                                    valueBoxOutput("ss_with_reloop"),
                             )
                           ),
                           actionButton("button4", "Next: Calculate Sample Sizes for Subgroups")),
                  position = 'after',
                  target = "Random Forest Variable Investigation"
        )
      }
      updateTabItems(session, "myNavbar", selected = "Initial Results")
    })

    effect_size <- reactive({
      if (input$effect_options == "Standard Deviations of Outcome") {
        return(input$effect_size * sd(Y()))
      }
      else {
        return(input$effect_size)
      }
    })
    alpha <- reactive({input$alpha})
    beta <- reactive({input$beta})

    # Get overall estimate using reloop
    samp_size_reloop <- reactive({
      #overall_mse <- sum((Y() - preds())^2)/length(Y())
      resid_var <- var(Y()-preds())
      return(samp_size(resid_var, effect_size = effect_size(), alpha = alpha(), beta = beta()))
    }) %>%
      bindEvent(input$run_everybody)

    output$ss_with_reloop <- renderValueBox({
      valueBox(
        value = round(samp_size_reloop()) ,
        subtitle = "Recommended overall sample size if entire remnant is representative of RCT, with reloop",
        color = "navy"
      )
    })

    # Get overall estimate without reloop
    samp_size_without <- reactive({
      sigma2 <- var(Y())
      return(samp_size(sigma2, effect_size = effect_size(), alpha = alpha(), beta = beta()))
    }) %>%
      bindEvent(input$run_everybody)

    output$ss_without_reloop <- renderValueBox({
      valueBox(
        value = round(samp_size_without()),
        subtitle = "Recommended overall sample size if entire remnant is representative of RCT, Without Reloop",
        color = "navy"
      )
    })

    # PAGE 5 - SUBGROUP RESULTS
    alpha2 <- reactive({input$alpha}) %>% bindEvent(input$button4)
    beta2 <- reactive({input$beta}) %>% bindEvent(input$button4)
    effect_size2 <- reactive({input$effect_size}) %>% bindEvent(input$button4)

    button4_counter <- reactiveVal(0)
    observeEvent(input$button4, {
      button4_counter(button4_counter()+1)
      if (button4_counter() == 1 ) {
        insertTab(inputId = "myNavbar",
                  tabPanel("Subgroup Results",
                           fluidRow(
                             column(6,
                                    wellPanel(
                                      fluidRow(radioButtons("effect_options2", label = "Effect Size Units", choices = c("Standard Deviations of Outcome", "Raw Number"))),
                                      fluidRow(numericInput("effect_size2", "Effect Size: ", value = effect_size2())),
                                      fluidRow(numericInput("alpha2", "Alpha: ", alpha2(), min = 0.01, max = 0.1)),
                                      fluidRow(numericInput("beta2", "Beta: ", value = beta2(), min = 0.01, max = 0.30)),
                                      fluidRow(actionLink("moreinfo_params2", "Click here for more information about these inputs"))

                                    )),
                             column(6,
                                    wellPanel(
                                      fluidRow(selectInput("method", "Which method would you like to use to select subgroups?",
                                                           choices = c('Defined by Factor Variable', 'Defined by Numeric Variable', 'Best-Worst Case Scenarios'), multiple = F)),
                                      fluidRow(conditionalPanel(condition = "input.method == 'Defined by Factor Variable'",
                                                                selectInput("grouping", "Select grouping variable (must be character or factor variable): ", c('Upload dataset to see column choices'), multiple = F)
                                      )),
                                      fluidRow(conditionalPanel(condition = "input.method == 'Best-Worst Case Scenarios'",
                                                                numericInput("num_groups", "How many subgroups?", value = 10)
                                      )),
                                      fluidRow(conditionalPanel(condition = "input.method == 'Defined by Numeric Variable'",
                                                                selectInput("num_var", "Select a numeric variable to split on: ", c('Upload dataset to see column choices'), multiple = F),
                                                                numericInput("max_groups", "Maximum Number of Subgroups", value = 10)
                                      )),
                                      fluidRow(actionButton("run_results", "Calculate Sample Sizes for Subgroups", style="background-color: #5F5CA3; color: #FFFFFF")),
                                      fluidRow(actionLink("moreinfo", "Click here for more information about subgroup methods"))
                                    )
                             )
                           ),
                           fluidRow(
                             tableOutput(outputId = "results")
                           ),
                           actionButton("button5", "Next: Investigate Subgroups of Interest")),
                  position = 'after',
                  target = "Initial Results"
        )
      }
      updateTabItems(session, "myNavbar", selected = "Subgroup Results")
    })

    observeEvent(input$moreinfo_params, {
      shinyalert(title = "", type = "info", html = T,
                 text = tagList(
                   h3("Effect Size"),
                   p("The minimum size that you want to detect. This can be entered either as a raw number,
                   or in terms of the standard deviation of the population. Note that the second option is always
                   in terms of the standard deviation of the population, even when calculating for subgroups."),
                   h3("Alpha"),
                   p("Alpha is the desired Type I error rate, or the probability of rejecting the null hypothesis when it is true. Alpha is typically set to 0.05. "),
                   h3("Beta"),
                   p("Beta represents the Type II error rate, or the probability of failing to reject the null hypothesis when the null hypothesis is false.
                  Beta is typically set at 0.10 or 0.20.")
                 ))
    })

    observeEvent(input$moreinfo_params2, {
      shinyalert(title = "", type = "info", html = T,
                 text = tagList(
                   h3("Effect Size"),
                   p("The minimum size that you want to detect. This can be entered either as a raw number,
                   or in terms of the standard deviation of the population. Note that the second option is always
                   in terms of the standard deviation of the population, even when calculating for subgroups."),
                   h3("Alpha"),
                   p("Alpha is the desired Type I error rate, or the probability of rejecting the null hypothesis when it is true. Alpha is typically set to 0.05. "),
                   h3("Beta"),
                   p("Beta represents the Type II error rate, or the probability of failing to reject the null hypothesis when the null hypothesis is false.
                  Beta is typically set at 0.10 or 0.20.")
                 ))
    })

    observeEvent(input$moreinfo, {
      shinyalert(title = "", type = "info", html = T,
                 text = tagList(
                   h3('Best-Worst Case Scenarios'),
                   p("This method splits observations into subgroups based on predicted error. Specifically, we use the initial random forest to cacluate the
                   absolute value of the error for each observation. We run a second random forest, attempting to predict the errors using all other variables.
                   Observations are then split into deciles (or a different number) based on the magnitude of their predicted error."),
                   h3("Defined by Factor Variable"),
                   p("This method defines subgroups by a single factor variable. One subgroup is created per level of the factor variable.
                 This method can also be used if you wish to define subgroups by more complex decision rules involving multiple variables.
                 In that case, create a factor variable according to the desired splits and re-upload the dataset with the new variable."),
                   h3("Defined by Numeric Variable"),
                   p("This method defines subgroups by splitting up a numeric variable into equally sized groups (default: 10).
                 However, if a variable has a large proportion of observations with an identical value, this will fail.
                 In that case, observations with the identical value will form one subgroup, and the remaining portion of observations
                 will be divided into 9 (or 1 minus the specified value) equally sized groups. If there are not enough observations left
                to divide into 9 equal groups, all other observations that do not have the identical value will form one subgroup.
                  Note: if you wish to create subgroups according to two or more numeric variables, you will need to create a factor variable
                  accordingly and use the previous method instead.")
                 ))
    })


    observeEvent(input$button4, {
      updateNumericInput(session, "effect_size2", value = effect_size2())
    })

    check_type <- function(x) {
      if (is.character(x) | is.factor(x)) {
        return(TRUE)
      }
      return(FALSE)
    }

    observeEvent(input$button4, {
      updateSelectInput(session, "grouping", choices=setdiff(names(data())[ sapply(data(), check_type)], input$out))
    })

    #should display numeric only
    observeEvent(input$button4, {
      updateSelectInput(session, "num_var", choices=setdiff(names(data())[ sapply(data(), is.numeric)], input$out))
    })

    method <- reactive({input$method})

    samp_table <- reactive({
      if (method() == 'Defined by Factor Variable') {
        grouping_col <- input$grouping
        validate(need((class(X()[[grouping_col]]) %in% c("factor", "character")), "Please select a factor or character column for grouping."))
        results <- reloop.samp(Y(), X(), grouping_col, preds(), effect_size(), alpha(), beta())
      }
      else if (method() == 'Best-Worst Case Scenarios') {
        subgroups <- error_subgroups(Y(), X(), preds(), num_groups = input$num_groups)
        X_temp <- cbind(X(), subgroups)
        results <- reloop.samp(Y(), X_temp, grouping_col = "subgroups", preds(), effect_size(), alpha(), beta())
      }
      else if (method() == "Defined by Numeric Variable") {
        max_groups <- input$max_groups
        grouping_col <- input$num_var
        validate(need((class(X()[[grouping_col]]) %in% c("numeric")), "Please select a numeric column"))
        #X()$subgroups <- numeric_subgroups(Y(), X(), grouping_col, preds(), max_groups)
        #print(X()$subgroups)

        subgroups <- numeric_subgroups(Y(), X(), grouping_col, preds(), max_groups)
        X_temp <- cbind(X(), subgroups)
        results <- reloop.samp(Y(), X_temp, grouping_col = "subgroups", preds(), effect_size(), alpha(), beta())
      }
      return(results)
    }) %>%
      bindEvent(input$run_results)

    renamed_table <- reactive({
      results <- samp_table()[[1]]
      results <- results %>%
        select(def, num, samp_size, samp_size_without) %>%
        rename("Subgroup Definition" = "def", "Recommended Sample Size With ReLoop" = "samp_size",
               "Number of Observations in Subgroup" = "num", "Recommended Sample Size Without ReLoop" = "samp_size_without") %>%
        relocate(`Recommended Sample Size With ReLoop`, .after = `Recommended Sample Size Without ReLoop`)
      return(results)
    })

    output$results <- renderTable(renamed_table(), digits = 0)


    #PAGE 6 - SUBGROUP INVESTIGATION
    button5_counter <- reactiveVal(0)
    observeEvent(input$button5, {
      button5_counter(button5_counter()+1)
      if (button5_counter() == 1 ) {
        insertTab(inputId = "myNavbar",
                  tabPanel("Subgroup Investigation",
                           fluidRow(
                             column(12,
                                    wellPanel(
                                      selectInput("subgroup", "Subgroup definition (from results table)", c("Run results to see column choices")),
                                      selectInput("subgroup_var", "Variable to compare", choices = c("Upload dataset to see choices")),
                                      plotOutput("sub_plot"),
                                      plotOutput("miss_plot")
                                    )
                             )
                           )),
                  position = 'after',
                  target = "Subgroup Results"
        )
      }
      updateTabItems(session, "myNavbar", selected = "Subgroup Investigation")
    })

    observeEvent(input$button5, {
      updateSelectInput(session, "subgroup", choices=samp_table()[[1]]$def)
    })

    observeEvent(input$button5, {
      updateSelectInput(session, "subgroup_var", choices=c(colnames(data())))
    })


    subgroup_def <- reactive({
      subgroup_def <- input$subgroup
      subgroup_def_index <- which(samp_table()[[1]]$def == subgroup_def)
      subgroup_index <- samp_table()[[2]][[subgroup_def_index]]

      vec <- seq(1, nrow(X()))
      subgroup_vec <- ifelse(vec %in% subgroup_index, "Subgroup", "General Population")
      return(subgroup_vec)
    })

    #TO DO - DOES THE GENERAL POPULATION ALSO INCLUDE THE SUBGROUP????
    #YES - DID I FIX THAT FOR THE NA PLOT TOO???
    #graphs comparing distribution of subgroup to rest of population
    #can this go into a function?
    # subgroup_tab <- reactive({
    #   df <- cbind(Y(), subgroup_def(), X())
    #   extra_subgroup <- filter(df, `subgroup_def()` == "Subgroup")
    #   extra_subgroup$`subgroup_def()` <- rep("General Population", nrow(extra_subgroup))
    #   df <- rbind(df, extra_subgroup)
    #
    #   var <- input$subgroup_var
    #
    #   if (var == input$out) {
    #     var <- "Y()"
    #     df[[var]] <- as.numeric(df[[var]])
    #   }
    #   if (class(df[[var]]) == "numeric") {
    #     p <- ggplot(df, aes(x = .data[[var]], fill = `subgroup_def()`)) +
    #       geom_boxplot() +
    #       labs(fill = "") +
    #       theme(legend.position = "bottom") +
    #       scale_fill_manual(values = c(rgb(95, 92, 163, max = 255), rgb(36, 45, 49, max = 255)))
    #   }
    #   else {
    #     p <- df %>%
    #       group_by(`subgroup_def()`, .data[[var]]) %>%
    #       summarize(n = n()) %>%
    #       group_by(`subgroup_def()`) %>% #group no longer exists here
    #       mutate(prop = n/sum(n)) %>%
    #       ggplot(aes(y = prop, fill = `subgroup_def()`, x = .data[[var]])) +
    #       geom_bar(position = "dodge", stat = "identity") +
    #       labs(fill = "", y = "Proportion of Group") +
    #       coord_flip() +
    #       theme(legend.position = "bottom") +
    #       scale_fill_manual(values = c(rgb(95, 92, 163, max = 255), rgb(36, 45, 49, max = 255)))
    #   }
    #   return(p)
    # })
    #
    # output$sub_plot <- renderPlot(subgroup_tab())
    output$sub_plot <- renderPlot(subgroup_plot(Y(), X(), subgroup_def(), var = input$subgroup_var, out = input$out))

    #density plot of number of NA values per observation for the subgroup compared to the rest of the population
    # miss_plot <- reactive({
    #   na_count <- rowSums(is.na(data()))
    #   if (sum(na_count) > 0) {
    #     na_df <- data.frame(cbind(na_count, subgroup_def()))
    #     na_df$na_count <- as.numeric(na_df$na_count)
    #     colnames(na_df) <- c("na_count", "Group Definition")
    #     p <- ggplot(na_df, aes(x = na_count, fill = `Group Definition`)) +
    #       geom_density(stat = "count", alpha = 0.5) +
    #       scale_fill_manual(values = c(rgb(95, 92, 163, max = 255), rgb(36, 45, 49, max = 255))) +
    #       labs(fill = "", x = "Number of NA Values", y = "Count") +
    #       theme(legend.position = "bottom") +
    #       ggtitle("Distribution of Number of NA Values per Observation")
    #   }
    #   else {p <- NULL}
    #   return(p)
    # })

    output$miss_plot <- renderPlot(na_plot(data(), subgroup_def()))

  }


  shinyApp(ui = ui, server = server)
}
