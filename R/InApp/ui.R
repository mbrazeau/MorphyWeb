shinyUI(fluidPage(
  titlePanel("Inapplicable data reconstruction"),



  # Create a row for additional information
  fluidRow(
      ## Tree input
      column(3,
        radioButtons("tree", label = h3("Tree input method"), choices = list("Random" = 1, "User" = 2), selected = 1),
        conditionalPanel(
          condition = "input.tree == 1",
            selectInput("tree_type", label = "Tree topology type", choices = list("Random", "Balanced", "Left", "Right", "Left-Right"), selected = "Random"),
            sliderInput("n_taxa", label = "Number of taxa:", min = 3, max = 100, value = 12),
            helpText("Note: if the tree type 'Balanced', the number of taxa must be a power of 2 (2,4,8, ...); if the tree type is `Left-Right` the number of taxa must be even. In both cases, if the number of taxa does not match, a random tree is used instead.")
        ),
        conditionalPanel(
          condition = "input.tree == 2",
            textInput("newick_tree", label = h5("Enter a newick tree:"), value = "((a,b),(c,d));")
        )
      ),
      
      ## Character input
      column(4,
        radioButtons("character", label = h3("Character input method"), choices = list("Random" = 1, "User" = 2), selected = 1),
        conditionalPanel(
          condition = "input.character == 2",
            textInput("character_string", label = h5("Enter a character string:"), value = "1?2-"),
            helpText("Note: the number of characters must match the size of the tree! Accepted states are any values from 0 to 9, - for the inapplicable token and ? for all states (missing data).")
        ),
        checkboxGroupInput("showPass", label = h3("Show passes"),  choices = list("1st Downpass" = 1, "1st Uppass" = 2, "2nd Downpass" = 3, "2nd Uppass" = 4), selected = c(1,2,3,4)),
        helpText("Tick the passes to be displayed on the nodes.")
      )
  ),

  
  fluidRow(
    ## Plots the algorithm results
    plotOutput("plot_out", width = "80%", "1000px")
  )

))

    # ## DEBUG
    # mainPanel(
    #   plotOutput("plot_out", width = "100%"),
    #   textOutput("tree"),
    #   textOutput("character"),
    #   textOutput("character_string"),
    #   textOutput("tree_type"),
    #   textOutput("n_taxa"),
    #   textOutput("newick_tree"),
    #   textOutput("showPass")
    # )