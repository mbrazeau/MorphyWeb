library(shiny)
library(ape)

# server.R
shinyServer(
  function(input, output) {

    source("helpers.R")
  
    # output$method <- renderText({ 
    #   is.numeric(input$method)
    #   if(as.numeric(input$method) == 1) {
    #     print("Should be fitch inapplicable")
    #     input$showPassInapp
    #   } 
    #   if(as.numeric(input$method) == 2) {
    #     print("Should be fitch normal")
    #     input$showPassFitch
    #   } 
    # })

    output$plot_out <- renderPlot({ 

      ## Setting up the tree
      if(input$tree == 1) {
        ## Random trees
        if(input$tree_type == "Balanced") {
          ## Balanced tree
          if(log2(input$n_taxa)%%1 == 0)  {
            tree <- ape::stree(input$n_taxa, type = "balanced")
          } else {
            tree <- ape::rtree(input$n_taxa)
          }
        }
        if(input$tree_type == "Left-Right") {
          ## Left-Right tree
          if(input$n_taxa%%2 == 0) {
            left <- ape::stree(input$n_taxa/2+1, type = "left")
            right <- ape::stree(input$n_taxa/2, type = "right")
            tree <- ape::bind.tree(left, right, where = 1)
          } else {
            tree <- ape::rtree(input$n_taxa)
          }
        }
        if(input$tree_type == "Left") {
          ## Left ladder tree
          tree <- ape::stree(input$n_taxa, type = "left")
        }
        if(input$tree_type == "Right") {
          ## Left ladder tree
          tree <- ape::stree(input$n_taxa, type = "right")
        }
        if(input$tree_type == "Random") {
          ## Left ladder tree
          tree <- ape::rtree(input$n_taxa)
        }
      } else {
        ## Tree is a newick format
        tree <- ape::read.tree(text = input$newick_tree)
      }

      ## Setting the character
      if(input$character == 1) {
        ## Generate a random character
        character <- paste(sample(c("0", "1", "2", "-", "?"), ape::Ntip(tree), prob = c(0.535, 0.135, 0.13, 0.1, 0.1), replace = TRUE))
      } else {
        ## User character
        character <- input$character_string
      }

      ## Setting up the number of passes to plot
      # print(input$showPass)

      ## Plot the inapplicable algorithm!
      if(as.numeric(input$method) == 1) {
        plot.inapplicable.algorithm(tree, character, passes = as.vector(as.numeric(input$showPassInapp)))
      }
      if(as.numeric(input$method) == 2) {
        plot.inapplicable.algorithm(tree, character, passes = as.vector(as.numeric(input$showPassFitch)))
      }
      # if(input$action) {
      #   plot.inapplicable.algorithm(tree, character, passes = as.vector(as.numeric(input$showPass)))
      # }

    })#, height = 1200, width = 600)

    # })


    ## DEBUG:
    # output$tree <- renderText({ 
    #   if(input$tree == 1) {
    #     ## The input tree is random
    #     paste("Tree type:", input$tree_type, "with", input$n_taxa, "taxa")
    #   } else {
    #     ## The input tree is user based
    #     paste("Tree:", input$newick_tree, "with", Ntip(ape::read.tree(text = input$newick_tree)), "taxa")
    #   }
    # })

    # output$character <- renderText({ 
    #   if(input$character == 1) {
    #     paste("Character input: Random")
    #   } else {
    #     paste("Character input: ", input$character_string)
    #   }
    # })

    # output$character_string <- renderText({ 
    #   paste("character_string:", input$character_string)
    # })

    # output$tree_type <- renderText({ 
    #   paste("tree_type:", input$tree_type)
    # })

    # output$n_taxa <- renderText({ 
    #   paste("n_taxa:", input$n_taxa)
    # })

    # output$newick_tree <- renderText({ 
    #   paste("newick_tree:", input$newick_tree)
    # })

  }
)