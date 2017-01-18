
## DEBUG
warning("DEBUG")
set.seed(1)
tree <- rtree(5)
character <- c(1,2,2,1,1) ## Inapplicable token will be converted to -1
plot(tree) ; nodelabels(); tiplabels(character, adj = -2)

states_matrix <- make.states.matrix(tree, character)

#' @title Sates matrix
#'
#' @description Creates a states matrix
#'
#' @param tree \code{phylo}, a tree
#' @param character \code{character}, a vector of character states
#'
#' @author Thomas Guillerme

make.states.matrix <- function(tree, character) {

    ## Check if the tree is a tree!
    if(class(tree) != "phylo") {
        stop("The tree must be of class 'phylo'.")
    }

    ## Check if the character is the same length as the tree
    if(Ntip(tree) != length(character)) {
        stop("The tree and character arguments don't match.")
    }

    ## Transform character
    if(class(character) != "numeric") {
        #character <- convert.char(character)
    } else {
        character <- as.list(character)
    }

    ## Set up the list of characters
    filling <- vector("list", Ntip(tree)+Nnode(tree))
    states_matrix <- list("Char" = filling, "Dp1" = filling, "Up1" = filling, "Dp2" = filling, "Up2" = filling)

    ## Add the character into the list
    states_matrix$Char[1:Ntip(tree)] <- character

    return(states_matrix)
}

#' @title Descendants and ancestors
#'
#' @description Returns the right and left descendants and the ancestor of one node
#'
#' @param node \code{numeric}, the number of a node
#' @param tree \code{phylo}, a tree
#'
#' @author Thomas Guillerme

desc.anc <- function(node, tree) {
    descendants <- tree$edge[which(tree$edge[,1] == node),2]
    ancestor <- tree$edge[which(tree$edge[,2] == node),1]
    return(c(descendants, ancestor))
}


#' @title Convert character
#'
#' @description Convert a character if it is not numeric (transforming - into -1 and ? into all characters (but - ))
#'
#' @param character any character vector
#'
#' @author Thomas Guillerme

convert.char <- function(character) {
    ## Convert - into -1
    ## Convert ? into 0:n
    ## Convert 01 into c(01)
    #as.numeric(strsplit(as.character(character), "")[[1]])
}



#' @title First downpass
#'
#' @description Applies a first down pass to a node
#'
#' @param states_matrix \code{matrix}, the matrix containing all the states
#' @param tree \code{phylo}, a tree
#'
#' @author Thomas Guillerme

first.downpass <- function(states_matrix, tree) {

    ## Transferring the characters in the right matrix column
    states_matrix$Dp1 <- states_matrix$Char

    ## Loop through the nodes
    for(node in rev(Ntip(tree)+1:Nnode(tree))) {

        ## Select the descendants and ancestors
        desc_anc <- desc.anc(node, tree)

        ## 
        if(states_matrix$Dp1[desc_anc[1]][[1]] && states_matrix$Dp1[desc_anc[2]][[1]]) {
            states_matrix$Dp1[[node]] <- states_matrix$Dp1[desc_anc[1]][[1]][which(states_matrix$Dp1[desc_anc[1]][[1]] == states_matrix$Dp1[desc_anc[2]][[1]])]
        } else {
            
        }

            

            states_matrix[desc_anc[1], 2] & states_matrix[desc_anc[2], 2]


            )


        cat(node)
    }

    return(states_matrix)
}

#' @title First uppass
#'
#' @description Applies a first uppass pass to a node
#'
#' @param node \code{numeric}, the focal node
#' @param tree \code{phylo}, a tree
#' @param character \code{character}, a vector of character states
#'
#' @author Thomas Guillerme

first.uppass <- function(states_matrix, tree) {

    ## Transferring the characters in the right matrix column
    states_matrix$Up1 <- states_matrix$Char
    
    return(states_matrix)
}

#' @title Second downpass
#'
#' @description Applies a second down pass to a node
#'
#' @param node \code{numeric}, the focal node
#' @param tree \code{phylo}, a tree
#' @param character \code{character}, a vector of character states
#'
#' @author Thomas Guillerme

second.downpass <- function(states_matrix, tree) {

    ## Transferring the characters in the right matrix column
    states_matrix$Dp2 <- states_matrix$Char
    
    return(states_matrix)
}

#' @title Second uppass
#'
#' @description Applies a second up pass to a node
#'
#' @param node \code{numeric}, the focal node
#' @param tree \code{phylo}, a tree
#' @param character \code{character}, a vector of character states
#'
#' @author Thomas Guillerme

second.uppass <- function(states_matrix, tree) {

    ## Transferring the characters in the right matrix column
    states_matrix$Up2 <- states_matrix$Char
    
    return(states_matrix)
}


#' @title Inapplicable algorithm
#'
#' @description Runs a full inapplicable algorithm
#'
#' @param tree \code{phylo}, a tree
#' @param character \code{character}, a vector of character states
#' @param passes \code{numeric}, the number of passes in the tree; from \code{1} to \code{4} (default)
#' 
#' @author Thomas Guillerme

inapplicable.algorithm <- function(tree, character, n_passes = 4) {

    ## Setting up the output state matrix
    states_matrix <- matrix(nrow = Ntip(tree)+Nnode(tree), ncol = 5)
    ## Adding the column names
    colnames(states_matrix) <- c("Char", "1Dp", "1Up", "2Dp", "2Up")
    ## Adding the row names
    if(!is.null(tree$node.label)) {
        nodes_list <- tree$node.label
        rownames(states_matrix) <- c(tree$tip.label, nodes_list)
    } else {
        nodes_list <- (Ntip(tree)+1):(Ntip(tree)+Nnode(tree))
        rownames(states_matrix) <- c(tree$tip.label, nodes_list)
    }

    ## Setting the list of passes
    passes <- list(first.downpass, first.uppass, second.downpass, second.uppass)

    ## Applying the passes for each node
    for (pass in 1:n_passes) {
        sates_matrix <- passes[[pass]](states_matrix, tree)
    }

    return(states_matrix)
}


#' @title Plot inapplicable algorithm
#'
#' @description Plots the inapplicable algorithm
#'
#' @param tree \code{phylo}, a tree
#' @param character \code{character}, a vector of character states
#' @param n_passes \code{numeric}, the number of passes in the tree; from \code{1} to \code{4} (default)
#' 
#' @author Thomas Guillerme

 
## Add which pass to plot (+ how to plot them)

plot.inapplicable.algorithm <- function(tree, character, n_passes = 4, ...) {
        
    ## Get the states matrix
    states_matrix <- inapplicable.algorithm(tree, character, n_passes)

    ## Plot the tree
    plot(tree, ...)




    return(invisible())
}

