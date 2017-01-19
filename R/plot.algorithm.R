
## DEBUG
warning("DEBUG")
set.seed(1)
tree <- rtree(5)
character <- c(1,2,2,1,1) ## Inapplicable token will be converted to -1
plot(tree) ; nodelabels(); tiplabels(character, adj = -2)

states_matrix <- make.states.matrix(tree, character)
states_matrix <- first.downpass(states_matrix, tree)

#' @title Convert character
#'
#' @description Convert a character if it is not numeric (transforming - into -1 and ? into all characters (but - ))
#'
#' @param character any character vector
#'
#' @author Thomas Guillerme

convert.char <- function(character) {

    convert.inappli <- function(X) {
        return(ifelse(X == "-", -1, X))
    }

    convert.missing <- function(X, all_states) {
        if(X == "?") {
            return(all_states)
        } else {
            return(X)
        }
    }

    ## Character is a list
    if(class(character) == "list") {
        if(unique(unlist(lapply(character, class))) != "numeric") {
            stop("Character list does not contain only numeric values.")
        } else {
            return(character)
        }
    }

    ## Character is a vector (numeric)
    if(class(character) == "numeric") {
        return(as.list(character))
    }
    
    ## Character is not numeric
    if(class(character) == "character") {
        if(length(character) == 1) {
            #Split the character chain
            character <- as.character(strsplit(as.character(character), "")[[1]])
        } 
        ## Convert into list
        character <- as.list(character)

        ## Get all states
        all_states <- as.numeric(character)
        all_states <- all_states[-c(which(is.na(all_states)), which(all_states == -1))]

        ## Convert inapplicables
        character <- lapply(character, convert.inappli)
        
        ## Convert missing
        character <- lapply(character, convert.missing, all_states)

        ## Convert into numeric
        return(lapply(character, as.numeric))
    }
}

#' @title Sates matrix
#'
#' @description Creates a states matrix
#'
#' @param tree \code{phylo}, a tree
#' @param character Either vector of character states (\code{"numeric"} or \code{"character"}) or a list of the same length of than the tips in the tree (see details)
#'
#' @details
#' If \code{character} argument is a list, each element of the list must be a \code{"numeric"} vector with \code{"?"} being all states and \code{"-"} being \code{-1}.
#' @author Thomas Guillerme

make.states.matrix <- function(tree, character) {

    ## Check if the tree is a tree!
    if(class(tree) != "phylo") {
        stop("The tree must be of class 'phylo'.")
    }

    ## Transform character
    if(class(character) != "list") {
        character <- convert.char(character)
    }

    ## Check if the character is the same length as the tree
    if(Ntip(tree) != length(character)) {
        stop("The tree and character arguments don't match.")
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

## Get an union (&)
get.union <- function(a, b) {
    options(warn = -1)
    if(!any(a == b)) {
        ## No union
        return(FALSE)
    } else {
        ## Union
        if(length(a) >= length(b)) {
            return(a[which(a == b)])
        } else {
            return(b[which(b == a)])
        }
    }
    options(warn = 0)
}

## Get an intersection inclusive (|)
get.intersection.incl <- function(a, b) {
    out <- unique(c(a,b))
    if(length(out) == 0) {
        return(FALSE)
    } else {
        return(out)
    }
}

## Get an intersection exclusive (^)
get.intersection.excl <- function(a, b) {
    out <- get.intersection.incl(a,b)
    out <- out[-which(out == get.union(a,b))]
    if(length(out) == 0) {
        return(FALSE)
    } else {
        return(out)
    }
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
        right <- states_matrix$Dp1[desc_anc[1]][[1]]
        left <- states_matrix$Dp1[desc_anc[2]][[1]]

        ## Get the union between the descendants
        union_desc <- get.union(left, right)

        if(union_desc != FALSE) {
            ## If there is a union, set the node to be the union
            states_matrix$Dp1[[node]] <- union_desc

            ## If the union is actually the inapplicable token, and both descendants have an applicable, set it to be the intersection between applicable
            if(union_desc == -1 && any(left != -1) && any(right != -1)) {
                states_matrix$Dp1[[node]] <- get.intersection.incl(left[-which(left == -1)], right[-which(right == -1)])
            }
        } else {
            ## Else set it to be the intersection
            states_matrix$Dp1[[node]] <- get.intersection.incl(left, right)

            ## If the node has inapplicable data
            if(any(states_matrix$Dp1[[node]] == -1)) {
                ## If both left and right have applicable, set the node to be applicable only
                if(any(left != -1) && any(right != -1)) {
                    states_matrix$Dp1[[node]] <- states_matrix$Dp1[[node]][-which(states_matrix$Dp1[[node]] == -1)]
                }
            }
        }
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

    