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
        options(warn = -1)
        all_states <- as.numeric(character)
        options(warn = 0)
        all_states <- unique(all_states[-c(which(is.na(all_states)), which(all_states == -1))])

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
options(warn = 0)
        ## No union
        return(FALSE)
    } else {
options(warn = 0)
        ## Union
        if(length(a) >= length(b)) {
            return(a[which(a == b)])
        } else {
            return(b[which(b == a)])
        }
    }
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
    states_matrix <- make.states.matrix(tree, character)

    ## Setting the list of passes
    passes <- list(first.downpass, first.uppass, second.downpass, second.uppass)

    ## Applying the passes for each node
    for (pass in 1:n_passes) {
        states_matrix <- passes[[pass]](states_matrix, tree)
    }

    return(states_matrix)
}


#' @title Plot inapplicable algorithm
#'
#' @description Plots the results of the inapplicable algorithm
#'
#' @param tree \code{phylo}, a tree
#' @param character \code{character}, a vector of character states
#' @param show.passes whether to display the ancestral reconstruction for each for passes (\code{TRUE}), any specific pass (\code{1}, \code{2}, \code{3} or \code{4}) or none (\code{FALSE}; default).
#' @param show.changes whether to display the character changes on the tree (\code{TRUE}) or not (\code{FALSE}; default).
#' @param use.pie if \code{show.passes} is not \code{FALSE}, whether to display the node states as a pie or not.
#' @param ... any optional arguments to be passed to \code{\link[ape]{plot.phylo}}
#' 
#' @author Thomas Guillerme

## DEBUG
warning("DEBUG")
set.seed(1)
tree <- rtree(5)
character <- c(1,2,2,1,1)
character <- "1??-?"

plot.inapplicable.algorithm(tree, character)

plot.inapplicable.algorithm <- function(tree, character, show.passes = FALSE, show.changes = FALSE, use.pie = FALSE, ...) { ## Add option of converting states into colors?
    
    plot.convert.state <- function(character, missing = FALSE) {

        plot.convert.inappli <- function(X) {
            return(ifelse(X == -1, "-", X))
        }

        plot.convert.missing <- function(X, all_states) {
            if(length(X) == length(all_states) && all(X == all_states)) {
                return("?")
            } else {
                return(X)
            }
        }

        if(missing) {
            ## Getting all states
            all_states <- unique(unlist(character))[-which(unique(unlist(character)) == -1)]
            ## Convert the missing states
            character <- lapply(character, plot.convert.missing, all_states)
        }

        ## Convert the inapplicables
        character <- lapply(character, plot.convert.inappli)

        ## Convert into character
        return(unlist(lapply(character, function(X) paste(as.character(X), collapse = ""))))
    }

    ## SANITIZING
    ## show.passes
    if(class(show.passes) != "logical") {
        if(class(show.passes) != "numeric") {
            stop("show.passes argument should be logical or any integer(s) between 1 and 4.")
        } else {
            if(any(is.na(match(show.passes, c(1,2,3,4))))) {
                stop("show.passes argument should be logical or any integer(s) between 1 and 4.")
            }
        }
    } else {
        if(show.passes) {
            show.passes <- c(1,2,3,4)
        }
    }
    ## show.changes
    if(class(show.changes) != "logical") {
        stop("show.changes argument must be logical.")
    }
    ## use.pie
    if(class(use.pie) != "logical") {
        stop("use.pie argument must be logical.")
    }

    ## Run the inapplicable algorithm
    warning("DEBUG, 1st dp only") ; n_passes = 1
    states_matrix <- inapplicable.algorithm(tree, character, n_passes)

    ## Remove Tree's branch length
    tree$edge.length <- NULL

    ## Plotting the tree
    plot(tree, show.tip.label = FALSE) ; warning("DEBUG, plot tree option")

    ## Add the tip states
    if(class(character) == "character" && length(character) == 1) {
        tiplabels(as.character(strsplit(as.character(character), "")[[1]]))
    } else {
        tips_labels <- plot.convert.state(states_matrix[[1]][1:Ntip(tree)], missing = TRUE)
        tiplabels(tips_labels)
    }

    ## Add the node labels
    node_labels <- plot.convert.state(states_matrix[[n_passes + 1]][-c(1:Ntip(tree))])
    nodelabels(node_labels)

    ## Add the changes



    return(invisible())
}

    
