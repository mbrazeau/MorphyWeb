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
        all_states <- unique(all_states[-c(which(is.na(all_states)))]) #, which(all_states == -1))]

        ## Convert inapplicable
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
get.common <- function(a, b) {
options(warn = -1)
    if(!any(a == b)) {
options(warn = 0)
        ## No union
        return(NULL)
    } else {
options(warn = 0)
        ## Union
        if(length(a) >= length(b)) {
            return(sort(a[which(a == b)]))
        } else {
            return(sort(b[which(b == a)]))
        }
    }
}

## Get an intersection inclusive (|)
get.union.incl <- function(a, b) {
    out <- unique(c(a,b))
    if(length(out) == 0) {
        return(NULL)
    } else {
        return(sort(out))
    }
}

## Get an intersection exclusive (^)
get.union.excl <- function(a, b) {
    out <- get.union.incl(a,b)
    out <- out[-which(out == get.common(a,b))]
    if(length(out) == 0) {
        return(NULL)
    } else {
        return(sort(out))
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
        right <- states_matrix$Dp1[desc_anc[1]][[1]] # The node's right descendant
        left <- states_matrix$Dp1[desc_anc[2]][[1]] # The node's left descendant

        ## Get the states in common between the descendants
        common_desc <- get.common(left, right)

        if(!is.null(common_desc)) {
            ## If there is any states in common, set the node to be that one
            states_matrix$Dp1[[node]] <- common_desc

            ## If state in common is actually the inapplicable token, but that both descendants have an applicable, set it to be the union between the applicable states
            if(common_desc == -1 && any(left != -1) && any(right != -1)) {
                states_matrix$Dp1[[node]] <- get.union.incl(left[which(left != -1)], right[which(right != -1)])
            }
        } else {
            ## Else set it to be the union of the descendants
            states_matrix$Dp1[[node]] <- get.union.incl(left, right)

            ## If the node has inapplicable data but that both descendants have also applicable states, remove the inapplicable state from the node
            if(any(states_matrix$Dp1[[node]] == -1) && any(left != -1) && any(right != -1)) {
                states_matrix$Dp1[[node]] <- states_matrix$Dp1[[node]][which(states_matrix$Dp1[[node]] != -1)]
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
        
    ## Pre-condition: if the root is inapplicable AND applicable, remove inapplicable (if there's more than 2 states and one -1)
    if(length(states_matrix$Dp1[[Ntip(tree)+1]]) > 1 && any(states_matrix$Dp1[[Ntip(tree)+1]] == -1)) {
        states_matrix$Up1[[Ntip(tree)+1]] <- states_matrix$Dp1[[Ntip(tree)+1]][-which(states_matrix$Up1[[Ntip(tree)+1]] == -1)]
    } else {
        states_matrix$Up1[[Ntip(tree)+1]] <- states_matrix$Dp1[[Ntip(tree)+1]]
    }

    ## For each node from the root
    for(node in (Ntip(tree)+2:Nnode(tree))) { ## Start past the root (+2)

        curr_node <- states_matrix$Dp1[[node]] # The current node
        ## Select the descendants and ancestors
        desc_anc <- desc.anc(node, tree)
        right <- states_matrix$Dp1[desc_anc[1]][[1]] # The node's right descendant
        left <- states_matrix$Dp1[desc_anc[2]][[1]] # The node's left descendant
        ancestor <- states_matrix$Dp1[desc_anc[3]][[1]]  # The node's ancestor

        if(any(curr_node == -1)) {
            ## If any of the states is inapplicable...
            if(any(curr_node != -1)) {
                ## If any of the states IS applicable
                if(ancestor == -1) {
                    ## If the ancestor state is inapplicable only, set to -1
                    states_matrix$Up1[[node]] <- -1
                } else {
                    ## Else remove the inapplicable
                    states_matrix$Up1[[node]] <- curr_node[which(curr_node != -1)]
                }
            } else {
                ## No state IS applicable
                if(ancestor == -1) {
                    ## If the ancestor state is inapplicable only, set to -1
                    states_matrix$Up1[[node]] <- -1
                } else {
                    ## If the union of left and right has an applicable
                    common_desc <- get.union.incl(left, right)
                    if(any(common_desc != -1)) {
                        ## Set to the union of applicable states
                        states_matrix$Up1[[node]] <- common_desc[which(common_desc != -1)]
                    } else {
                        ## Set to inapplicable
                        states_matrix$Up1[[node]] <- -1
                    }
                }
            }
        } else {
            ## No inapplicable states so don't change any
            states_matrix$Up1[[node]] <- curr_node
        }
    }

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
    
    ## Loop through the nodes
    for(node in rev(Ntip(tree)+1:Nnode(tree))) {

        ## Select the descendants and ancestors
        desc_anc <- desc.anc(node, tree)
        right <- states_matrix$Up1[desc_anc[1]][[1]] # The node's right descendant
        left <- states_matrix$Up1[desc_anc[2]][[1]] # The node's left descendant

        ## If any state on the node is applicable
        if(any(states_matrix$Up1[[node]] != -1)) {
            
            ## Get the states in common between the descendants
            common_desc <- get.common(left, right)

            ## If there is any applicable state in this common, set the node to be that state
            if(any(common_desc != -1)) {
                states_matrix$Dp2[[node]] <- common_desc[which(common_desc != -1)]
            } else {
            ## Else set the node state to be the union of the descendants without the inapplicable tokens
                union_desc <- get.union.incl(left, right)
                states_matrix$Dp2[[node]] <- union_desc[which(union_desc != -1)]
            }
        } else {
            ## Else, leave the state as it was after the first uppass
            states_matrix$Dp2[[node]] <- states_matrix$Up1[[node]]
        }
    }

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

    ## Root state is inherited from the second downpass
    states_matrix$Up2[[Ntip(tree)+1]] <- states_matrix$Dp2[[Ntip(tree)+1]]

    ## For each node from the root
    for(node in (Ntip(tree)+2:Nnode(tree))) { ## Start past the root (+2)

        ## Current node
        curr_node <- states_matrix$Dp2[[node]] # The current node
        ## Select the descendants and ancestors
        desc_anc <- desc.anc(node, tree)
        right <- states_matrix$Dp2[desc_anc[1]][[1]] # The node's right descendant
        left <- states_matrix$Dp2[desc_anc[2]][[1]] # The node's left descendant
        ancestor <- states_matrix$Dp2[desc_anc[3]][[1]] # The node's ancestor

        if(any(curr_node != -1)) { # If any state in the previous pass is not inapplicable
            if(any(ancestor != -1)) { # If any state in the ancestor is not inapplicable
                common_anc_node <- get.common(ancestor, curr_node)
                if(all(common_anc_node == ancestor)) { # If the common state between the ancestor and the final is the ancestor
                    states_matrix$Up2[[node]] <- common_anc_node
                } else { # If the common state between the ancestor and the final is not the ancestor
                    if(get.common(left, right) != FALSE) { # If there is a state in common between left and right
                        states_matrix$Up2[[node]] <- get.union.incl(curr_node, get.common(ancestor, get.union.incl(left, right)))
                    } else { # If there is no state in common between left and right
                        union_desc <- get.union.incl(left, right)
                        if(any(union_desc == -1)) { # If the union of left and right has the inapplicable character
                            if(get.common(ancestor, union_desc) != FALSE) { # If the union of left and right has a state in common with the ancestor
                                states_matrix$Up2[[node]] <- get.union.incl(get.common(ancestor, union_desc), ancestor)
                            } else { # If the union of left and right has no state in common with the ancestor
                                union_all <- get.union.incl(union_desc, ancestor)
                                states_matrix$Up2[[node]] <- union_all[which(union_all != -1)]
                            }
                        } else { # If the union of left and right has no inapplicable character
                            union_node_anc <- get.union.incl(curr_node, ancestor)
                            states_matrix$Up2[[node]] <- union_node_anc
                            if(all(union_node_anc == ancestor)) { # If the state in common between the node and the ancestor is the ancestor
                                states_matrix$Up2[[node]] <- get.common(ancestor, curr_node)
                            }
                        }
                    }
                }
            } else { # If the ancestor has no applicable state
                union_desc <- get.union.incl(left, right)
                if(all(union_desc != FALSE)) { # If there is a state in common between left and right
                    states_matrix$Up2[[node]] <- union_desc
                } else { # If there is no state in common between left and right
                    states_matrix$Up2[[node]] <- curr_node
                }
            }
        } else { # If there is no applicable state in the previous pass
            states_matrix$Up2[[node]] <- curr_node
        }
    }
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

## Internal plot utility: converts characters (-1,0,n,c(-1,0,n)) into character ("-0n?")
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


#' @title Plot inapplicable algorithm
#'
#' @description Plots the results of the inapplicable algorithm
#'
#' @param tree \code{phylo}, a tree
#' @param character \code{character}, a vector of character states
#' @param passes \code{numeric}, the number of passes to plot (default = \code{c(1,2,3,4)})
#' @param show.tip.label \code{logical}, whether to display tip labels (default = \code{FALSE}).
#' @param col.tips.nodes \code{character}, a vector of one or two colors to be used for displaying respectively the tips and the nodes.
##' @param col.character \code{logical} or \code{character}, whether to display the characters states as colors.
#' @param ... any optional arguments to be passed to \code{\link[ape]{plot.phylo}}
#' 
#' @author Thomas Guillerme

## DEBUG
warning("DEBUG")
set.seed(1)
tree <- read.tree(text = "((((((1,2),3),4),5),6),(7,(8,(9,(10,(11,12))))));")
character <- "1100----1100"

tree <- rtree(50)
character <- sample(c(1,0), 50, replace = TRUE)

plot.inapplicable.algorithm(tree, character, n_passes = 1)

plot.inapplicable.algorithm <- function(tree, character, passes = c(1,2,3,4), show.tip.label = FALSE, col.tips.nodes = c("orange", "blue"), ...) { ## Add option of converting states into colors?

    ## SANITIZING
    ## tree character
   

    ## Run the inapplicable algorithm
    states_matrix <- inapplicable.algorithm(tree, character, passes = 4)

    # ## Remove Tree's branch length
    # tree$edge.length <- NULL

    ## Plotting the tree
    plot(tree, show.tip.label = FALSE, type = "phylogram", use.edge.length = FALSE, ...)

    ## Get the node plotting size
    cex = ifelse(Ntip(tree) >= 10, 1/Ntip(tree)*7, 0.7)

    ## Add the tip states
    if(class(character) == "character" && length(character) == 1) {
        tiplabels(as.character(strsplit(as.character(character), "")[[1]]))
    } else {
        tips_labels <- plot.convert.state(states_matrix[[1]][1:Ntip(tree)], missing = TRUE)
        tiplabels(tips_labels)
    }

    ## Get the node labels
    if(n_passes > 1) {
        node_labels <- plot.convert.state(states_matrix[[2]][-c(1:Ntip(tree))])
        node_labels <- paste("1:", node_labels)
        for(pass in 2:n_passes) {
            node_labels <- paste(node_labels, paste(pass, ": ", plot.convert.state(states_matrix[[pass + 1]][-c(1:Ntip(tree))]), sep = ""), sep = "\n")
        }
    } else {
        node_labels <- plot.convert.state(states_matrix[[2]][-c(1:Ntip(tree))])
    }


    ## Plot the node labels
    nodelabels(node_labels, cex = cex)

    ## Add the changes
    # if(show.changes) {

    #     silent <- apply(tree$edge, 1, plot.change)

    #     plot.change <- function(edge, states_matrix, n_passes) {
            
    #         ## Function for collapsing edge states
    #         collapse <- function(edge, states_matrix, n_passes) {
    #             sort(states_matrix[[n_passes + 1]][[edge[1]]])
    #         }


    #         ## Check if the two ends of the edge are different


    #         if(states_matrix[[n_passes + 1]][[edge[1]]] states_matrix[[n_passes + 1]][[edge[2]]])

    #     }
    # }

    return(invisible())
}
