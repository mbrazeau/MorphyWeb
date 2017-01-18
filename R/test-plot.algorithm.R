library(testthat)

context("utilities")

test_that("make.states.matrix works", {
    ## Simple case
    tree <- rtree(5)
    character <- rep(1,5)

    ## Good sanitizing
    expect_error(make.states.matrix(1, character))
    expect_error(make.states.matrix(tree, rep(1,4)))
    expect_error(make.states.matrix(tree, rep(1,6)))

    ## Good functioning
    test <- make.states.matrix(tree, character)
    expect_is(test, "list")
    expect_equal(length(test), 5)
    expect_equal(unique(unlist(lapply(test, length))), 9)
    expect_equal(unlist(test$Char), rep(1, 5))
    expect_equal(names(test), c("Char", "Dp1", "Up1", "Dp2", "Up2"))
})

test_that("desc.anc works", {
    ## Simple case
    set.seed(1)
    tree <- rtree(5)

    ## Good functioning
    expect_equal(desc.anc(6, tree), c(7,8))
    expect_equal(desc.anc(7, tree), c(1,2,6))
    expect_equal(desc.anc(8, tree), c(3,9,6))
    expect_equal(desc.anc(9, tree), c(4,5,8))
})
