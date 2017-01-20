library(testthat)

context("utilities")

test_that("convert.char works", {
    ## Simple case
    character <- list(1,2,3,-1,c(1,2,3))
    expect_equal(convert.char(character), character)

    ## Numeric case
    char_num <- c(1,2,3,4,5)
    expect_equal(convert.char(char_num), as.list(char_num))

    ## Character case 1
    char_char <- c("1", "2", "3", "-", "?")
    expect_equal(convert.char(char_char), character)

    ## Character case 2
    char_char <- c("123-?")
    expect_equal(convert.char(char_char), character)
})

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

    ## Good functioning with characters
    test <- make.states.matrix(tree, "123-?")
    expect_is(test, "list")
    expect_equal(length(test), 5)
    expect_equal(unique(unlist(lapply(test, length))), 9)
    expect_equal(unlist(test$Char), c(1,2,3,-1,1,2,3))
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


test_that("logicals work", {
    ## Union
    expect_false(get.common(1,2))
    expect_equal(get.common(1,1), 1)
    expect_equal(get.common(c(1),c(1,2,3)), 1)

    ## Intersection.incl
    expect_equal(get.union.incl(1,1), 1)
    expect_equal(get.union.incl(c(1),c(1,2,3)), c(1,2,3))

    ## Intersection.excl
    expect_false(get.union.excl(1,1))
    expect_equal(get.union.excl(c(1),c(1,2,3)), c(2,3))
})
