library(testthat)

## convert.char
context("convert.char")
test_that("convert.char works", {
    character <- "01---?010101"
    list_out <- list(0,1,-1,-1,-1,c(0,1),0,1,0,1,0,1)

    ## Not working (return NULL)
    expect_null(convert.char(NULL))
    expect_null(convert.char(matrix(1)))
    expect_null(convert.char(data.frame(1)))

    ## Right output format (returns a list)
    expect_is(convert.char(1), "list")
    expect_is(convert.char("1"), "list")
    expect_is(convert.char(c(1)), "list")
    expect_is(convert.char(list(1)), "list")

    ## Right conversion
    expect_true(all(mapply(function(X,Y) {all(X == Y)}, convert.char(list_out), list_out)))
    expect_true(unlist(convert.char("01-?")[1]) == 0)
    expect_true(unlist(convert.char("01-?")[2]) == 1)
    expect_true(unlist(convert.char("01-?")[3]) == -1)
    expect_true(all(unlist(convert.char("01-?")[4]) == c(-1,0,1)))
})

context("make.states.matrix")
test_that("make.states.matrix works", {
    tree <- rtree(4)
    character <- "01-?"

    ## Not working (error)
    expect_error(make.states.matrix("tree", character, inapplicable = NULL))
    expect_error(make.states.matrix(tree, 1, inapplicable = NULL))

    ## Right output style
    matrix <- make.states.matrix(tree, character)
    expect_is(matrix, "list")
    expect_equal(unique(unlist(lapply(matrix, class))), "list")
    expect_equal(unique(unlist(lapply(matrix, length))), 7)
    expect_equal(length(matrix), 5)
    expect_equal(names(matrix), c("Char", "Dp1", "Up1", "Dp2", "Up2"))

    ## Right output values
    expect_equal(unlist(make.states.matrix(tree, character)[[1]]), c(0,1,-1,-1,0,1))
    expect_equal(unlist(make.states.matrix(tree, character, inapplicable = 1)[[1]]), c(0,1,0,1,0,1))
    expect_equal(unlist(make.states.matrix(tree, character, inapplicable = 2)[[1]]), c(0,1,2,0,1,2))
})

context("desc.anc")
test_that("desc.anc works", {
    tree <- rtree(4)

    ## Not working (wrong input)
    expect_error(desc.anc(1, "tree"))
    expect_error(desc.anc(1, matrix(0)))
    expect_equal(desc.anc("a", tree), numeric(0))
    expect_equal(desc.anc(88, tree), numeric(0))

    ## Right output
    expect_is(desc.anc(6, tree), "integer")

    ## Right answers
    answers <- list(c(5), c(7), c(7), c(6), c(1, 6), c(7, 4, 5), c(2, 3, 6))
    for(test in 1:7) {
        expect_equal(desc.anc(test, tree), answers[[test]])
    }
})

context("get.common")
test_that("get.common works", {

})

context("get.union.incl")
test_that("get.union.incl works", {

})

context("get.union.excl")
test_that("get.union.excl works", {

})

context("fitch.downpass")
test_that("fitch.downpass works", {

})

context("fitch.uppass")
test_that("fitch.uppass works", {

})

context("first.downpass")
test_that("first.downpass works", {

})

context("first.uppass")
test_that("first.uppass works", {

})

context("second.downpass")
test_that("second.downpass works", {

})

context("second.uppass")
test_that("second.uppass works", {

})

context("inapplicable.algorithm")
test_that("inapplicable.algorithm works", {
    ## Test the three cases
})

context("plot.convert.statem")
test_that("plot.convert.statem works", {

})

