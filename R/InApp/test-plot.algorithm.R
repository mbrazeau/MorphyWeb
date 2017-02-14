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
    expect_equal(unlist(convert.char(list_out)), unlist(list_out))
})

context("make.states.matrix")
test_that("make.states.matrix works", {

})

context("desc.anc")
test_that("desc.anc works", {

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

