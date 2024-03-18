test_that("get kinase PWMs", {
  expect_type(getKinasePWM(), "list")
})

test_that("process phospho-peptides", {
  expect_s3_class(processPhosphopeptides(c("SAGLLS*DEDC")), "data.frame")
})

test_that("score phosphosites", {
  expect_type(scorePhosphosites(getKinasePWM(), 'TGRRHTLAEV'), "double")
})
