test_that("get kinase PWMs", {
  expect_type(getKinasePWM(), "list")
})

test_that("process phosphosites", {
  expect_s3_class(processPhosphosites(c("SAGLLS*DEDC")), "data.frame")
})

test_that("score phosphosites", {
  expect_type(scorePhosphosites(getKinasePWM(), 'TGRRHTLAEV'), "double")
})
