test_that("get kinase PWMs", {
  expect_type(get_kinase_pwms(), "list")
})

test_that("process phosphosites", {
  expect_s3_class(process_phosphosites(c("SAGLLS*DEDC")), "data.frame")
})

test_that("score phosphosites", {
  expect_type(score_phosphosites(get_kinase_pwms(), 'TGRRHTLAEV'), "double")
})
