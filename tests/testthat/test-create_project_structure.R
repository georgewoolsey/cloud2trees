testthat::test_that("create_project_structure() returns a data.frame", {
  testthat::expect_s3_class(
    object = create_project_structure(
      output_dir = tempdir()
      , input_las_dir = system.file(package = "lidR", "extdata/")
    )
    , class = "data.frame"
  )
})
