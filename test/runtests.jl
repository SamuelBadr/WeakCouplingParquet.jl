using TestItemRunner

@run_package_tests filter = ti -> (contains(ti.filename, "test/"))
