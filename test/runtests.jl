using TestItemRunner

@run_package_tests filter = ti -> (contains(ti.filename, "test/"))
# @run_package_tests filter = ti -> (contains(ti.filename, "benchmarks"))
