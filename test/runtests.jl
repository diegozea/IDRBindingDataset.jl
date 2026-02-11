using TestItemRunner
using TestItems

include("quality_tests.jl")
include("raw_data_download_tests.jl")

@run_package_tests
