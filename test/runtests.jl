using TestItemRunner
using TestItems

include("quality_tests.jl")
include("raw_data_download_tests.jl")
include("caid_data_download_tests.jl")
include("sifts_data_download_tests.jl")
include("pdb_release_dates_tests.jl")
include("pdb_release_dataset_tests.jl")
include("sequence_clustering_tests.jl")

@run_package_tests
