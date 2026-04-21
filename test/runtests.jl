using TestItemRunner
using TestItems

include("quality_tests.jl")
include("shared_helpers_tests.jl")
include("raw_data_download_tests.jl")
include("caid_data_download_tests.jl")
include("caid3_invalid_regions_tests.jl")
include("go_data_download_tests.jl")
include("disprot_challenge_terms_tests.jl")
include("sifts_data_download_tests.jl")
include("pdb_release_dates_tests.jl")
include("pdb_release_dataset_tests.jl")
include("sequence_clustering_tests.jl")
include("fixed_labels_tests.jl")
include("elm_data_download_tests.jl")
include("elm_motifs_tests.jl")
include("elm_motif_clustering_tests.jl")
include("dataset_stats_tests.jl")
include("dataset_finalization_tests.jl")
include("dataset_json_export_tests.jl")
include("caid3_binding_labels_tests.jl")

@run_package_tests
