using IDRBindingDataset
using Test
using Aqua
using JET

@testset "IDRBindingDataset.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(IDRBindingDataset)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(IDRBindingDataset; target_defined_modules = true)
    end
    # Write your tests here.
end
