using TestItems

@testitem "Code quality (Aqua.jl)" begin
    using Aqua
    using IDRBindingDataset

    Aqua.test_all(IDRBindingDataset)
end

@testitem "Code linting (JET.jl)" begin
    using JET
    using IDRBindingDataset

    JET.test_package(IDRBindingDataset; target_modules = (IDRBindingDataset,))
end
