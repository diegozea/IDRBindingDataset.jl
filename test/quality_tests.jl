using TestItems

@testitem "Code quality (Aqua.jl)" begin
    using Test
    using Aqua
    using IDRBindingDataset

    function aqua_persistent_tasks_pass(pkgid::Base.PkgId; attempts::Int = 3)
        for attempt in 1:attempts
            if !Aqua.has_persistent_tasks(pkgid)
                return true
            end

            if attempt < attempts
                @info "Retrying Aqua persistent tasks check after a failed attempt" attempt attempts
            end
        end

        return false
    end

    Aqua.test_all(IDRBindingDataset; persistent_tasks = false)

    # Aqua can report a false positive in this check when wrapper precompilation
    # exits early for reasons unrelated to persistent tasks.
    # Ref: https://github.com/JuliaTesting/Aqua.jl/issues/315
    withenv("JULIA_NUM_PRECOMPILE_TASKS" => "1") do
        @test aqua_persistent_tasks_pass(Base.PkgId(IDRBindingDataset))
    end
end

@testitem "Code linting (JET.jl)" begin
    using JET
    using IDRBindingDataset

    JET.test_package(IDRBindingDataset; target_modules = (IDRBindingDataset,))
end
