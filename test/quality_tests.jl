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
    # Refs: https://github.com/JuliaTesting/Aqua.jl/issues/315
    #       https://github.com/JuliaTesting/Aqua.jl/issues/269
    run_persistent_tasks = get(ENV, "IDRBINDINGDATASET_RUN_PERSISTENT_TASKS", "")
    if isempty(run_persistent_tasks)
        run_persistent_tasks = get(ENV, "GITHUB_ACTIONS", "") == "true" ? "false" : "true"
    end
    run_persistent_tasks = lowercase(strip(run_persistent_tasks)) in ("true", "1", "yes")

    if run_persistent_tasks
        withenv("JULIA_NUM_PRECOMPILE_TASKS" => "1") do
            @test aqua_persistent_tasks_pass(Base.PkgId(IDRBindingDataset))
        end
    else
        @info "Skipping Aqua persistent tasks check in this environment"
        @test_skip "Aqua persistent tasks check skipped"
    end
end

@testitem "Code linting (JET.jl)" begin
    using JET
    using IDRBindingDataset

    JET.test_package(IDRBindingDataset; target_modules = (IDRBindingDataset,))
end
