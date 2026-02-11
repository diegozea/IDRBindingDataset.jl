# 11-02-2026

## Goal

Migrate package tests to the Test Item framework so individual tests can run in VS Code while preserving `Pkg.test()` behavior.

## Assumptions

- The requested migration targets the official Test Item ecosystem (`TestItems.jl` + `TestItemRunner.jl`).
- Existing quality/lint checks (Aqua + JET) should stay as test coverage.

## What I tried

- Inspected current test setup in `test/runtests.jl` and `test/Project.toml`.
- Verified current usage from official Test Item documentation and upstream package metadata.
- Rewrote `test/runtests.jl` to:
  - define tests with `@testitem`,
  - run with `@run_package_tests`.
- Added `TestItems` and `TestItemRunner` to `test/Project.toml`.
- Resolved and instantiated test dependencies to update `test/Manifest.toml`.
- Ran full package tests with `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Updated deprecated JET keyword from `target_defined_modules` to `target_modules=(IDRBindingDataset,)`, then re-ran full tests.

## Results and observations

- Final test run passed:
  - `Package | Pass 12 / Total 12`
- Test-item runner works with package tests through `@run_package_tests`.
- Added test dependencies resolved to:
  - `TestItemRunner v1.1.4`
  - `TestItems v1.0.0`
- Environment details observed in test manifest:
  - Julia `1.12.5`
  - Aqua `0.8.14`
  - JET `0.11.3`

## Decisions

- Kept the same two logical checks (Aqua and JET), only changed execution model to test items.
- Applied the JET deprecation fix now to keep CI/local runs warning-clean and future-compatible.
- Did not add `@testitem` tags yet because current suite has only two top-level checks; tags can be added once functional tests are split into multiple files.

## Sources

- Julia VS Code docs, "Test Item Framework":
  - https://www.julia-vscode.org/docs/stable/userguide/testitems/
  - Accessed on 11-02-2026.
  - Used for:
    - `@testitem` structure,
    - package test integration with `@run_package_tests`,
    - default imports behavior in package tests.
- TestItemRunner.jl upstream metadata:
  - https://raw.githubusercontent.com/julia-vscode/TestItemRunner.jl/main/Project.toml
  - Accessed on 11-02-2026.
  - Used to verify package UUID (`f8b46487-2199-4994-9208-9a1283c18c0a`) and compatibility relationship with TestItems.
- TestItems.jl upstream metadata:
  - https://raw.githubusercontent.com/julia-vscode/TestItems.jl/main/Project.toml
  - Accessed on 11-02-2026.
  - Used to verify package UUID (`1c621080-faea-4a02-84b6-bbd5e436b8fe`).
