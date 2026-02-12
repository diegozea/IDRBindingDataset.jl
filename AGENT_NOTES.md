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

# 11-02-2026

## Goal

Create `scripts/01_download_data.jl` to fetch DisProt raw inputs into `data/raw`, including:
- regions JSON for release `2025_12` with ambiguous entries included and obsolete entries excluded,
- IDPO ontology v0.3.0,
- DisProt/DO/GO/IDPO mapping table v0.2 endpoint.

## Assumptions

- The `release` query parameter should be set from a script constant defaulting to `"2025_12"`.
- The provided URL with `release=2025_12025_12` is a typo; the intended release token is `2025_12`.

## What I tried

- Created `scripts/01_download_data.jl` as an executable orchestration script.
- Implemented constants for release, destination directory, output filenames, and static source URLs.
- Constructed the DisProt search URL from `RELEASE` and fixed query flags:
  - `show_ambiguous=true`
  - `show_obsolete=false`
  - `format=json`
- Added atomic download logic (`Downloads.download` to temporary file then `mv`) to avoid partial outputs.
- Added optional CLI override for release:
  - `julia scripts/01_download_data.jl [release]`
- Executed the script to verify files were downloaded.
- Ran `Pkg.test()` afterward to confirm no regressions in package tests.

## Results and observations

- Script created: `scripts/01_download_data.jl`
- Downloaded files:
  - `data/raw/disprot_regions.json` (`26,821,566` bytes)
  - `data/raw/IDPO_v0.3.0.json` (`9,329` bytes)
  - `data/raw/mapping_ontology_Disprot-DO_GO_IDPO.tsv` (`9,699` bytes)
- `Pkg.test()` passed after the change:
  - `Package | Pass 12 / Total 12`

## Decisions

- Kept download implementation in `scripts/` (not `src/`) because this is pipeline orchestration, not reusable package logic yet.
- Used deterministic output names and a fixed default release constant to support reproducible regeneration.
- Included one positional release override for convenience without changing the default behavior.

## Sources

- DisProt search API endpoint (release-specific regions export):
  - https://disprot.org/api/search?release=2025_12&show_ambiguous=true&show_obsolete=false&format=json
  - Accessed on 11-02-2026.
  - Used for regions JSON download URL and query flags.
- IDPontology JSON (v0.3.0):
  - https://disprot.org/assets/data/IDPO_v0.3.0.json
  - Accessed on 11-02-2026.
  - Used for ontology download URL.
- DisProt/DO/GO/IDPO mappings TSV:
  - https://disprot.org/assets/data/mapping_ontology_Disprot-DO_GO_IDPO.tsv
  - Accessed on 11-02-2026.
  - Used for mapping table download URL.

# 11-02-2026

## Goal

Rename output filenames in `scripts/01_download_data.jl` to shorter, version-agnostic names for readability and generality.

## Assumptions

- Requested filename changes are intentional and should replace the previous output names in future downloads.

## What I tried

- Updated constants in `scripts/01_download_data.jl`:
  - `IDPO_FILENAME` from `IDPO_v0.3.0.json` to `idpo.json`.
  - `IDPO_MAPPINGS_FILENAME` from `mapping_ontology_Disprot-DO_GO_IDPO.tsv` to `mapping_ontology.tsv`.
- Kept `DISPROT_REGIONS_FILENAME = "disprot_regions.json"` unchanged, as requested.
- Re-ran the downloader script to confirm files are written with new names.
- Re-ran `Pkg.test()` to confirm no regressions.

## Results and observations

- Script now writes:
  - `data/raw/disprot_regions.json`
  - `data/raw/idpo.json`
  - `data/raw/mapping_ontology.tsv`
- Previous files with old names remain in `data/raw` from earlier runs:
  - `IDPO_v0.3.0.json`
  - `mapping_ontology_Disprot-DO_GO_IDPO.tsv`
- `Pkg.test()` passed:
  - `Package | Pass 12 / Total 12`

## Decisions

- Applied only constant renames to keep behavior stable and change minimal.
- Did not auto-delete old files to avoid unexpected data loss in `data/raw`.

## Sources

- No new external sources used for this change; implementation is a local refactor of file naming constants.

# 11-02-2026

## Goal

Enhance `scripts/01_download_data.jl` to:
- make IDPO version explicit as a constant,
- clean `data/raw` before each run,
- write `data/raw/README.md` with last download status, versions, and exact download date,
- fail loudly on errors.

## Assumptions

- Cleaning `data/raw` should remove all prior files to prevent cross-run contamination.
- A run report is useful for both success and failure states.

## What I tried

- Added constants:
  - `IDPO_VERSION = "0.3.0"`
  - `IDPO_MAPPINGS_VERSION = "0.2"` (to document versions in report output)
- Parameterized IDPO URL using `IDPO_VERSION`.
- Added strict pre-run cleanup:
  - remove `data/raw` recursively if present,
  - recreate directory before downloads.
- Added UTC timestamp helper and report writer that creates `data/raw/README.md` with:
  - run status (`SUCCESS`/`FAILED`),
  - exact UTC timestamp,
  - DisProt release, IDPO version, IDPO mappings version,
  - file presence/sizes and source URLs.
- Wrapped download loop in `try/catch`:
  - on failure, write `FAILED` report,
  - rethrow exception (no silent failures).
- Added bounded download retries (`3` attempts, fixed backoff) per file to mitigate transient network slowdowns while still failing if all retries fail.
- Ran downloader script and verified output.
- Ran full package tests via `Pkg.test()`.

## Results and observations

- Script now cleans `data/raw` before downloads and repopulates only current outputs.
- `data/raw/README.md` now records a machine-readable run summary for the latest run.
- Successful run produced:
  - `disprot_regions.json`
  - `idpo.json`
  - `mapping_ontology.tsv`
  - `README.md`
- `Pkg.test()` passed:
  - `Package | Pass 12 / Total 12`

## Decisions

- Added failure report generation before rethrowing to preserve diagnostics while maintaining strict error behavior.
- Kept retry count modest to balance robustness and fail-fast behavior.

## Sources

- No new external sources were required for this implementation change.

# 11-02-2026

## Goal

Make `scripts/01_download_data.jl` use `ArgParse.jl` for command-line argument handling.

## Assumptions

- The script should remain runnable from the repo root while keeping package QA checks (Aqua stale deps) clean.
- Existing optional positional release usage should remain supported for backward compatibility.

## What I tried

- Added `ArgParse` to the root package environment first, then refactored script parsing.
- Implemented ArgParse-based parsing with:
  - `--release` option,
  - optional positional release argument.
- Validation showed Aqua failure (`stale_deps`) because `ArgParse` was not used by package `src/`.
- Reworked dependency placement:
  - removed `ArgParse` from root `Project.toml`,
  - added `scripts/Project.toml` with `ArgParse` compat,
  - instantiated `scripts/Manifest.toml`.
- Updated script startup to activate/instantiate script-local environment:
  - `Pkg.activate(@__DIR__)`
  - `Pkg.instantiate()`
- Added conflict guard so users cannot pass both `--release` and positional release simultaneously.

## Results and observations

- `scripts/01_download_data.jl` now uses `ArgParse` for CLI parsing.
- Script is runnable from root (`julia --project=. scripts/01_download_data.jl --release 2025_12`) because it activates its own `scripts/` environment internally.
- CLI conflict handling works:
  - passing both release forms now errors explicitly.
- Full package tests pass after the environment split:
  - `Package | Pass 12 / Total 12`

## Decisions

- Kept `ArgParse` out of root package dependencies to preserve package hygiene and Aqua checks.
- Chose script-local environment as the cleanest way to support script-only dependencies.
- Retained optional positional release support to avoid breaking prior invocation style.

## Sources

- ArgParse.jl README (installation and API conventions):
  - https://raw.githubusercontent.com/carlobaldassi/ArgParse.jl/master/README.md
  - Accessed on 11-02-2026.
- ArgParse.jl examples (authoritative usage of `ArgParseSettings`, `@add_arg_table!`, `parse_args`):
  - https://raw.githubusercontent.com/carlobaldassi/ArgParse.jl/master/examples/argparse_example1.jl
  - https://raw.githubusercontent.com/carlobaldassi/ArgParse.jl/master/examples/argparse_example3.jl
  - Accessed on 11-02-2026.

# 11-02-2026

## Goal

Refactor pipeline logic so reusable functions live in `src/`, split implementation/tests by feature with `include`, add docstrings/comments, and add unit tests for the moved functions.

## Assumptions

- Download/reporting helpers from `scripts/01_download_data.jl` are reusable enough to become module APIs.
- `scripts/01_download_data.jl` should remain a thin orchestration entrypoint.

## What I tried

- Added feature module file `src/raw_data_download.jl` with reusable utilities:
  - `build_disprot_regions_url`
  - `utc_now_iso8601`
  - `clean_directory!`
  - `download_to_file_with_retries!`
  - `write_download_report`
- Added docstrings for all utilities and a focused inline comment around atomic download writes.
- Updated `src/IDRBindingDataset.jl` to export the new API and include the feature file.
- Refactored `scripts/01_download_data.jl` to:
  - include the module from `src/`,
  - keep only argument parsing + orchestration,
  - call module utilities instead of duplicated local functions.
- Split tests by feature:
  - `test/quality_tests.jl` for Aqua/JET checks,
  - `test/raw_data_download_tests.jl` for new utility unit tests,
  - `test/runtests.jl` now includes both files and runs package tests.
- Added deterministic unit tests for URL building, timestamp formatting, directory cleanup, retry behavior, failure behavior, and README report generation.
- Added `Dates` and `Downloads` to package deps and compat because they are now imported by the package module.
- Added `Dates` to `test/Project.toml` because one test imports `Dates` directly.

## Results and observations

- `scripts/01_download_data.jl` still runs successfully and now depends on module utilities.
- Full package tests pass after refactor:
  - `Package | Pass 29 / Total 29`
- Test architecture is now split by functionality, reducing growth pressure in `src/IDRBindingDataset.jl` and `test/runtests.jl`.

## Decisions

- Kept script-specific argument parsing in the script; moved only generally reusable data-download/reporting logic to the module.
- Used one feature file (`raw_data_download.jl`) rather than one file per function, following requested granularity.

## Sources

- No new external sources were required for this refactor; changes are based on repository-local architecture and behavior.

# 11-02-2026

## Goal

Make `scripts/01_download_data.jl` fail with a clear message when downloaded files exceed Git/GitHub file-size limits.

## Assumptions

- The key enforced GitHub constraint for regular Git repositories is the single-file hard limit (100 MiB).
- The requested behavior should trigger after download and before reporting success.

## What I tried

- Verified current GitHub limits from official GitHub Docs.
- Added a reusable validation utility in `src/raw_data_download.jl`:
  - `ensure_files_within_size_limit!(paths; max_bytes, limit_name, remediation)`
- Wired `scripts/01_download_data.jl` to enforce:
  - `GITHUB_SINGLE_FILE_HARD_LIMIT_BYTES = 100 * 1024 * 1024`
  - validation runs after all downloads and before success report writing.
- Ensured failures are explicit and actionable:
  - error lists oversized file paths and exact sizes,
  - error includes configured limit and remediation suggestion (Git LFS/release assets/object storage),
  - catch block continues to write `FAILED` status in `data/raw/README.md` and rethrows.
- Added/updated tests in `test/raw_data_download_tests.jl`:
  - acceptance case for files within limit,
  - failure case asserting clear error text.
- Ran downloader script and full package tests.

## Results and observations

- `scripts/01_download_data.jl` now enforces the GitHub hard file-size limit and fails fast with a clear message when violated.
- Full test suite passed:
  - `Package | Pass 35 / Total 35`

## Decisions

- Enforced the 100 MiB hard limit (blocking threshold) rather than the 50 MiB warning threshold, to align with the “too large for Git/GitHub” failure requirement.

## Sources

- GitHub Docs, "About large files on GitHub":
  - https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github
  - Accessed on 11-02-2026.
  - Used for:
    - 50 MiB Git warning behavior,
    - 100 MiB blocked-file hard limit.
- GitHub Docs, "Repository limits":
  - https://docs.github.com/en/repositories/creating-and-managing-repositories/repository-limits
  - Accessed on 11-02-2026.
  - Used to cross-check:
    - single object hard enforcement at 100 MB.

# 11-02-2026

## Goal

Apply three updates to the raw-data download pipeline:
1. remove `IDPO_MAPPINGS_VERSION` reporting,
2. store downloaded artifacts as `.gz`,
3. report final compressed size in MB in `data/raw/README.md`.

## Assumptions

- Compression should apply to all downloaded artifacts produced by `scripts/01_download_data.jl`.
- The README size summary should reflect compressed files actually saved to disk.

## What I tried

- Removed `IDPO_MAPPINGS_VERSION` from script logic and README report generation path.
- Added gzip download helper in `src/raw_data_download.jl`:
  - `download_to_gzip_file_with_retries!`
- Added reusable byte-size validation/report helpers updates:
  - kept GitHub size-limit validation,
  - updated report writer to include final compressed total size in MB.
- Updated script outputs to gzip names:
  - `disprot_regions.json.gz`
  - `idpo.json.gz`
  - `mapping_ontology.tsv.gz`
- Added `CodecZlib` dependency:
  - root project (`Project.toml` / `Manifest.toml`) because compression is in `src/`,
  - `scripts/Project.toml` / `scripts/Manifest.toml` to satisfy script-local environment.
- Updated tests:
  - report writer expectations (no mappings version line, includes compressed MB summary),
  - new unit test for gzip helper output,
  - existing suite still runs with split test-file layout.
- Ran downloader and full package tests.

## Results and observations

- `data/raw/README.md` now includes:
  - `Final size (compressed, total): <MB>`
  - per-file compressed sizes in bytes and MB.
- `IDPO_MAPPINGS_VERSION` is no longer used by pipeline code or written to `data/raw/README.md`.
- Current compressed outputs observed:
  - `disprot_regions.json.gz` ~ `5.16 MB`
  - `idpo.json.gz` ~ `0.00 MB`
  - `mapping_ontology.tsv.gz` ~ `0.00 MB`
- Full test suite passed:
  - `Package | Pass 36 / Total 36`

## Decisions

- Implemented compression in reusable `src/` code (not script-only) to keep pipeline logic testable and composable.
- Used decimal MB (`bytes / 1_000_000`) for README size reporting, matching the “megabytes” request.

## Sources

- No new external sources were needed beyond previously verified GitHub limits.

# 11-02-2026

## Goal

Update compressed-size reporting to use 4 decimal places in MB.

## What I tried

- Added a reusable formatter in `src/raw_data_download.jl`:
  - `format_mb(bytes) = round(bytes_to_mb(bytes); digits = 4)`
- Applied formatter to:
  - gzip download success log message,
  - total compressed size line in `data/raw/README.md`,
  - per-file compressed size lines in `data/raw/README.md`.
- Re-ran `scripts/01_download_data.jl` to regenerate `data/raw/README.md`.
- Re-ran full test suite.

## Results and observations

- Compressed sizes now print with 4 decimal places, e.g.:
  - `5.1616 MB` total
  - `5.1583 MB` for `disprot_regions.json.gz`
  - `0.0015 MB` and `0.0018 MB` for the smaller files.
- Full test suite passed:
  - `Package | Pass 36 / Total 36`

## Sources

- No new external sources were needed for this formatting-only change.

# 11-02-2026

## Goal

Simplify retry logic in `src/raw_data_download.jl` by using Julia's built-in `retry` with exponential backoff.

## What I tried

- Replaced manual retry loops in:
  - `download_to_file_with_retries!`
  - `download_to_gzip_file_with_retries!`
- New default behavior now uses:
  - `delays = ExponentialBackOff(n = 5)`
- Removed custom retry-specific parameters no longer needed:
  - `retries`
  - `retry_backoff_seconds`
  - `sleep_fn`
- Updated `scripts/01_download_data.jl` call site to pass:
  - `delays = ExponentialBackOff(n = 5)`
- Updated tests to use deterministic custom delay iterators (e.g. `[0.0, 0.0]`) for fast retry-path testing without sleeps.

## Results and observations

- Code path is simpler and uses standard Julia retry utilities.
- Script run succeeded after refactor.
- Full test suite passed:
  - `Package | Pass 35 / Total 35`

## Decisions

- Kept retry logging via `check` callback so failures remain visible during backoff retries.
- Retained explicit per-attempt temp-file cleanup to avoid partial artifacts.

## Sources

- Verified API behavior locally in Julia 1.12.5:
  - `Base.retry`
  - `Base.ExponentialBackOff`

# 11-02-2026

## Goal

Align naming/style in `src/raw_data_download.jl` with Julia conventions requested by the user:
- remove `!` from functions that do not mutate input objects,
- move retry-attempt helper outside the parent function,
- use a private underscore-prefixed helper for the attempt function.

## What I tried

- Refactored retry helper naming and placement:
  - renamed nested `attempt_download!` to top-level `_attempt_download`,
  - moved it outside `download_to_file_with_retries`.
- Renamed public helpers to remove unnecessary bang suffixes:
  - `clean_directory!` -> `clean_directory`
  - `download_to_file_with_retries!` -> `download_to_file_with_retries`
  - `download_to_gzip_file_with_retries!` -> `download_to_gzip_file_with_retries`
  - `ensure_files_within_size_limit!` -> `ensure_files_within_size_limit`
- Updated all call sites/imports:
  - `src/IDRBindingDataset.jl` exports,
  - `scripts/01_download_data.jl`,
  - `test/raw_data_download_tests.jl`.
- Re-ran validation:
  - `julia --project=. -e 'using Pkg; Pkg.test()'`
  - `julia --project=. scripts/01_download_data.jl --release 2025_12`

## Results and observations

- Naming is now consistent with the requested convention and helper visibility preference.
- Retry logic remains equivalent, now with a clearer private helper boundary.
- Tests passed:
  - `Package | Pass 35 / Total 35`
- End-to-end downloader run succeeded (with one transient network retry handled by backoff) and produced expected `.gz` outputs.

## Decisions

- Kept `write_download_report` and other non-mutating names unchanged where they were already convention-aligned.
- Left private numeric format helpers (`bytes_to_mb`, `format_mb`) unexported and internal, without over-applying underscore prefixes.

## Sources

- No new external sources were needed; this was a repository-local naming/style refactor.

# 11-02-2026

## Goal

Switch `scripts/01_download_data.jl` from local `include(...)` loading to package import (`using IDRBindingDataset`) by wiring the `scripts/` environment to the current repo as a local dev dependency, and verify this follows recommended Julia/Pkg practice.

## Assumptions

- The script should keep its own `scripts/Project.toml` environment (for script-only deps like `ArgParse`).
- We want a workflow that works for local iterative package development without duplicating source loading.

## What I tried

- Checked official Julia/Pkg documentation for:
  - non-interactive package operations from scripts,
  - local-path development (`develop`) behavior,
  - when `[sources]` should be used.
- Added the current repo as a dev dependency in `scripts/` using:
  - `julia --project=scripts -e 'using Pkg; Pkg.activate("scripts"); Pkg.develop(path="."); Pkg.resolve(); Pkg.instantiate()'`
- Updated `scripts/01_download_data.jl`:
  - removed `include(joinpath(@__DIR__, "..", "src", "IDRBindingDataset.jl"))`
  - replaced `using .IDRBindingDataset: ...` with `using IDRBindingDataset: ...`
- Removed unused direct script dependency:
  - `CodecZlib` removed from `scripts/Project.toml` (still available transitively via `IDRBindingDataset`).
- Re-ran validations:
  - `julia --project=scripts -e 'using IDRBindingDataset, ArgParse; println("ok")'`
  - `julia --project=. scripts/01_download_data.jl --release 2025_12`
  - `julia --project=. -e 'using Pkg; Pkg.test()'`

## Results and observations

- `scripts/Project.toml` now includes:
  - `IDRBindingDataset` in `[deps]`.
- `scripts/Manifest.toml` now tracks the local package path:
  - `[[deps.IDRBindingDataset]]`
  - `path = ".."`
- `scripts/01_download_data.jl` now imports package APIs directly with `using IDRBindingDataset`.
- Validation passed:
  - script environment imports succeed,
  - downloader script runs successfully,
  - package test suite passed (`35/35`).

## Decisions

- Chose `Pkg.develop(path=...)` as the primary mechanism because it is the official local-development workflow and works cleanly with script environments.
- Did not rely on `include(...)` anymore, since package import is clearer and closer to production use.
- Kept script-local environment activation (`Pkg.activate(@__DIR__)`) so script dependencies remain isolated from the root package environment.

## Sources

- Julia stdlib `Pkg` API docs (functional API recommended for non-interactive/script use):
  - https://docs.julialang.org/en/v1/stdlib/Pkg/
  - Accessed on 11-02-2026.
- Pkg.jl docs, "Managing Packages" (`develop` on local paths and manifest path tracking):
  - https://pkgdocs.julialang.org/v1/managing-packages/
  - Accessed on 11-02-2026.
- Pkg.jl docs, `Project.toml`/`Manifest.toml` guide (`[sources]` section behavior and usage):
  - https://pkgdocs.julialang.org/v1/toml-files/#The-[sources]-section
  - Accessed on 11-02-2026.

# 11-02-2026

## Goal

Simplify downloader CLI to only use flagged options and add explicit configurability for IDPO ontology version.

## What I tried

- Updated `scripts/01_download_data.jl` argument parsing to:
  - keep `--release`,
  - add `--idpo-version`,
  - remove positional release support (`release_positional`).
- Refactored `parse_command_line(args)` return value to include both values:
  - `release`
  - `idpo_version`
- Replaced static `IDPO_URL` constant with runtime URL construction from `idpo_version` so CLI option is effective.
- Added startup log line for the selected IDPO version.
- Kept report metadata aligned by passing parsed `idpo_version` into `write_download_report(...)`.
- Validated CLI behavior with:
  - `julia --project=. scripts/01_download_data.jl --help`
- Ran full test suite:
  - `julia --project=. -e 'using Pkg; Pkg.test()'`

## Results and observations

- Script usage is now:
  - `01_download_data.jl [--release RELEASE] [--idpo-version IDPO-VERSION]`
- Positional release argument is no longer accepted.
- The selected IDPO version now directly controls the downloaded ontology URL.
- Full package tests passed:
  - `Package | Pass 35 / Total 35`

## Decisions

- Kept defaults as constants (`RELEASE = "2025_12"`, `IDPO_VERSION = "0.3.0"`) to preserve deterministic baseline behavior while enabling explicit overrides.

## Sources

- No new external sources were required for this CLI refactor.

# 11-02-2026

## Goal

Improve CLI help text in `scripts/01_download_data.jl` so users see a concrete example of correct execution context (where to run, and environment invocation).

## What I tried

- Updated `ArgParseSettings(description=...)` in `scripts/01_download_data.jl` to include:
  - a repo-root execution example using `julia --project=. scripts/01_download_data.jl ...`,
  - an explicit note that the script activates/instantiates `scripts/Project.toml` internally.
- Validated rendered CLI help with:
  - `julia --project=. scripts/01_download_data.jl --help`
- Ran full package tests:
  - `julia --project=. -e 'using Pkg; Pkg.test()'`

## Results and observations

- `--help` now documents the expected invocation context and a full command example.
- Full package tests passed:
  - `Package | Pass 35 / Total 35`

## Decisions

- Kept the change localized to help text only; no behavioral changes to argument parsing.

## Sources

- No new external sources were required for this documentation/help-text update.

# 11-02-2026

## Goal

Re-check whether `scripts/01_download_data.jl` actually uses the `scripts/` environment and remove ambiguity in CLI guidance.

## What I tried

- Re-ran help in two invocation modes:
  - `julia --project=. scripts/01_download_data.jl --help`
  - `julia scripts/01_download_data.jl --help`
- Confirmed both runs print:
  - `Activating project at .../scripts`
- Updated help example in `ArgParseSettings(description=...)` to avoid confusing `--project=.` and use:
  - `julia scripts/01_download_data.jl --release 2025_12 --idpo-version 0.3.0`
- Re-ran full test suite:
  - `julia --project=. -e 'using Pkg; Pkg.test()'`

## Results and observations

- Script activation is correct: it uses `scripts/Project.toml` because the code executes:
  - `Pkg.activate(@__DIR__)` at startup.
- This remains true regardless of outer invocation style (`--project=.` or not), since startup activation switches to `scripts/`.
- Full package tests passed:
  - `Package | Pass 35 / Total 35`

## Decisions

- Kept internal `Pkg.activate(@__DIR__)` as the source of truth for script environment selection.
- Simplified help text to reduce user confusion around startup project vs script project.

## Sources

- No new external sources were required for this verification and wording update.

# 11-02-2026

## Goal

Rename the orchestration script from `scripts/01_download_data.jl` to `scripts/01_download_disprot_data.jl`.

## What I tried

- Renamed the file:
  - `scripts/01_download_data.jl` -> `scripts/01_download_disprot_data.jl`
- Updated the CLI help description example inside the script to use the new filename.
- Searched for remaining live references to the old filename in `scripts/`, `src/`, `test/`, and `data/`.
- Re-ran package tests:
  - `julia --project=. -e 'using Pkg; Pkg.test()'`

## Results and observations

- New script path is active and `--help` now shows:
  - `usage: 01_download_disprot_data.jl ...`
- No remaining old-name references were found in active code/docs directories (`scripts`, `src`, `test`, `data`).
- Full test suite passed:
  - `Package | Pass 35 / Total 35`

## Decisions

- Historical references in prior `AGENT_NOTES.md` entries were left unchanged (append-only notebook policy).

## Sources

- No external sources were required for this refactor.

# 11-02-2026

## Goal

Require Julia 1.12 across the repository, update CI to match, and adopt Julia-version-specific manifest files.

## Assumptions

- The package baseline should become Julia `1.12` (drop `1.10` support).
- CI should test the supported baseline explicitly rather than a mixed-version matrix.
- Version-specific manifests should be used for root and sub-environments (`scripts/`, `test/`, `docs/`).

## What I tried

- Reviewed official recommendations for:
  - Julia version compatibility in `Project.toml`,
  - version-specific manifest naming/behavior,
  - CI version pinning practices with `setup-julia`.
- Updated compatibility bounds:
  - `Project.toml`: `julia = "1.12"`
  - `scripts/Project.toml`: `julia = "1.12"`
  - `test/Project.toml`: added `[compat] julia = "1.12"`
  - `docs/Project.toml`: added `[compat] julia = "1.12"`
- Updated CI workflow:
  - `.github/workflows/CI.yml` test matrix now runs only `1.12`.
  - `.github/workflows/CI.yml` docs job now uses `version: '1.12'`.
- Adopted Julia-version-specific manifests by renaming:
  - `Manifest.toml` -> `Manifest-v1.12.toml`
  - `scripts/Manifest.toml` -> `scripts/Manifest-v1.12.toml`
  - `test/Manifest.toml` -> `test/Manifest-v1.12.toml`
  - `docs/Manifest.toml` -> `docs/Manifest-v1.12.toml`
- Verified environments and behavior:
  - `julia --project=. -e 'using Pkg; Pkg.test()'`
  - `julia scripts/01_download_disprot_data.jl --help`
  - `julia --project=scripts -e 'using Pkg; Pkg.instantiate(); using ArgParse, IDRBindingDataset'`
  - `julia --project=test -e 'using Pkg; Pkg.instantiate()'`
  - `julia --project=docs -e 'using Pkg; Pkg.instantiate(); using IDRBindingDataset'`

## Results and observations

- The repository now enforces Julia 1.12 at package and sub-environment level.
- CI is aligned to Julia 1.12 for both tests and docs.
- All manifests are now version-specific for Julia 1.12.
- Validation succeeded:
  - package tests passed (`35 / 35`),
  - script help path works,
  - scripts/test/docs environments instantiate and load correctly.

## Decisions

- Pinned CI to `1.12` only to match the stated support baseline and ensure deterministic manifest usage (`Manifest-v1.12.toml`).
- Added explicit julia compat entries in `test/` and `docs/` projects for clarity and consistency with root/scripts policy.

## Sources

- Julia Pkg docs, `Project.toml` / `Manifest.toml`:
  - https://pkgdocs.julialang.org/v1/toml-files/
  - Accessed on 11-02-2026.
  - Used for:
    - `julia` compatibility via `[compat]`,
    - version-specific manifest naming (`Manifest-v{major}.{minor}.toml`) and precedence behavior.
- `julia-actions/setup-julia` README:
  - https://github.com/julia-actions/setup-julia
  - Accessed on 11-02-2026.
  - Used for CI version configuration recommendations (quoted version strings and explicit pinning).

# 11-02-2026

## Goal

Add Julia formatting guardrails to repository automation and instructions:
- add a `Format` workflow patterned after MIToS but updated to current recommendations,
- format all Julia source files in the repository,
- add explicit formatter instructions in `AGENTS.md`.

## What I tried

- Retrieved MIToS workflow reference:
  - `https://github.com/diegozea/MIToS.jl/blob/master/.github/workflows/Format.yml`
  - confirmed it uses `julia-actions/julia-format` with suggestion comments on PRs.
- Verified current upstream formatter action recommendations and release:
  - `julia-actions/julia-format` README usage for `@v4`,
  - latest release `v4.1.0`.
- Added `.github/workflows/Format.yml` using MIToS-style PR suggestion workflow, updated to:
  - `uses: julia-actions/julia-format@v4`
  - `version: '2'`
  - `suggestion-label: 'format-suggest'`
  - explicit permissions for PR suggestions (`pull-requests: write`).
- Ran repository-wide formatting with JuliaFormatter:
  - `julia -e 'import Pkg; Pkg.activate(temp=true); Pkg.add(name="JuliaFormatter", version="2"); using JuliaFormatter; format(".")'`
- Updated `AGENTS.md` to require running JuliaFormatter after code changes (with concrete command).
- Ran full test suite:
  - `julia --project=. -e 'using Pkg; Pkg.test()'`

## Results and observations

- New workflow added:
  - `.github/workflows/Format.yml`
- Julia files were reformatted (style-only diffs) in:
  - `docs/make.jl`
  - `scripts/01_download_disprot_data.jl`
  - `src/raw_data_download.jl`
  - `test/raw_data_download_tests.jl`
- Full package tests passed:
  - `Package | Pass 35 / Total 35`

## Decisions

- Kept the workflow behavior close to MIToS (PR suggestions, label-gated option) for consistency.
- Upgraded from MIToS's `julia-format@v3` to `@v4` based on current upstream recommendation.
- Used JuliaFormatter major version `2` in workflow and local formatting command for consistency.

## Sources

- MIToS Format workflow reference:
  - https://github.com/diegozea/MIToS.jl/blob/master/.github/workflows/Format.yml
  - Accessed on 11-02-2026.
- `julia-actions/julia-format` repository and README (usage of `@v4`):
  - https://github.com/julia-actions/julia-format
  - Accessed on 11-02-2026.
- `julia-actions/julia-format` releases (latest `v4.1.0`):
  - https://github.com/julia-actions/julia-format/releases
  - Accessed on 11-02-2026.

# 11-02-2026

## Goal

Add explicit guidance in `AGENTS.md` to prefer `julia-mcp` when available, and validate it in this environment for routine Julia tasks (tests/formatting context).

## What I tried

- Updated `AGENTS.md` operational conventions to state:
  - prefer `julia-mcp` for Julia execution when available (especially tests and formatting),
  - use shell `julia` only as fallback,
  - minimum validation should run full tests (prefer via `julia-mcp`).
- Verified `julia-mcp` availability and execution by running full package tests through MCP:
  - `Pkg.test()` in env `IDRBindingDataset.jl` via persistent MCP Julia session.

## Results and observations

- `julia-mcp` is available and functional in this environment.
- Full test suite run via MCP passed:
  - `Package | Pass 35 / Total 35`

## Decisions

- Kept formatter command examples in AGENTS generic (policy-level) while making MCP preference explicit.
- Used test execution as MCP proof-of-use for this change, which directly exercises the requested low-TTFX path.

## Sources

- No external sources were required for this repository-internal policy update and runtime verification.

# 11-02-2026

## Goal

Improve `AGENTS.md` Julia MCP guidance by adding concrete usage examples instead of only policy text.

## What I tried

- Updated the `Repo operational conventions (agent)` section in `AGENTS.md` to include explicit Julia MCP examples for:
  - running tests in the repo environment,
  - running formatting in the repo environment.
- Kept shell-Julia fallback guidance intact.
- Validated Julia MCP operation in this environment by running full tests via MCP:
  - `Pkg.test()` with `env_path = /home/diego.zea/.julia/dev/IDRBindingDataset.jl`.

## Results and observations

- `AGENTS.md` now includes actionable Julia MCP examples (`env_path` + `code` pairs).
- MCP validation passed:
  - `Package | Pass 35 / Total 35`

## Decisions

- Used minimal example payloads so the pattern is clear and reusable without over-specifying tool internals.

## Sources

- No external sources were required for this documentation/policy update.

# 12-02-2026

## Goal

Update the `scripts/01_download_disprot_data.jl` cleanup step so it removes only the files that script downloads (plus its generated report), instead of deleting the entire `data/raw` directory.

## Assumptions

- The intended cleanup scope is limited to files produced by script 01 in `data/raw`.
- Keeping unrelated files in `data/raw` is desired to avoid accidental data loss.

## What I tried

- Added a new helper `remove_files(paths)` in `src/raw_data_download.jl` that deletes only listed files if they exist.
- Exported `remove_files` from `src/IDRBindingDataset.jl`.
- Updated `scripts/01_download_disprot_data.jl` to:
  - import `remove_files` instead of `clean_directory`,
  - build the list of expected downloaded file paths from `sources`,
  - remove only those files plus `data/raw/README.md` before download.
- Added a unit test in `test/raw_data_download_tests.jl` to verify `remove_files` does not delete unlisted files.

## Results and observations

- The script now performs targeted cleanup and preserves unrelated files in `data/raw`.
- Added test coverage for selective file deletion behavior.
- Could not run Julia formatter/tests in this environment because `julia` is not installed on PATH.

## Decisions

- Kept existing `clean_directory` helper unchanged for compatibility and existing tests.
- Implemented targeted cleanup as a separate reusable helper to keep script logic concise and testable.

## Sources

- No external sources used; change is based on repository-local requirements and code.
