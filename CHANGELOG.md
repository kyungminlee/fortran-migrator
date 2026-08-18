# Changelog

Release notes live on the [GitHub releases page](https://github.com/kyungminlee/eplinalg/releases);
this file summarizes each tagged version. The current version is in
[VERSION](VERSION).

## v0.22.1

- Release archives link against glibc 2.28 again. v0.22.0 was built on
  glibc 2.39 and, although static archives bake in no glibc version of
  their own, two of the symbols they reference exist only in far newer
  glibc, so the *consumer's* link failed on anything older. Both came
  from the build host's headers rather than from any source change:
  `__isoc23_sscanf` (glibc >= 2.38), because GKlib's `io.c` reaches for
  `_GNU_SOURCE` to declare `getline` and that macro also switches on the
  C23 `scanf` redirects, and `stat64` (glibc >= 2.33), because
  `_FILE_OFFSET_BITS=64` redirects `stat` there. `io.c` now asks for
  `_POSIX_C_SOURCE=200809L`, which declares `getline` and nothing else,
  and large-file support is requested only on 32-bit targets, where it
  is not already the default. No other symbol in the release needed a
  glibc newer than 2.28.
- `scripts/check_glibc_floor.sh` is new, and the release workflow runs it
  on a real glibc 2.28 (a `manylinux_2_28` container) against every
  built archive before publishing, then links a program against the
  whole METIS archive there. A symbol that raises the floor now fails the
  release instead of reaching users.
- MUMPS's out-of-core reopen path calls `open()` with runtime flags and
  no mode. Under the fortify headers of the glibc the release is built
  on, that leaves `mumps_io_basic.c.o` referencing
  `__open_missing_mode` — a diagnostic marker no library anywhere
  defines — so any consumer whose link pulled that object in got an
  undefined reference, on every glibc, not just old ones. The vendored
  source now passes the same `0666` the file's other `open()` already
  does; it is ignored unless the flags really ask to create. Found by
  the new floor check, which is where it was first visible.
- The METIS ordering archives pin `C_STANDARD_REQUIRED` at the target
  rather than inheriting it, so a compiler defaulting to C23 (GCC 15)
  cannot reintroduce the `__isoc23_*` redirects.

## v0.22.0

- Vendored METIS upgraded 5.1.0 -> 5.2.1. From 5.2.1 on METIS no longer
  bundles GKlib — it is a separate upstream repository with a `src/` +
  `include/` split — so the vendored tree keeps that split and the staged
  layout gains a fourth directory, `_mumps_metis_gklib_inc`. 5.2.1
  comments the `IDXTYPEWIDTH` / `REALTYPEWIDTH` self-defaults out of
  `metis.h`; the vendored copy restores them to 32, so translation units
  including the installed header still need no `-D`.
- Vendored MUMPS upgraded 5.9.0 -> 5.9.1.
- Both vendored trees are now reproducible: `scripts/vendor_metis.sh` and
  `scripts/vendor_mumps.sh` fetch pinned upstreams, prune, apply the
  mechanical private rename, and verify their own output. Each was
  previously pristine-upstream-plus-edits recoverable only by
  archaeology, which is how three METIS internals stayed unrenamed for
  the life of the vendored copy. The next upgrade is editing two
  constants.
- The METIS archive no longer exports a single generic symbol. A census
  found 28 globals with no prefix at all (three of them live on the
  `ICNTL(7)=5` ordering path, where a stock system METIS defines the same
  bare names and link order decided which copy won), and a further ~640
  under GKlib's own `gk_` prefix — GKlib's, not ours, and shared with any
  other GKlib in the link. All of them now carry one of four private
  prefixes: `METIS_MUMPS_`, `metis_mumps_`, `libmetis_MUMPS_`, `MUMPS_`.
  The rename is verified on every regeneration, and an unprefixed name a
  future upstream introduces aborts the run and names itself.

## v0.21.0

- Behavior-preserving cleanup across the codebase from an 81-finding
  code-smell audit, following the same pass done for v0.17.0 over the
  four modules and the CMake infrastructure it never opened. No library
  ABI or build-interface changes; staged migrator output is byte-identical
  to v0.20.0 for all five targets, with one intended exception below.
- The `scalapack_tools` recipe is gone, inlined into `pbblas.yaml`. Its
  migrated output stopped being built when NUMROC / ICEIL / ILCM were
  folded into `scalapack_common`, and its only remaining job — naming
  `extern/scalapack-2.2.3/TOOLS` as a symbol context — is now done by the
  one recipe that needed it. Staged trees therefore lose a
  `scalapack_tools/` directory (three files byte-identical to
  `scalapack_common`'s own copies) and two entries from
  `STAGED_LIBRARIES` / `PRIVATIZED_LIBRARIES`. Both staged kind10 trees
  configure identically, so no build consumes any of it.
- Migrator fix: a `REAL(` / `CMPLX(` / `DBLE(` call whose arguments
  contain a *nested* `KIND=` is now migrated instead of skipped. The two
  paren-call rewriters had disagreed — one tracked paren depth, one did
  not — and the depth-blind one would leave `REAL(f(x, KIND=1))`
  unconverted, since the `KIND=` it saw belongs to `f`. No call in the
  current corpus reaches this.
- Target YAMLs gain `param_nuke_extra:`, the names the free-form
  PARAMETER-nuking pass must drop that are not `la_constants_map`
  members. The pass previously carried a hardcoded copy of that map,
  which no second module-based target could have shared.

## v0.20.0

- Installed packages are consumable from single-language projects. The
  generated Config no longer requests `find_package(MPI COMPONENTS C)`
  unconditionally — it asks only for components whose language the
  consuming `project()` enabled, and falls back to reading the MPI flavor
  out of `mpi.h` when no compiled probe can run. A Fortran-only consumer
  could not previously find any MPI-bound package at all.
- Archives that contain Fortran objects now carry the Fortran runtime on
  their exported interface for consumers that have no Fortran driver. A
  C or C++ project resolving the ABI tag with `-DEPLINALG_FORTRAN_TAG`
  links without naming `-lgfortran` itself.
- The pure-C ordering archives (`metis`, `scotch`) export `m` in their
  link interface, and `eplinalgOrdering` is produced by the shared
  package generator rather than a hand-written Config, gaining the
  `find_dependency` hook it previously lacked.
- Both failure paths in the generated Config set a `_NOT_FOUND_MESSAGE`
  naming the consumer's enabled languages and the `-DEPLINALG_MPI_TAG` /
  `-DEPLINALG_FORTRAN_TAG` escape hatches, instead of failing opaquely.
- CI gains an install-and-consume guard: three throwaway projects,
  each enabling exactly one language, configured against an installed
  prefix and linking through imported targets. See
  `test/consume/README.md`.

## v0.19.0

- Documentation converged on the scheme shared with epblas-parallel and
  multifloats: `doc/build.sh` is the canonical build entry, and the
  furo/MyST Sphinx site now covers both the user and developer guides.
- Rendered HTML is published to GitHub Pages
  (https://kyungminlee.github.io/eplinalg/) instead of shipping as a
  release asset: the `eplinalg-vX.Y.Z-docs-html.tar.gz` asset and the
  HTML in the combined archive are gone. Markdown documentation still
  ships in every install tree under `share/doc/eplinalg/`.

## v0.18.0

- User documentation ships with the binary release: the Sphinx HTML
  docs are published as an `eplinalg-vX.Y.Z-docs-html.tar.gz` release
  asset and included in the combined archive under
  `share/doc/eplinalg/html/`. Every per-combo archive carries
  `share/doc/eplinalg/` with the license, changelog, and a pointer
  README.
- CI gains a documentation build check mirroring the release docs job.

## v0.17.0

- Behavior-preserving cleanup across the codebase from a 79-finding
  code-smell audit (migrator internals, recipe YAML, runtime bridges,
  test harness, CMake infrastructure, scripts, examples). Staged
  migrator output is byte-identical to v0.16.0 for all targets;
  no library ABI or build-interface changes.
- New push/PR CI workflow: migrator unit tests plus a stage +
  library-only build of the kind10 / gfortran-15 / Intel MPI slot.

## v0.16.0

- Removed the deprecated `MPI_DD_*` / `MPI_ZZ_*` multifloats MPI op
  aliases (and the matching `mf_mpi_dd_*` / `mf_mpi_zz_*` Fortran
  bindings), as scheduled after the one-release deprecation cycle
  announced in v0.14.0. Use `MPI_MM_*` / `MPI_WW_*`.

## v0.15.0

- Repository restructured to the standard package layout: `extern/`
  (was `external/`), `codegen/{migrator,recipes,targets}/`,
  `src/` (first-party runtime, was `runtime/`),
  `test/{unit,integration}/`, `example/`, root `CMakePresets.json`.
  The staged-tree layout is unchanged, so consumers of the release
  archives are unaffected.
- Version single-sourced from the root `VERSION` file (read by
  `pyproject.toml` and the staged CMake build).
- New root files: `LICENSE` (MIT), `CHANGELOG.md`, `CONTRIBUTING.md`,
  thin root `CMakeLists.txt` with an optional Sphinx `doc` target,
  `.editorconfig`, `.clang-format`, `.clang-tidy`.

## v0.14.0

- Multifloats MPI reduce ops renamed `MPI_MM_*` / `MPI_WW_*`
  (`MPI_DD_*` / `MPI_ZZ_*` kept as deprecated aliases for one release
  cycle).
- Documentation reorganized into `doc/user/` and `doc/dev/`.

## v0.13.x

- v0.13.1: shared-library repackaging requires `-Wl,-z,now`; reference
  repackaging script.
- v0.13.0: position-independent archives — static libraries can be
  repackaged as shared objects.

## v0.12.x

- v0.12.1: extended-precision PBLAS typeset SIGSEGV fix
  (hidden-CHARACTER-length ABI trampolines).
- v0.12.0: MKL-coexistence `ep_` symbol privatization + Intel MPI
  non-commutative reduce-op guard.

## v0.11.0

- Netlib-pristine public C headers via precision-prefixed siblings.

## v0.10.0

- Multifloats acquisition defaults to the upstream binary release.

## v0.9.x

- v0.9.1: `float64x2` diagnostic-print robustness + seq datatype
  defaults.
- v0.9.0: vendored MUMPS upgraded 5.8.2 → 5.9.0; all 10 arithmetics ×
  openmpi/intelmpi/seq validated.

## v0.8.0

- 10-precision linalg stacks, seq consumer path, clean module install,
  ordering-library headers.

## Earlier

v0.1.0 – v0.7.1: incremental build-out of the migrator, the kind16 /
kind10 / multifloats stacks, packaging, and `find_package` support. See
the tag annotations (`git tag -n9`) and the releases page.
