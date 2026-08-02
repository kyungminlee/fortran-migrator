# Installed-package consumer guards

Three throwaway projects that are configured **against an installed prefix**,
not as part of the staged tree. They exist to catch defects that are invisible
to the in-tree build, because in-tree everything is compiled by one `project()`
that enables Fortran, C and C++ at once and links through the gfortran driver.

Each directory deliberately enables **one** language:

| Directory        | `project(... LANGUAGES ...)` | Guards |
|------------------|------------------------------|--------|
| `fortran/`       | `Fortran`                    | consumer-side MPI component detection, and the `MPI::MPI_C` forwarding shim |
| `c/`             | `C`                          | libm in the exported interface of the pure-C ordering archives, and that `eplinalgOrdering` is a normally-generated package |
| `c-fortran-tag/` | `C`                          | the `-DEPLINALG_FORTRAN_TAG` escape hatch, and that the package carries the producer's Fortran runtime for a consumer with no Fortran driver |

**The single-language split is the whole point.** Every defect these guard is
about a language being *absent* from the consumer's `project()`, so one
consumer that enables both Fortran and C would catch neither: with C enabled
`find_package(MPI COMPONENTS C)` succeeds, and with Fortran enabled the link
pulls in `libgfortran`, whose `NEEDED libm.so.6` silently satisfies METIS's
`sqrt`/`log`/`pow` — and which supplies `_gfortran_*` whether or not the
package declares it.

`c/` and `c-fortran-tag/` are both C-only but guard different things. `c/`
consumes a package whose archives are pure C, so a C project resolves its tag
unaided. `c-fortran-tag/` consumes a package whose archives carry Fortran
objects, so their tag depends on a Fortran compiler this project does not
have: the consumer names the tag itself, and the package must then supply
`libgfortran` et al. Two things must never be added to `c-fortran-tag/`, both
of which would make it pass while testing nothing — **Fortran in `LANGUAGES`**
(restores the driver, which supplies the runtime on its own) and **`-lgfortran`
on the link line** (exactly the manual step this path exists to remove).

They also link **imported targets** (`eplinalg::scalapack`, `eplinalg::metis`,
`eplinalg::blas`), never `$<TARGET_FILE:...>` paths. Naming an archive by path
bypasses the package Config, the tag resolution and the exported link interface
— i.e. every mechanism under test. `example/mmsolve` links by path and
hand-appends `m quadmath stdc++`, which is why it cannot serve as this guard.

Driven by the consumer steps in `.github/workflows/ci.yml`. To run by hand
against a prefix you have already installed into:

```sh
cmake -S test/consume/fortran -B /tmp/cf -DCMAKE_PREFIX_PATH=<prefix> \
      -DCMAKE_Fortran_COMPILER=gfortran-15
cmake --build /tmp/cf && /tmp/cf/consume_fortran

cmake -S test/consume/c -B /tmp/cc -DCMAKE_PREFIX_PATH=<prefix> \
      -DCMAKE_C_COMPILER=gcc-15
cmake --build /tmp/cc && /tmp/cc/consume_c

cmake -S test/consume/c-fortran-tag -B /tmp/ct -DCMAKE_PREFIX_PATH=<prefix> \
      -DCMAKE_C_COMPILER=gcc-15 -DEPLINALG_FORTRAN_TAG=gfortran-15
cmake --build /tmp/ct && /tmp/ct/consume_c_fortran_tag
```
