orthomclx vendored MCL note
===========================

This repository currently includes upstream MCL version `14-137` in `src/mcl/`.

Users are welcome to replace it with a newer upstream release from:

- https://micans.org/mcl/

When updating the vendored copy:

1. Preserve the upstream `COPYING` and other license files.
2. Keep the upstream source layout intact as much as possible.
3. Rebuild the bundled executable with `make build-mcl`.
