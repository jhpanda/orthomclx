# orthomcl_c

Modern C/Python rewrite of the OrthoMCL workflow.

## Rationale

This project replaces the old Perl + MySQL OrthoMCL implementation with a
file-based Python/C pipeline that is easier to run on HPC systems.

The main goals are:

- keep the OrthoMCL algorithm and output structure
- remove the MySQL dependency and hidden database state
- replace legacy Perl scripts with maintainable Python code
- move heavy pair-building work into C for better runtime and memory behavior
- keep BLAST-like tabular input and MCL clustering support

This rewrite is meant to preserve the useful biology while making the workflow
simpler to install, easier to inspect, and more scalable for large datasets.

## What It Does

Current implemented workflow:

1. adjust FASTA headers into OrthoMCL-compatible IDs
2. filter compliant FASTA files
3. parse BLAST m8 into OrthoMCL-style similarity rows
4. compile parsed similarities into compact binary files
5. compute:
   - orthologs
   - inparalogs
   - co-orthologs
   - strict `1:1` reciprocal best hits (`rbh_1to1`)
6. generate `mclInput`
7. convert MCL output into `groups.txt`

External tools still expected:

- `blastp` or another search tool that produces BLAST-like tabular output
- `mcl` for clustering, either installed on `PATH` or vendored under `src/mcl/`

This repository currently vendors MCL version `14-137` under `src/mcl/`.
Users are also welcome to download and use the latest upstream MCL release from
[https://micans.org/mcl/](https://micans.org/mcl/).

## Why Python And C

Python is used for:

- CLI and workflow orchestration
- FASTA and BLAST parsing
- output handling
- tests and validation helpers

C is used for:

- indexed numeric computation on large similarity tables
- ortholog / inparalog / co-ortholog generation
- faster and more memory-efficient large-dataset processing

## Repository Layout

```text
src/
  orthomcl/   Python package
  c/          C source files
  mcl/        Vendored upstream MCL source (optional)
scripts/      Useful wrapper scripts
tests/        Automated tests and validation datasets
```

## Main CLI

The package exposes the `orthomclx` command:

```bash
orthomclx adjust-fasta
orthomclx filter-fasta
orthomclx parse-blast
orthomclx compile-similarities
orthomclx pairs
orthomclx indexed-orthologs
orthomclx indexed-inparalogs
orthomclx indexed-coorthologs
orthomclx indexed-pairs
orthomclx mcl
orthomclx mcl-to-groups
orthomclx run
```

## Typical Usage

### 1. Integrated one-shot run

For the common case, you can run the whole file-based indexed workflow from the
top level:

```bash
orthomclx --input compliantFasta \
  --blast blast.tsv \
  --out run_dir \
  --pcut 30 \
  --ecut -3 \
  --jobs 2 \
  --run-mcl
```

This parses BLAST output, compiles similarities, runs the indexed pair stages,
optionally runs `mcl`, and writes:

- `run_dir/similarSequences.txt`
- `run_dir/compiled/`
- `run_dir/pairs/orthologs.txt`
- `run_dir/pairs/inparalogs.txt`
- `run_dir/pairs/coorthologs.txt`
- `run_dir/pairs/rbh_1to1.txt`
- `run_dir/mclInput`
- `run_dir/mclOutput` when `--run-mcl` is used
- `run_dir/groups.txt` when `--run-mcl` or `--mcl-output` is used

### 1. Parse BLAST output

```bash
PYTHONPATH=src python3 -m orthomcl.cli parse-blast blast.tsv compliantFasta -o similarSequences.txt
```

### 2. Compile similarities into binary form

```bash
PYTHONPATH=src python3 -m orthomcl.cli compile-similarities similarSequences.txt compiled_dir
```

This produces:

- `similarities.bin`
- `proteins.tsv`
- `taxa.tsv`

### 3. Run the indexed pipeline

```bash
PYTHONPATH=src python3 -m orthomcl.cli indexed-pairs compiled_dir out_dir \
  --percent-match-cutoff 30 \
  --evalue-exp-cutoff -3 \
  --jobs 2
```

This writes:

- `out_dir/pairs/orthologs.txt`
- `out_dir/pairs/inparalogs.txt`
- `out_dir/pairs/coorthologs.txt`
- `out_dir/pairs/rbh_1to1.txt`
- `out_dir/mclInput`

### 4. Convert MCL output to groups

```bash
PYTHONPATH=src python3 -m orthomcl.cli mcl-to-groups mclOutput OG 1000 -o groups.txt
```

### 5. Run MCL directly

```bash
PYTHONPATH=src python3 -m orthomcl.cli mcl mclInput -o mclOutput --inflation 1.5 --threads 4
```

If you vendor upstream MCL into `src/mcl/`, you can build a bundled binary with:

```bash
make build-mcl
```

After that, `orthomclx mcl ...` and `--run-mcl` will automatically prefer
`build/mcl` over a system `mcl` on `PATH`.

This repo currently includes MCL version `14-137`, but users are welcome to
replace it with a newer upstream release from
[https://micans.org/mcl/](https://micans.org/mcl/).

## Building C Components

```bash
make build-c-engine
make build-indexed-orthologs
make build-indexed-inparalogs
make build-indexed-coorthologs
make build-indexed-rbh
```

Or build all current binaries:

```bash
make
```

## Testing

Run the test suite:

```bash
python3 -m unittest discover -s tests -v
```

Useful helper wrappers:

```bash
PYTHONPATH=src python3 scripts/compare_fixture.py --help
PYTHONPATH=src python3 scripts/compare_real_dataset.py --help
```

## Notes

The rewrite is already usable.

It also adds a strict `1:1` RBH output:

- `pairs/rbh_1to1.txt`

The rewrite does not try to exactly imitate every tiny floating-point score
difference from the old Perl/MySQL pipeline. The focus is correct relationships,
better runtime, simpler I/O, and easier HPC use.

`build/` is generated output and can be removed and rebuilt.
`src/mcl/` is reserved for vendored upstream MCL source.
`scripts/` contains wrappers; reusable code lives in `src/`.
