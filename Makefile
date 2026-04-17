PYTHON ?= python3
CC ?= cc
CFLAGS ?= -O3 -std=c11 -Wall -Wextra

.PHONY: all clean test build-cli build-c-engine build-indexed-orthologs build-indexed-inparalogs build-indexed-coorthologs build-indexed-rbh build-mcl

all: build/orthomclx build/pairs_engine build/indexed_orthologs build/indexed_inparalogs build/indexed_coorthologs build/indexed_rbh

build/orthomclx: src/orthomcl/cli.py
	mkdir -p build
	printf '%s\n' '#!/usr/bin/env bash' 'set -euo pipefail' 'ROOT=$$(cd "$$(dirname "$$0")/.." && pwd)' 'export PYTHONPATH="$$ROOT/src$${PYTHONPATH:+:$$PYTHONPATH}"' 'exec "$(PYTHON)" -m orthomcl.cli "$$@"' > $@
	chmod +x $@

build-cli: build/orthomclx

build/pairs_engine: src/c/pairs_engine.c
	mkdir -p build
	$(CC) $(CFLAGS) $< -lm -o $@

build-c-engine: build/pairs_engine

build/indexed_orthologs: src/c/indexed_orthologs.c
	mkdir -p build
	$(CC) $(CFLAGS) $< -lm -o $@

build-indexed-orthologs: build/indexed_orthologs

build/indexed_inparalogs: src/c/indexed_inparalogs.c
	mkdir -p build
	$(CC) $(CFLAGS) $< -lm -o $@

build-indexed-inparalogs: build/indexed_inparalogs

build/indexed_coorthologs: src/c/indexed_coorthologs.c
	mkdir -p build
	$(CC) $(CFLAGS) $< -lm -o $@

build-indexed-coorthologs: build/indexed_coorthologs

build/indexed_rbh: src/c/indexed_rbh.c
	mkdir -p build
	$(CC) $(CFLAGS) $< -lm -o $@

build-indexed-rbh: build/indexed_rbh

build-mcl:
	mkdir -p build
	@if [ -x src/mcl/configure ] || [ -f src/mcl/Makefile ] || [ -f src/mcl/makefile ]; then \
		echo "Building vendored MCL from src/mcl"; \
		$(MAKE) -C src/mcl; \
		if [ -x src/mcl/bin/mcl ]; then cp src/mcl/bin/mcl build/mcl; \
		elif [ -x src/mcl/mcl ]; then cp src/mcl/mcl build/mcl; \
		else echo "Built src/mcl but could not find the mcl executable."; exit 1; fi; \
	else \
		echo "No vendored MCL source found in src/mcl."; \
		echo "Place the upstream MCL source tree in src/mcl/ and rerun make build-mcl."; \
		exit 1; \
	fi

test:
	$(PYTHON) -m unittest discover -s tests -v

clean:
	rm -rf build tmp/compat
