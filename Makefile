# Native build/install and browser demo packaging.
CARGO       ?= cargo
CC          ?= cc
INSTALL     ?= install
PREFIX      ?= /usr/local
DESTDIR     ?=
BINDIR      ?= $(PREFIX)/bin
LIBDIR      ?= $(PREFIX)/lib
INCLUDEDIR  ?= $(PREFIX)/include
PKGCONFIGDIR ?= $(LIBDIR)/pkgconfig

HOST_OS      := $(shell uname -s)
SHARED_EXT   ?= $(if $(filter Darwin,$(HOST_OS)),dylib,so)
NATIVE_BIN   := target/release/qs-factor
NATIVE_STATIC := target/release/librusqsieve.a
NATIVE_SHARED := target/release/librusqsieve.$(SHARED_EXT)
VERSION      := $(shell sed -n 's/^version = "\(.*\)"/\1/p' Cargo.toml)

WASM_TARGET := wasm32-unknown-unknown
WASM_SCALAR := target/wasm-scalar/$(WASM_TARGET)/release/rusqsieve.wasm
WASM_SIMD   := target/wasm-simd/$(WASM_TARGET)/release/rusqsieve.wasm
DOCS        := docs
WEB         := web
# Derived from the directory, not hand-maintained: this list silently lost
# coordinator.js when the relation coordinator moved into its own Worker, which
# 404s on GitHub Pages and leaves boot() awaiting a "ready" message that never
# arrives. `serve.mjs` is the local preview server and is deliberately not
# published; it is excluded by extension (`*.js` does not match `*.mjs`).
ASSETS      := $(notdir $(wildcard $(WEB)/*.html) $(wildcard $(WEB)/*.css) \
	$(wildcard $(WEB)/*.js))
DOCS_FILES  := $(addprefix $(DOCS)/,$(ASSETS)) $(DOCS)/rusqsieve.wasm \
	$(DOCS)/rusqsieve-simd.wasm $(DOCS)/.nojekyll

.DEFAULT_GOAL := native
.PHONY: native install docs docs-verify wasm serve test c-api-smoke clean

native:
	$(CARGO) build --release

install: native
	$(INSTALL) -d "$(DESTDIR)$(BINDIR)" "$(DESTDIR)$(LIBDIR)" \
		"$(DESTDIR)$(INCLUDEDIR)" "$(DESTDIR)$(PKGCONFIGDIR)"
	$(INSTALL) -m 0755 "$(NATIVE_BIN)" "$(DESTDIR)$(BINDIR)/qs-factor"
	$(INSTALL) -m 0755 "$(NATIVE_SHARED)" \
		"$(DESTDIR)$(LIBDIR)/$(notdir $(NATIVE_SHARED))"
	$(INSTALL) -m 0644 "$(NATIVE_STATIC)" \
		"$(DESTDIR)$(LIBDIR)/$(notdir $(NATIVE_STATIC))"
	$(INSTALL) -m 0644 rusqsieve.h "$(DESTDIR)$(INCLUDEDIR)/rusqsieve.h"
	sed -e 's|@PREFIX@|$(PREFIX)|g' \
		-e 's|@LIBDIR@|$(LIBDIR)|g' \
		-e 's|@INCLUDEDIR@|$(INCLUDEDIR)|g' \
		-e 's|@VERSION@|$(VERSION)|g' \
		rusqsieve.pc.in > "$(DESTDIR)$(PKGCONFIGDIR)/rusqsieve.pc"
	chmod 0644 "$(DESTDIR)$(PKGCONFIGDIR)/rusqsieve.pc"

# Every same-directory asset the published frontend references must exist in
# docs/. Deriving ASSETS from the directory prevents the omission that broke
# GitHub Pages; this catches the reverse case, a reference to a file that was
# renamed or never added to web/ at all.
# Intentionally has no prerequisites: it must audit `docs/` exactly as it stands,
# so it also catches a tree that was committed without re-running `make docs`.
docs-verify:
	@missing=0; \
	for src in $(DOCS)/*.js $(DOCS)/*.html; do \
		for ref in $$(grep -oE '"\./[A-Za-z0-9_.-]+"' "$$src" | tr -d '"'); do \
			if [ ! -f "$(DOCS)/$${ref#./}" ]; then \
				echo "$$src references $$ref, absent from $(DOCS)/"; \
				missing=1; \
			fi; \
		done; \
	done; \
	if [ $$missing -ne 0 ]; then echo "docs/ is incomplete"; exit 1; fi; \
	echo "docs/: all $$(ls $(DOCS) | wc -l | tr -d ' ') files present, every reference resolves."

docs: $(DOCS_FILES) docs-verify
	@echo "docs/ ready for GitHub Pages (scalar $$(ls -lh $(DOCS)/rusqsieve.wasm | awk '{print $$5}'), SIMD $$(ls -lh $(DOCS)/rusqsieve-simd.wasm | awk '{print $$5}'))."
	@echo "  Local preview:  make serve"
	@echo "  Publish:        Settings > Pages > Deploy from branch > /docs"

wasm: $(WASM_SCALAR) $(WASM_SIMD)

$(WASM_SCALAR): $(shell find src -name '*.rs') Cargo.toml
	$(CARGO) build --release --target $(WASM_TARGET) --target-dir target/wasm-scalar --lib --no-default-features

$(WASM_SIMD): $(shell find src -name '*.rs') Cargo.toml
	RUSTFLAGS="$(RUSTFLAGS) -C target-feature=+simd128" \
		$(CARGO) build --release --target $(WASM_TARGET) --target-dir target/wasm-simd --lib --no-default-features --features wasm-simd128

# Keep LLVM's speed-optimized output intact. Binaryen 120's -O3 and -Oz both
# regress the measured 192-bit sieve by roughly 50%, despite saving about 20 KiB.
$(DOCS)/rusqsieve.wasm: $(WASM_SCALAR)
	@mkdir -p $(DOCS)
	cp $< $@

$(DOCS)/rusqsieve-simd.wasm: $(WASM_SIMD)
	@mkdir -p $(DOCS)
	cp $< $@

$(DOCS)/%: $(WEB)/%
	@mkdir -p $(DOCS)
	cp $< $@

$(DOCS)/.nojekyll:
	@mkdir -p $(DOCS)
	@touch $@

serve: docs
	node $(WEB)/serve.mjs $(DOCS) 8000

# Native + wasm correctness checks. The browser architecture check drives the
# coordinator-Worker/sieve-Worker protocol on node worker threads; the focused
# glue check covers malformed packets and invalid command sequences.
c-api-smoke: native
	$(CC) -std=c11 -Wall -Wextra -Werror -I. tests/c_api_smoke.c \
		-Ltarget/release -lrusqsieve -Wl,-rpath,'$$ORIGIN/release' \
		-o target/c-api-smoke
	target/c-api-smoke

test:
	$(CARGO) test
	$(MAKE) c-api-smoke
	$(MAKE) wasm
	$(MAKE) docs-verify
	node tools/browser-arch-check.mjs
	node tools/browser-glue-failure-check.mjs

clean:
	rm -rf $(DOCS)
