# Native build/install and browser demo packaging.
CARGO       ?= cargo
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
ASSETS      := index.html index.css abi.js numtheory.js worker.js index.js
DOCS_FILES  := $(addprefix $(DOCS)/,$(ASSETS)) $(DOCS)/rusqsieve.wasm \
	$(DOCS)/rusqsieve-simd.wasm $(DOCS)/.nojekyll

.DEFAULT_GOAL := native
.PHONY: native install docs wasm serve test clean

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

docs: $(DOCS_FILES)
	@echo "docs/ ready for GitHub Pages (scalar $$(ls -lh $(DOCS)/rusqsieve.wasm | awk '{print $$5}'), SIMD $$(ls -lh $(DOCS)/rusqsieve-simd.wasm | awk '{print $$5}'))."
	@echo "  Local preview:  make serve"
	@echo "  Publish:        Settings > Pages > Deploy from branch > /docs"

wasm: $(WASM_SCALAR) $(WASM_SIMD)

$(WASM_SCALAR): $(shell find src -name '*.rs') Cargo.toml
	$(CARGO) build --release --target $(WASM_TARGET) --target-dir target/wasm-scalar --lib --no-default-features

$(WASM_SIMD): $(shell find src -name '*.rs') Cargo.toml
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

# Native + wasm correctness checks.
test:
	$(CARGO) test
	$(MAKE) wasm

clean:
	rm -rf $(DOCS)
