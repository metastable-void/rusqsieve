# Build the browser demo (GitHub Pages) into docs/. `make docs` is the default.
WASM_TARGET := wasm32-unknown-unknown
WASM_SCALAR := target/wasm-scalar/$(WASM_TARGET)/release/rusqsieve.wasm
WASM_SIMD   := target/wasm-simd/$(WASM_TARGET)/release/rusqsieve.wasm
DOCS        := docs
WEB         := web
ASSETS      := index.html index.css abi.js numtheory.js worker.js index.js
DOCS_FILES  := $(addprefix $(DOCS)/,$(ASSETS)) $(DOCS)/rusqsieve.wasm \
	$(DOCS)/rusqsieve-simd.wasm $(DOCS)/.nojekyll

.DEFAULT_GOAL := docs
.PHONY: docs wasm serve test clean

docs: $(DOCS_FILES)
	@echo "docs/ ready for GitHub Pages (scalar $$(ls -lh $(DOCS)/rusqsieve.wasm | awk '{print $$5}'), SIMD $$(ls -lh $(DOCS)/rusqsieve-simd.wasm | awk '{print $$5}'))."
	@echo "  Local preview:  make serve"
	@echo "  Publish:        Settings > Pages > Deploy from branch > /docs"

wasm: $(WASM_SCALAR) $(WASM_SIMD)

$(WASM_SCALAR): $(shell find src -name '*.rs') Cargo.toml
	cargo build --release --target $(WASM_TARGET) --target-dir target/wasm-scalar --lib --no-default-features

$(WASM_SIMD): $(shell find src -name '*.rs') Cargo.toml
	cargo build --release --target $(WASM_TARGET) --target-dir target/wasm-simd --lib --no-default-features --features wasm-simd128

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
	cargo test
	$(MAKE) wasm

clean:
	rm -rf $(DOCS)
