# Build the browser demo (GitHub Pages) into docs/. `make docs` is the default.
WASM   := target/wasm32-unknown-unknown/release/rusqsieve.wasm
DOCS   := docs
WEB    := web
ASSETS := index.html index.css abi.js numtheory.js worker.js index.js
DOCS_FILES := $(addprefix $(DOCS)/,$(ASSETS)) $(DOCS)/rusqsieve.wasm $(DOCS)/.nojekyll

.DEFAULT_GOAL := docs
.PHONY: docs wasm serve test clean

docs: $(DOCS_FILES)
	@echo "docs/ ready for GitHub Pages (wasm $$(ls -lh $(DOCS)/rusqsieve.wasm | awk '{print $$5}'))."
	@echo "  Local preview:  make serve"
	@echo "  Publish:        Settings > Pages > Deploy from branch > /docs"

wasm: $(WASM)

$(WASM): $(shell find src -name '*.rs') Cargo.toml
	cargo build --release --target wasm32-unknown-unknown --lib --no-default-features

# Copy the wasm, running wasm-opt -Oz when available to shrink it.
$(DOCS)/rusqsieve.wasm: $(WASM)
	@mkdir -p $(DOCS)
	@if command -v wasm-opt >/dev/null 2>&1; then \
	  wasm-opt -Oz --enable-bulk-memory --enable-nontrapping-float-to-int -o $@ $< && echo "wasm-opt -Oz -> $@"; \
	else cp $< $@; fi

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
	cargo build --release --target wasm32-unknown-unknown --lib --no-default-features

clean:
	rm -rf $(DOCS)
