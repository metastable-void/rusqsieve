#!/bin/sh
#
# Build and package rusqsieve releases. With no arguments every supported
# target is built. Pass one or more target triples to build a subset.
#
# Archives are written to $OUT_DIR (the repository root by default).

set -eu

SCRIPT_DIR=$(
    CDPATH='' cd -P "$(dirname "$0")" >/dev/null 2>&1
    pwd
)
MANIFEST="$SCRIPT_DIR/Cargo.toml"
HEADER="$SCRIPT_DIR/rusqsieve.h"
PC_TEMPLATE="$SCRIPT_DIR/rusqsieve.pc.in"
RELEASE_FILES="$SCRIPT_DIR/release"
OUT_DIR="${OUT_DIR:-$SCRIPT_DIR}"
BUILD_DIR="${CARGO_TARGET_DIR:-$SCRIPT_DIR/target}"
case "$OUT_DIR" in
    /*) ;;
    *) OUT_DIR="$SCRIPT_DIR/$OUT_DIR" ;;
esac
case "$BUILD_DIR" in
    /*) ;;
    *) BUILD_DIR="$SCRIPT_DIR/$BUILD_DIR" ;;
esac
readonly SCRIPT_DIR MANIFEST HEADER PC_TEMPLATE RELEASE_FILES OUT_DIR BUILD_DIR

ALL_TARGETS='
    x86_64-unknown-linux-gnu
    x86_64-unknown-linux-musl
    aarch64-unknown-linux-musl
    aarch64-unknown-linux-gnu
    x86_64-unknown-freebsd
    x86_64-pc-windows-msvc
    aarch64-apple-darwin
    wasm32-unknown-unknown
'
readonly ALL_TARGETS

die() {
    printf 'build-release.sh: %s\n' "$*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
Usage: ./build-release.sh [TARGET ...]

Build release archives for all supported targets, or only the listed targets:
  x86_64-unknown-linux-gnu
  x86_64-unknown-linux-musl
  aarch64-unknown-linux-musl
  aarch64-unknown-linux-gnu
  x86_64-unknown-freebsd
  x86_64-pc-windows-msvc
  aarch64-apple-darwin
  wasm32-unknown-unknown

Environment:
  SDKROOT           Required for aarch64-apple-darwin.
  OUT_DIR           Archive output directory (default: repository root).
  CARGO_TARGET_DIR  Cargo build directory (default: ./target).
  SOURCE_DATE_EPOCH Timestamp used for reproducible archives.
EOF
}

is_supported_target() {
    candidate=$1
    for supported_target in $ALL_TARGETS; do
        [ "$candidate" = "$supported_target" ] && return 0
    done
    return 1
}

need_command() {
    command -v "$1" >/dev/null 2>&1 ||
        die "required command '$1' was not found in PATH"
}

need_file() {
    [ -f "$1" ] || die "expected build artifact is missing: $1"
}

run_with_rustflags() (
    release_rustflags=$1
    shift
    unset CARGO_ENCODED_RUSTFLAGS
    RUSTFLAGS=$release_rustflags
    export RUSTFLAGS
    exec "$@"
)

copy_file() {
    mode=$1
    source=$2
    destination=$3
    need_file "$source"
    mkdir -p "$(dirname "$destination")"
    install -m "$mode" "$source" "$destination"
}

manifest_value() {
    key=$1
    sed -n "s/^${key}[[:space:]]*=[[:space:]]*\"\\([^\"]*\\)\".*/\\1/p" "$MANIFEST" |
        head -n 1
}

generate_pc() {
    pc_target=$1
    pc_prefix=$2
    pc_destination=$3

    case "$pc_target" in
        *-linux-gnu | *-linux-musl)
            private_libs='-ldl -lpthread -lm'
            ;;
        *-freebsd)
            private_libs='-lpthread -lm'
            ;;
        *-apple-darwin)
            private_libs='-framework Security -framework CoreFoundation -liconv -lSystem'
            ;;
        *-pc-windows-msvc)
            private_libs='-lws2_32 -luserenv -lbcrypt -lntdll'
            ;;
        *)
            die "no pkg-config system-library mapping for $pc_target"
            ;;
    esac

    # The single quotes intentionally preserve pkg-config's ${prefix}.
    # shellcheck disable=SC2016
    sed \
        -e "s|@PREFIX@|$pc_prefix|g" \
        -e 's|@LIBDIR@|${prefix}/lib|g' \
        -e 's|@INCLUDEDIR@|${prefix}/include|g' \
        -e "s|@VERSION@|$VERSION|g" \
        -e "s|^Libs\\.private:.*|Libs.private: $private_libs|" \
        "$PC_TEMPLATE" >"$pc_destination"
    chmod 0644 "$pc_destination"
}

add_common_files() {
    package_dir=$1
    copy_file 0644 "$HEADER" "$package_dir/include/rusqsieve.h"
    copy_file 0644 "$SCRIPT_DIR/LICENSE-APACHE" "$package_dir/LICENSE-APACHE"
    copy_file 0644 "$SCRIPT_DIR/LICENSE-MPL" "$package_dir/LICENSE-MPL"
    copy_file 0644 "$SCRIPT_DIR/README.md" "$package_dir/README.md"
}

add_posix_installer() {
    package_dir=$1
    copy_file 0755 "$RELEASE_FILES/install.sh" "$package_dir/install.sh"
}

add_windows_installer() {
    package_dir=$1
    copy_file 0755 "$RELEASE_FILES/install.bat" "$package_dir/install.bat"
    copy_file 0644 "$RELEASE_FILES/install.ps1" "$package_dir/install.ps1"
}

create_archive() {
    package_dir=$1
    archive_name="$(basename "$package_dir").tar.gz"
    archive_tmp="$WORK_DIR/$archive_name"

    mkdir -p "$OUT_DIR"
    find "$package_dir" -type d -exec chmod 0755 {} +
    if tar --version 2>/dev/null | grep -q 'GNU tar'; then
        tar \
            --sort=name \
            --mtime="@$SOURCE_DATE_EPOCH" \
            --owner=0 \
            --group=0 \
            --numeric-owner \
            -C "$WORK_DIR" \
            -czf "$archive_tmp" \
            "$(basename "$package_dir")"
    else
        tar \
            -C "$WORK_DIR" \
            -czf "$archive_tmp" \
            "$(basename "$package_dir")"
    fi
    mv -f "$archive_tmp" "$OUT_DIR/$archive_name"
    printf 'Created %s\n' "$OUT_DIR/$archive_name"
}

package_posix() {
    target=$1
    dynamic_suffix=$2
    include_dynamic=$3
    artifacts="$BUILD_DIR/$target/release"
    package_dir="$WORK_DIR/rusqsieve-$VERSION-$target"

    mkdir -p "$package_dir/bin" "$package_dir/lib/pkgconfig"
    copy_file 0755 "$artifacts/qs-factor" "$package_dir/bin/qs-factor"
    copy_file 0644 "$artifacts/librusqsieve.a" "$package_dir/lib/librusqsieve.a"
    if [ "$include_dynamic" = yes ]; then
        copy_file 0755 \
            "$artifacts/librusqsieve.$dynamic_suffix" \
            "$package_dir/lib/librusqsieve.$dynamic_suffix"
    fi
    add_common_files "$package_dir"
    generate_pc "$target" /usr/local "$package_dir/lib/pkgconfig/rusqsieve.pc"
    add_posix_installer "$package_dir"
    create_archive "$package_dir"
}

package_windows() {
    target=x86_64-pc-windows-msvc
    artifacts="$BUILD_DIR/$target/release"
    package_dir="$WORK_DIR/rusqsieve-$VERSION-$target"

    mkdir -p "$package_dir/bin" "$package_dir/lib/pkgconfig"
    copy_file 0755 "$artifacts/qs-factor.exe" "$package_dir/bin/qs-factor.exe"
    copy_file 0755 "$artifacts/rusqsieve.dll" "$package_dir/bin/rusqsieve.dll"
    copy_file 0644 "$artifacts/rusqsieve.lib" "$package_dir/lib/rusqsieve.lib"
    copy_file 0644 "$artifacts/rusqsieve.dll.lib" "$package_dir/lib/rusqsieve.dll.lib"
    add_common_files "$package_dir"
    generate_pc \
        "$target" \
        'C:/Program Files/rusqsieve' \
        "$package_dir/lib/pkgconfig/rusqsieve.pc"
    add_windows_installer "$package_dir"
    create_archive "$package_dir"
}

package_wasm() {
    target=wasm32-unknown-unknown
    scalar="$BUILD_DIR/release-wasm-scalar/$target/release/rusqsieve.wasm"
    simd="$BUILD_DIR/release-wasm-simd128/$target/release/rusqsieve.wasm"
    package_dir="$WORK_DIR/rusqsieve-$VERSION-$target"

    mkdir -p "$package_dir/web"
    for web_file in abi.js index.css index.html index.js numtheory.js serve.mjs worker.js; do
        copy_file 0644 "$SCRIPT_DIR/web/$web_file" "$package_dir/web/$web_file"
    done
    copy_file 0644 "$scalar" "$package_dir/web/rusqsieve.wasm"
    copy_file 0644 "$simd" "$package_dir/web/rusqsieve-simd.wasm"
    copy_file 0644 "$SCRIPT_DIR/LICENSE-APACHE" "$package_dir/LICENSE-APACHE"
    copy_file 0644 "$SCRIPT_DIR/LICENSE-MPL" "$package_dir/LICENSE-MPL"
    copy_file 0644 "$SCRIPT_DIR/README.md" "$package_dir/README.md"
    create_archive "$package_dir"
}

build_target() {
    target=$1

    printf '\nBuilding %s\n' "$target"
    case "$target" in
        x86_64-unknown-linux-gnu | aarch64-unknown-linux-gnu | x86_64-unknown-freebsd)
            run_with_rustflags "" \
                cross build --locked --release --target "$target" --target-dir "$BUILD_DIR"
            case "$target" in
                *-freebsd) package_posix "$target" so yes ;;
                *) package_posix "$target" so yes ;;
            esac
            ;;
        x86_64-unknown-linux-musl | aarch64-unknown-linux-musl)
            run_with_rustflags "" \
                cross build --locked --release --target "$target" --target-dir "$BUILD_DIR"
            package_posix "$target" so no
            ;;
        x86_64-pc-windows-msvc)
            run_with_rustflags "-C target-feature=+crt-static" \
                cargo xwin build \
                    --locked \
                    --release \
                    --target "$target" \
                    --target-dir "$BUILD_DIR"
            package_windows
            ;;
        aarch64-apple-darwin)
            run_with_rustflags "" \
                cargo zigbuild --locked --release --target "$target" --target-dir "$BUILD_DIR"
            package_posix "$target" dylib yes
            ;;
        wasm32-unknown-unknown)
            run_with_rustflags "" \
                cargo build \
                    --locked \
                    --release \
                    --target "$target" \
                    --target-dir "$BUILD_DIR/release-wasm-scalar" \
                    --lib \
                    --no-default-features
            run_with_rustflags "-C target-feature=+simd128" \
                cargo build \
                    --locked \
                    --release \
                    --target "$target" \
                    --target-dir "$BUILD_DIR/release-wasm-simd128" \
                    --lib \
                    --no-default-features \
                    --features wasm-simd128
            package_wasm
            ;;
        *)
            die "internal error: unhandled target $target"
            ;;
    esac
}

if [ "$#" -eq 1 ] && { [ "$1" = --help ] || [ "$1" = -h ]; }; then
    usage
    exit 0
fi

if [ "$#" -eq 0 ]; then
    TARGETS=$ALL_TARGETS
else
    TARGETS=$*
fi

for target in $TARGETS; do
    is_supported_target "$target" || {
        usage >&2
        die "unsupported target: $target"
    }
done

need_command cargo
need_command find
need_command install
need_command mv
need_command sed
need_command tar
need_command mktemp
need_file "$MANIFEST"
need_file "$HEADER"
need_file "$PC_TEMPLATE"

for target in $TARGETS; do
    case "$target" in
        x86_64-unknown-linux-* | aarch64-unknown-linux-* | x86_64-unknown-freebsd)
            need_command cross
            ;;
        x86_64-pc-windows-msvc)
            need_command cargo-xwin
            ;;
        aarch64-apple-darwin)
            need_command cargo-zigbuild
            [ -n "${SDKROOT:-}" ] ||
                die "SDKROOT must name an Apple SDK when building aarch64-apple-darwin"
            case "$SDKROOT" in
                /*) ;;
                *) die "SDKROOT must be an absolute path: $SDKROOT" ;;
            esac
            [ -d "$SDKROOT" ] || die "SDKROOT is not a directory: $SDKROOT"
            ;;
    esac
done

VERSION=$(manifest_value version)
readonly VERSION
[ -n "$VERSION" ] || die "could not read the package version from Cargo.toml"

if [ -z "${SOURCE_DATE_EPOCH:-}" ]; then
    if command -v git >/dev/null 2>&1 &&
        SOURCE_DATE_EPOCH="$(git -C "$SCRIPT_DIR" log -1 --format=%ct 2>/dev/null)" &&
        [ -n "$SOURCE_DATE_EPOCH" ]; then
        :
    else
        SOURCE_DATE_EPOCH=0
    fi
fi
readonly SOURCE_DATE_EPOCH

mkdir -p "$OUT_DIR"
WORK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/rusqsieve-release.XXXXXXXX")
readonly WORK_DIR
trap 'rm -rf "$WORK_DIR"' EXIT HUP INT TERM

cd "$SCRIPT_DIR"
for target in $TARGETS; do
    build_target "$target"
done
