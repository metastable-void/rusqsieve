#!/bin/sh

set -eu

usage() {
    cat <<'EOF'
Usage: ./install.sh [--prefix PATH] [--destdir PATH]

Install rusqsieve. PREFIX defaults to /usr/local and may also be set in the
environment. DESTDIR is supported for staged packaging. If the destination is
not writable, the installer asks sudo to re-run it.
EOF
}

prefix=${PREFIX:-/usr/local}
destdir=${DESTDIR:-}
elevated=no

while [ "$#" -gt 0 ]; do
    case "$1" in
        --prefix)
            [ "$#" -ge 2 ] || {
                printf '%s\n' 'install.sh: --prefix requires a path' >&2
                exit 2
            }
            prefix=$2
            shift 2
            ;;
        --prefix=*)
            prefix=${1#*=}
            shift
            ;;
        --destdir)
            [ "$#" -ge 2 ] || {
                printf '%s\n' 'install.sh: --destdir requires a path' >&2
                exit 2
            }
            destdir=$2
            shift 2
            ;;
        --destdir=*)
            destdir=${1#*=}
            shift
            ;;
        --elevated)
            elevated=yes
            shift
            ;;
        -h | --help)
            usage
            exit 0
            ;;
        *)
            printf 'install.sh: unknown argument: %s\n' "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

case "$prefix" in
    /*) ;;
    *)
        printf '%s\n' 'install.sh: prefix must be an absolute path' >&2
        exit 2
        ;;
esac

case "$destdir" in
    '' | /*) ;;
    *)
        printf '%s\n' 'install.sh: destdir must be empty or an absolute path' >&2
        exit 2
        ;;
esac

script_dir=$(
    CDPATH='' cd -P "$(dirname "$0")" >/dev/null 2>&1
    pwd
)
self=$script_dir/$(basename "$0")
install_root=$destdir$prefix

writable_ancestor=$install_root
while [ ! -e "$writable_ancestor" ]; do
    parent=$(dirname "$writable_ancestor")
    [ "$parent" != "$writable_ancestor" ] || break
    writable_ancestor=$parent
done

if [ ! -w "$writable_ancestor" ] && [ "$(id -u)" -ne 0 ]; then
    [ "$elevated" = no ] || {
        printf 'install.sh: destination remains unwritable: %s\n' "$install_root" >&2
        exit 1
    }
    command -v sudo >/dev/null 2>&1 || {
        printf 'install.sh: %s is not writable and sudo is unavailable\n' \
            "$install_root" >&2
        exit 1
    }
    printf 'Administrator access is required to install under %s.\n' "$prefix"
    exec sudo -- "$self" \
        --elevated \
        --prefix "$prefix" \
        --destdir "$destdir"
fi

install -d \
    "$install_root/bin" \
    "$install_root/lib" \
    "$install_root/lib/pkgconfig" \
    "$install_root/include"

for source in "$script_dir"/bin/*; do
    [ -f "$source" ] || continue
    install -m 0755 "$source" "$install_root/bin/$(basename "$source")"
done

for source in "$script_dir"/lib/*; do
    [ -f "$source" ] || continue
    case "$source" in
        *.so | *.dylib)
            mode=0755
            ;;
        *)
            mode=0644
            ;;
    esac
    install -m "$mode" "$source" "$install_root/lib/$(basename "$source")"
done

install -m 0644 "$script_dir/include/rusqsieve.h" \
    "$install_root/include/rusqsieve.h"

escaped_prefix=$(printf '%s' "$prefix" | sed 's/[\\&|]/\\&/g')
pc_tmp=$(mktemp "${TMPDIR:-/tmp}/rusqsieve.pc.XXXXXXXX")
trap 'rm -f "$pc_tmp"' EXIT HUP INT TERM
sed "s|^prefix=.*|prefix=$escaped_prefix|" \
    "$script_dir/lib/pkgconfig/rusqsieve.pc" >"$pc_tmp"
install -m 0644 "$pc_tmp" "$install_root/lib/pkgconfig/rusqsieve.pc"

printf 'Installed rusqsieve under %s\n' "$install_root"
case "$(uname -s)" in
    Linux | FreeBSD)
        printf '%s\n' \
            "If $prefix/lib is not already on the loader path, configure it or set LD_LIBRARY_PATH."
        ;;
esac
