#!/bin/sh
#
# Teach fortls -- the Fortran language server behind editor plugins such as the
# VS Code "Modern Fortran" extension -- about PETSc, so that PETSc types and
# procedures get completion, hover and go-to-definition inside mef90 sources.
#
# Two artifacts are generated.  Both are machine- and PETSC_ARCH-specific, both
# are covered by .gitignore, and both should be regenerated after reconfiguring
# or updating PETSc:
#
#   ${MEF90_DIR}/${PETSC_ARCH}/petsc-ftn-mod/*.F90
#       cpp-expanded copies of PETSc's Fortran module sources.  PETSc assembles
#       petscsys/petscvec/.../petscsnes out of #include'd .h90 fragments, and
#       fortls resolves #include only to harvest #define's -- it does not splice
#       the included Fortran back in.  Without this it sees the module *names*
#       but none of their contents.  The directory sits inside the repository,
#       where fortls's own recursive scan finds it, so no source_dirs setting is
#       required.  `make clean` does not remove it.
#
#   ${MEF90_DIR}/.fortlsrc
#       pp_defs for the PETSc Fortran type macros (PetscReal, PetscInt, Vec,
#       Mat, SNES, ...), so that variables declared with them are typed rather
#       than unknown.  fortls (as of 3.2.2) matches include paths with the
#       regex [\"\w\.]*, which accepts neither <angle brackets> nor a '/', so it
#       can never open petsc/finclude/petsc.h no matter how include_dirs is set.
#       Handing it the already-expanded macros sidesteps that.
#
# Requires MEF90_DIR, PETSC_DIR and PETSC_ARCH.  When PETSc was installed with
# --prefix, its Fortran module sources are not part of the installation; set
# PETSC_SRC to the source tree that was configured with this PETSC_ARCH.

set -e

: "${MEF90_DIR:?set MEF90_DIR}"
: "${PETSC_DIR:?set PETSC_DIR}"
: "${PETSC_ARCH:?set PETSC_ARCH}"

# in-place builds keep petscvariables under ${PETSC_ARCH}, prefix installs do not
PETSC_VARIABLES=
for d in "${PETSC_DIR}/${PETSC_ARCH}" "${PETSC_DIR}"; do
    if [ -r "${d}/lib/petsc/conf/petscvariables" ]; then
        PETSC_VARIABLES=${d}/lib/petsc/conf/petscvariables
        break
    fi
done
if [ -z "${PETSC_VARIABLES}" ]; then
    echo "no petscvariables under ${PETSC_DIR} -- is PETSC_DIR/PETSC_ARCH correct?" >&2
    exit 1
fi

# preprocess with the very compilers PETSc was configured with
FC=${FC:-$(sed -n 's/^FC = //p' "${PETSC_VARIABLES}")}
CC=${CC:-$(sed -n 's/^PCC = //p' "${PETSC_VARIABLES}")}
: "${FC:?could not read FC from ${PETSC_VARIABLES}}"
: "${CC:?could not read PCC from ${PETSC_VARIABLES}}"

# the ftn-mod sources only exist in the source tree, never in a --prefix install
if [ -z "${PETSC_SRC}" ]; then
    if [ -d "${PETSC_DIR}/src/sys/ftn-mod" ]; then
        PETSC_SRC=${PETSC_DIR}
    else
        echo "${PETSC_DIR} has no src/*/ftn-mod (--prefix install?)" >&2
        echo "set PETSC_SRC to the PETSc source tree configured with PETSC_ARCH=${PETSC_ARCH}" >&2
        exit 1
    fi
fi
if [ ! -d "${PETSC_SRC}/${PETSC_ARCH}/include" ]; then
    echo "no ${PETSC_SRC}/${PETSC_ARCH}/include -- ${PETSC_SRC} was not configured with PETSC_ARCH=${PETSC_ARCH}" >&2
    exit 1
fi

# petscconf.h lives under ${PETSC_ARCH} in an in-place build, in include/ otherwise
INCLUDES="-I${PETSC_DIR}/include"
if [ -d "${PETSC_DIR}/${PETSC_ARCH}/include" ]; then
    INCLUDES="-I${PETSC_DIR}/${PETSC_ARCH}/include ${INCLUDES}"
fi

# 1. self-contained PETSc Fortran modules
OUT=${MEF90_DIR}/${PETSC_ARCH}/petsc-ftn-mod
rm -rf "${OUT}"
mkdir -p "${OUT}"
for f in "${PETSC_SRC}"/src/*/ftn-mod/*.F90; do
    ${FC} -cpp -E -P \
        -I"${PETSC_SRC}/${PETSC_ARCH}/include" \
        -I"${PETSC_SRC}/include" \
        "${f}" -o "${OUT}/$(basename "${f}")"
done
echo "wrote $(ls "${OUT}" | wc -l | tr -d ' ') modules to ${OUT}"

# 2. .fortlsrc holding the PETSc Fortran type macros
CC="${CC}" python3 - "${MEF90_DIR}/.fortlsrc" ${INCLUDES} <<'PYEOF'
import json, os, re, subprocess, sys

out_path, includes = sys.argv[1], sys.argv[2:]
cc = os.environ["CC"]
header = '#include "petsc/finclude/petsc.h"\n'


def cpp(flags, stdin):
    args = [cc, "-E", *flags, "-w", *includes, "-x", "c", "-"]
    return subprocess.run(args, input=stdin, capture_output=True, text=True).stdout


# every object-like macro reachable from petsc/finclude/petsc.h
dumped = cpp(["-dM"], header)
names = sorted(
    {
        m.group(1)
        for m in (re.match(r"#define\s+([A-Za-z_]\w*)(\()?", l) for l in dumped.splitlines())
        if m and not m.group(2)
    }
)
if not names:
    sys.exit(f"{cc} found no macros in petsc/finclude/petsc.h with {' '.join(includes)}")

# expand each one; the marker is quoted so that it is not substituted itself
probe = header + "".join(f'"@@{n}@@" {n}\n' for n in names)
expanded = cpp(["-P"], probe)

is_type = re.compile(r"^(type\s*\(|class\s*\(|integer|real|logical|complex|character|double)", re.I)
pp_defs = {}
for line in expanded.splitlines():
    m = re.match(r'"@@(\w+)@@"\s+(.*)$', line.strip())
    if not m:
        continue
    name, value = m.group(1), m.group(2).strip()
    if value and value != name and is_type.match(value):
        pp_defs[name] = value

with open(out_path, "w") as f:
    json.dump({"pp_defs": pp_defs}, f, indent=4, sort_keys=True)
    f.write("\n")
print(f"wrote {len(pp_defs)} pp_defs to {out_path}")
PYEOF
