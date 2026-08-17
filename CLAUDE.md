# mef90 / vDef

Fortran + PETSc reference implementation of the variational approach to fracture. See `README.md` for the
scientific references and `INSTALL.md` for a full from-scratch install.

## Environment

Every build and run needs these set (loading the `petsc/main-g` module does it):

`PETSC_DIR` and `PETSC_ARCH` are needed. Ask the user for the proper value if they are not set

`MEF90_DIR` is required — the Makefiles and the test option files reference it directly. PETSc source, if
needed for reference, is at `/opt/HPC/src/petsc/main`.

Build products land in `${PETSC_ARCH}/objs` (`.o`/`.mod`) and `${PETSC_ARCH}/bin` (executables).

## Building

```sh
make vDef              # main program
make ThermoElasticity  # second program to keep working
```

`make` with no target builds everything in `DIRS` (MEF90 m_HeatXfer HeatXfer m_DefMech vDef Utils).

**`ThermoElasticity`, `vDefP`, `vDefUpa` and `WorkControlled` are currently broken** — do not build or test them, and do not report their failures as regressions. They are not in the default `vDef` target (see the commented-out line in `Makefile`); they only build
if invoked explicitly.

## Documentation

`make ford` generates the browsable source documentation into `doc/html` (config in `ford.md`),
creating the `.venv-ford` virtualenv on first use; `make fordclean` removes the output. `make doc`
is unrelated — it builds the LaTeX user manual `doc/vDef.pdf`.

## Testing

Standard smoke test — delete the result file first, since vDef will not overwrite an existing one:

```sh
rm -f junk.exo
mpirun -np 4 ${MEF90_DIR}/${PETSC_ARCH}/bin/vDef \
  -geometry ${MEF90_DIR}/TestMeshes/Beam1x10-tri.gen \
  -result junk.exo \
  -options_file vDef/data/test7.yaml
```

Regression suites (each runs the binary and diffs against stored output in `<dir>/results/`):

```sh
make -C vDef test              # 2D + 3D, test7.yaml and test8.yaml, NP=2
make -C ThermoElasticity test  # 2D + 3D, order 1 and 2, NP=4
```

Top-level `make test` also runs `HeatXfer`. `make runtests` runs the unit tests under `Tests/`.

## Layout

- `MEF90/` — core library: DMPlex wrappers, elements, Hooke's law, materials, EXO I/O, contexts.
- `m_DefMech/` — deformation mechanics: AT models (AT1/AT2/KKL/LinSoft), energy splits, plasticity.
- `m_HeatXfer/` — heat transfer module.
- `vDef/`, `ThermoElasticity/`, `HeatXfer/` — programs; each has `data/` (YAML option files) and `results/`
  (expected test output).
- `TestMeshes/` — shared `.gen`/`.msh` meshes used by all test targets.
- `Utils/`, `python/`, `bin/` — helper tools and scripts.

Runtime configuration is PETSc options plus YAML option files passed with `-options_file`.

## Style

Fortran free-form. `fortitude.toml` sets a 132-column limit for linting; `fprettify.cfg` holds the formatter
settings (4-space whitespace, lowercase keywords). Match the surrounding file.
