---
title: Maintaining these pages
---

How this documentation is written, built, and published. None of it is needed
to *read* the documentation — see the source files, modules, and procedures in
the navigation bar for that.

## Doc-comment convention

The `!!!` banner blocks placed immediately *before* a module, type, or
procedure definition are picked up as its documentation (`predocmark: !!`), as
are `!!!` comments preceding a declaration. Trailing doc comments, if ever
needed, use `!>` (`docmark: >`). Comments that must *not* appear in the
documentation — section headers inside procedure bodies, implementation notes —
use `!!` or `!`, which FORD ignores.

Authorship is recorded as FORD metadata rather than a copyright banner, one
line per author, and must come first in the docstring:

```fortran
!!! author: Blaise Bourdin (2020, bourdin@lsu.edu)
!!! author: Blaise Bourdin (2026, bourdin@mcmaster.ca)
!!!
!!!  aAT1: the "a" function of the standard AT1 model, i.e. $a(\alpha) = (1-\alpha)^2$
!!!
```

FORD only reads metadata from the opening lines of a docstring, so an
`author:` line placed after the description is silently treated as prose.

## Building

```
make ford
```

creates `.venv-ford` on first use, installs FORD and the `doi_link` extension
into it, generates `doc/html`, and cross-references the result. `make fordclean`
removes the output. Note that `make doc` is a different target: it builds the
LaTeX user manual `doc/vDef.pdf`.

Pass `PETSC_SRC=<petsc source tree>` to link the PETSc calls in the source
listings to their manual pages; the default is `PETSC_DIR`, which only holds
sources when PETSc was built in place.

The equivalent by hand is

```
python3 -m venv .venv-ford && .venv-ford/bin/pip install ford ./doc/doi-link
.venv-ford/bin/ford ford.md
```

## Publishing

```
make fordpublish
```

rebuilds the pages and pushes them to the repository that serves them, checked
out as an ordinary clone next to this one (`MEF90_DOCS`, cloned from
`MEF90_DOCS_URL`). The site is served from that repository by GitHub Pages.

`MEF90_DOCS_KEEP` holds back the files belonging to that repository rather than
to the generated site — its `README`, `LICENSE`, and the `CNAME` that GitHub
writes there when a custom domain is set. Without it the publish step would
delete them and take the custom domain down.

## Preprocessor templates

Several files are compiled multiple times by the Makefiles with different macros
(`MEF90_DIM=2/3`, `MEF90_ELEMENTTYPE=...`). FORD preprocesses each file once, so
these pages document a single representative instantiation: `MEF90_DIM=3` with
`MEF90_ELEMENTTYPE=MEF90Element3DVect`. References to the other instantiations
(e.g. `*_MEF90Element2DScal`) appear unlinked.

## Cross-referencing the source listings

FORD highlights the source files it renders with Pygments and does nothing else
to them: a call goes nowhere, and Pygments gives user-defined names no colour at
all. So `make ford` finishes by running `bin/ford-xref`, which

  - links each identifier naming a documented module, type, or procedure to its
    page, and re-tags it so the stylesheet already shipped with FORD colours it
    as a function, class, or namespace;
  - links PETSc calls to their manual pages on petsc.org. The section of those
    URLs is not derivable from the routine name — `PetscSectionGetOffset` is
    filed under `PetscSection`, `PetscViewerASCIIPrintf` under `Viewer`,
    `PetscPrintf` under `Sys` — so it is read out of the PETSc sources the way
    PETSc's own `doc/build_man_pages.py` reads it, from `SUBMANSEC`/`MANSEC`;
  - re-tags the PETSc error-checking macros (`PetscCall`, `SETERRQ`, ...), which
    are preprocessor macros rather than Fortran keywords and so are invisible to
    the lexer even though they open nearly every statement in this codebase.

Because Fortran is case insensitive and the listings are only syntax
highlighted, a local variable sharing a name with an entity is linked too;
definition sites are left alone so that a procedure does not link to itself.

`bin/ford-xref` also resets the underline FORD's line-number anchors would
otherwise paint across every cross-referenced line, in `doc/user.css`.

## DOI identifiers

DOIs written as plain text (`doi:10.xxxx/...`) are turned into links to
`https://doi.org/...` automatically by the `doi_link` markdown extension in
`doc/doi-link/`.
