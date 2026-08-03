# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

SKIRT is a Monte Carlo radiative transfer C++ code for dusty astrophysical systems (`SKIRT/`). It's built on top of SMILE (`SMILE/`), a generic, SKIRT-independent framework for declaring a class hierarchy's properties via metadata and auto-deriving a schema from it — used both to parse/write `.ski` parameter files and to drive `MakeUp` (`MakeUp/`), the Qt desktop Q&A wizard for building `.ski` files interactively. All three share the same schema, so the CLI parser and the wizard enforce identical rules.

## Build

The build lives in a sibling directory `../release` (i.e. next to `git/`, not inside it). From the repo root:

```bash
./configSKIRT.sh              # configure the cmake build (only needed once, or after changing build options)
./makeSKIRT.sh <threads>      # build, e.g. ./makeSKIRT.sh 8
```

`makeSKIRT.sh` builds `framework`, the `skirt` CLI (`release/SKIRT/main/skirt`), and, if enabled, `MakeUp.app`. It also regenerates the SMILE schema file (`release/schemas/skirt.smile`) as part of the build — a malformed SMILE conditional-value-expression (e.g. in an `ATTRIBUTE_TYPE_INSERT`/`ATTRIBUTE_INSERT`/`ATTRIBUTE_TYPE_ALLOWED_IF` string) surfaces as a build failure at that step, so a clean build is a reasonable sanity check for schema-metadata edits, not just C++ changes.

To adjust a build option: `./configSKIRT.sh <OPTION_NAME>=<value>`. Relevant options (see top-level `CMakeLists.txt`, `SKIRT/mpi/CMakeLists.txt`): `BUILD_MAKE_UP` (requires Qt5/Qt6), `BUILD_WITH_MPI`, `WARNINGS_AS_ERRORS` (CI builds with this on), `BUILD_DOX_STYLE`. SMILE-only options `BUILD_SMILE_TOOL`/`BUILD_SMILE_SHAPES` build small standalone SMILE example/console tools.

Simulation resource files (dust properties, SED libraries, etc.) are not checked into the repo; `./downloadResources.sh` fetches them into a sibling `resources` directory next to `release`.

## Formatting

```bash
./formatSourceCode.sh
```

Requires clang-format **version 18.1** exactly (looks for `clang-format-18` then `clang-format` on `PATH`; refuses to run against any other version — CI uses `clang-format-18` on Ubuntu 24.04 and fails the build on any diff). Only `.hpp`/`.cpp` are reformatted (not `.h`/`.c`/`.cc`, to leave vendored third-party code alone).

## Typical edit workflow

1. Make the code change.
2. `./makeSKIRT.sh <threads>` — confirm it builds cleanly.
3. `./formatSourceCode.sh` — confirm clang-format needed no changes (or accept the ones it makes).
4. Re-run `./makeSKIRT.sh <threads>` if formatting changed anything, to reconfirm.
5. `git diff` / `git status` to review before committing.

## Running a simulation

```bash
release/SKIRT/main/skirt -o <output-dir> <path-to-ski-file>
```

The output directory must already exist (`skirt` will not create it). The `.ski` file is validated against the SMILE schema at parse time — this enforces the same `ATTRIBUTE_TYPE_ALLOWED_IF`/`ATTRIBUTE_RELEVANT_IF` metadata rules as the interactive wizard, so a schema-metadata bug can make an otherwise-valid `.ski` file unparseable (or vice versa, let through a combination the runtime then rejects) independently of the C++ runtime logic — when changing behavior that depends on item metadata, both layers may need updating (see Architecture below).

There is no in-repo unit test suite or `ctest` target. Correctness is checked with a separate functional/regression suite of 800+ `.ski` test cases with reference output (in a sibling `Functional9` directory, see below), run via the Python Toolkit for SKIRT (PTS): `pts test_fun <pattern>` (e.g. `pts test_fun *Disk*`) or `pts test_fun .` for the full suite, `pts endorse_fun <pattern>` to accept new reference output after verifying it. Ask before running or modifying anything under `Functional9` — it's not part of this repository.

## Documentation and project layout

The full user/developer/administrator documentation published on the SKIRT project web site is maintained as Doxygen-style source text in a separate sibling repository, normally checked out at `../../Web9/git`. It's worth consulting while working in this repo, in particular:

- `root/text/33-DeveloperGuide/DevCoding.txt` — design principles and coding style (naming conventions, formatting, preprocessor/`using`/pointer/container/exception-handling conventions, etc.) beyond what `clang-format` enforces automatically.
- `root/text/33-DeveloperGuide/DevGitHubFlow.txt` — the fork-and-pull contribution workflow. Key point: pull requests must target `SKIRT/SKIRT9:master` (the upstream repo), never the personal fork's own `master` — opening one against the fork by mistake still merges it there instead of upstream and requires manually re-opening against upstream plus resyncing the fork afterward.
- `root/text/33-DeveloperGuide/DevItems.txt` and `DevSmile.txt` — the canonical description of the `SimulationItem`/SMILE metadata system summarized below.
- `root/text/34-AdministratorGuide/AdminFunTests.txt` — the functional/regression test procedure referenced above.

On a full administrator checkout, this repo (`SKIRT9/git`) is one of several sibling project directories under `~/SKIRT`, alongside `PTS9` (the Python toolkit, providing the `pts` command used for functional tests), `Web9` (this documentation), and `Functional9` (the regression test cases) — see `AdminStruct.txt` in the Web9 repo for the full layout.

## Architecture

### Repository layout

- `SKIRT/core` — groups `utils`, `mpi`, and `framework`; `framework` is the library with several hundred `SimulationItem` subclasses implementing the actual physics (sources, media, spatial grids, material mixes, SEDs, instruments, probes, ...).
- `SKIRT/utils` — the `utils` library: plain (non-`SimulationItem`) math/geometry helper classes used throughout `core` (`Box`, `Direction`, `Position`, `PhotonPacket`, `Quadrics`, `PathSegmentGenerator`, ...).
- `SKIRT/main` — the `skirt` CLI entry point and command-line handling.
- `SKIRT/mpi`, `SKIRT/fitsio`, `SKIRT/tetgen`, `SKIRT/voro` — MPI support and vendored third-party libraries (CFITSIO, TetGen, Voro++).
- `SMILE/schema` — the `ItemInfo`/`NameManager`/`PropertyHandler`/`SchemaDef` metadata engine (see below).
- `SMILE/serialize` — reads/writes SMILE data sets to/from XML (i.e. `.ski` files).
- `SMILE/wizard` — the Q&A engine `MakeUp` is built on.
- `SMILE/fundamentals` — low-level shared utilities (e.g. `BooleanExpression`, used to evaluate the conditional-expression attributes below).

### The SimulationItem / SMILE metadata system

This pattern spans nearly the whole codebase and matters for adding or modifying almost any class.

- Every configurable class derives (directly or indirectly) from `SimulationItem` (itself derived from SMILE's `Item`) and declares its shape via macros from `SMILE/schema/ItemInfo.hpp`: `ITEM_CONCRETE`/`ITEM_ABSTRACT` for the class itself, `PROPERTY_*` for each discoverable property, `ATTRIBUTE_*` for metadata on the class or the most-recently-declared property (default values, min/max, visibility/relevance conditions, etc.). These macros generate the private data member, public getter, and private SMILE-only setter for each property, plus the metadata used to derive the SMILE schema.
- Conditional attributes (`ATTRIBUTE_RELEVANT_IF`, `ATTRIBUTE_DISPLAYED_IF`, `ATTRIBUTE_TYPE_ALLOWED_IF`, ...) are boolean expressions evaluated against a `NameManager` (`SMILE/schema/NameManager.*`) that accumulates named flags while walking the item tree in property-declaration order: **uppercase** names are global (persist for the whole configuration session, e.g. `"Emission"`, `"IteratePrimary"`, `"RadiationField"`), **lowercase** names are local (scoped to the current item, discarded when its subtree finishes). `ATTRIBUTE_INSERT`/`ATTRIBUTE_TYPE_INSERT` are how a class or property contributes new flags into these sets. This is how, e.g., choosing a `DustEmission` simulation mode changes which options are later offered or which defaults apply — both in the wizard and when parsing a `.ski` file, since both use the exact same schema evaluation.
- This flag mechanism is one-way (flags only ever get added, never removed) and strictly ordered by property declaration — a global flag set by one branch of the tree (e.g. while processing `sourceSystem`) is visible to everything processed afterward (e.g. `mediumSystem`, since it's declared later), but there's no way to later isolate "which branch contributed this flag." Distinguishing per-branch contributions (as needed e.g. for the per-source vs. per-medium symmetry tracking below) requires an explicit marker flag inserted at the boundary between branches, not something the mechanism gives you for free.
- **Every new `SimulationItem` subclass must be registered** in `SKIRT/core/SimulationItemRegistry.cpp`: add its header `#include` (alphabetized) and an `ItemRegistry::add<YourClass>()` call in `SimulationItemRegistry::SimulationItemRegistry()` (registration order there is the order it appears in the wizard's choice lists). A class that compiles but isn't registered is invisible to both the parser and the wizard.

### Top-level simulation structure

A `MonteCarloSimulation` (`SKIRT/core/MonteCarloSimulation.*`) is the root of the runtime object hierarchy. It holds four "system" compartments, each its own `SimulationItem` subtree: `SourceSystem` (primary emitting sources), `MediumSystem` (transfer media, the spatial grid, and radiation-field storage), `InstrumentSystem`, and `ProbeSystem`. It also owns a non-discoverable `Configuration` instance (`SKIRT/core/Configuration.*`), set up very early in the process, that other items query during their own setup for cross-cutting derived state (which simulation phases are active, dimension/symmetry requirements, wavelength grids, ...) rather than reaching across the hierarchy directly.

Setup runs depth-first: `SimulationItem::setup()` calls `setupSelfBefore()`, recurses `setup()` into children, then calls `setupSelfAfter()` — the before/after split exists because some state can only be computed correctly before vs. after children have configured themselves (and `setupSelfAfter()` is generally where cross-item queries and logging happen).

### Symmetry ("dimension") bookkeeping

Both at runtime (`Configuration::sourceDimension()`/`mediaDimension()`/`gridDimension()`/`modelDimension()`) and in the SMILE schema (`SourceDimensionN`/`MediaDimensionN`/`RequiredDimensionN`/`DimensionN` flags inserted throughout the `Geometry`/`VectorField`/`Source`/`Medium` class hierarchies), the code tracks source and media symmetry separately: primary sources aren't discretized onto the spatial grid, so the grid's required dimension is normally just the media's, expanding to the combined source+media dimension only when the radiation field is actually stored. `Geometry`/`VectorField` (and their decorators) are shared between source and medium contexts, so their schema-level dimension contribution is tagged by whether `MediumSystem`'s `InMedia` marker has already been inserted (it always follows `SourceSystem`). See `SKIRT/core/Configuration.cpp` and the `grid`/`samplingOptions` properties of `MediumSystem` for the current implementation.
