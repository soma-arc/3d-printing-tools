# Repository Guidelines

## Project Structure & Module Organization
- Root contains two C++ console apps: `GenMesh/` (OpenVDB/TBB mesh generator) and `GenTex/` (OpenGL-based texture generator).
- Source: `GenMesh/src/`, `GenTex/src/` (headers alongside sources).
- Build output: `GenMesh/build/<Config>/bin/`, `GenTex/` project target in current dir unless otherwise configured.
- Assets/Shaders: `GenTex/src/*.vert|*.frag` and single-file deps under `src/`.

## Build, Test, and Development Commands
- Generate makefiles (GNU Make):
  - `cd GenMesh && premake5 gmake`
  - `cd GenTex && premake5 gmake`
- Build Release (example):
  - `make -C GenMesh/build config=Release`
  - `make -C GenTex config=Release`
- Run locally (examples):
  - `GenMesh/build/Release/bin/GenMesh --help`
  - `./GenTex` (or `./genTex_debug` in Debug)
- Requirements: Linux toolchain, Premake5, and system libs: `openvdb`, `tbb`, `GL`, `GLEW`, `glfw`. Configure runtime linker paths if needed (e.g., `LD_LIBRARY_PATH`).

## Coding Style & Naming Conventions
- Language: C++17.
- Indentation: 4 spaces; brace style: same-line (K&R).
- Naming: Classes `PascalCase` (e.g., `Sphairahedron`), functions/methods `camelCase` (e.g., `getUniLocations`), files `snake_case.cpp/h` when adding new.
- Headers may be `.h`/`.hpp`; keep lightweight and include-local where possible.
- Prefer standard library + single-header deps already used in `src/`.

## Testing Guidelines
- No formal tests yet. If adding tests:
  - Framework: Catch2 or GoogleTest.
  - Location: `tests/` mirroring module structure; name `*_test.cpp`.
  - Integrate via Premake and add a `test` make target.

## Commit & Pull Request Guidelines
- Commits: short, imperative, scoped messages (e.g., `Add sphairahedron class`, `Fix shader uniforms`). Group related changes; one logical unit per commit.
- PRs: clear description of intent, key changes, how to build/run, and any screenshots of visual output. Link issues and note platform/toolchain tested.
- CI not configured; ensure local Release build succeeds for both projects before requesting review.

## Security & Configuration Tips
- Ensure matching versions of `openvdb` and `tbb`; verify linker paths (`ldconfig` or rpath) to avoid runtime errors.
- For `GenTex`, confirm GPU drivers and OpenGL support; run on X11/Wayland session with valid context.

