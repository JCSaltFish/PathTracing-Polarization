# Polarizer

**Polarizer** is a desktop path tracing renderer written in **C++17**. It combines:

- A GPU **compute-shader path tracer** (progressive accumulation)
- A real-time **previewer** for interactive scene editing/selection
- An experimental **polarization** pipeline (Stokes-like vector + Fresnel s/p terms)
- A custom binary scene format (`*.pls`) backed by a lightweight **DB** module
- Cross-backend **GFX** abstraction (**Vulkan** + **OpenGL**) and an ImGui-based **GUI** layer

---

## Build

Polarizer is a **cross-platform** CMake project targeting **Windows** and **Linux**.

### Requirements

- **CMake ≥ 3.15**
- A **C++17** compiler
  - Windows: Visual Studio (MSVC) recommended
  - Linux: GCC or Clang
- **Vulkan headers/libs** (recommended)
  - The project can use a Vulkan SDK installed system-wide, or a vendored SDK path (see below).
- GPU drivers supporting:
  - **Vulkan** compute (recommended path; the app config defaults to Vulkan)
  - and/or **OpenGL** (there is an OpenGL backend as well)

### Third-party dependencies

Most dependencies are expected under `third_party/` and are wired by `CMakeLists.txt`, including:
- GLFW
- glad
- ImGui
- stb
- tinyfiledialogs
- nlohmann/json (header-only)

### Vulkan SDK path

`CMakeLists.txt` uses `VULKAN_SDK_PATH` to locate Vulkan headers/libraries:
- If `VULKAN_SDK_PATH` is **not** set, it defaults to: `third_party/vulkan`.

If you have an installed Vulkan SDK, point CMake at it, e.g.:

```bash
cmake -S . -B build -DVULKAN_SDK_PATH="C:/VulkanSDK/1.xx.x.x"   # Windows example
# or
cmake -S . -B build -DVULKAN_SDK_PATH="$VULKAN_SDK"            # Linux example
```

### Build commands

#### Windows (Visual Studio)

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
```

Artifacts are placed under:
- `build/bin` (executables)
- `build/lib` (libraries)

#### Linux (Makefiles / Ninja)

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

### Notes / troubleshooting

- If you see missing Vulkan headers/libs at configure time, set `VULKAN_SDK_PATH` (or ensure the vendored `third_party/vulkan` is present).
- If Vulkan is not available on your system, you may still be able to run with the OpenGL backend (the repo contains both backends), but Vulkan is the default in app initialization.

---

## Run

Polarizer accepts an optional scene file (`*.pls`) as the first argument:

```bash
./Polarizer path/to/scene.pls
```

Drag & drop:
- Drop a `.pls` file to open it
- Drop `.obj` files to load models into the current scene

---

## UI overview

The main application is `PathTracerApp` (`inc/app/PathTracerApp.h`, `src/app/PathTracerApp.cpp`). The main window is composed of:

- **Menu bar** (`UiMenuBar`)
- **Tool bar** (`UiToolBar`) with render controls (Start/Pause/Stop/Restart) and display-mode switching
- **Main viewport** (`UiMainViewport`): shows preview output or path tracer output depending on display mode
- **Right panel** (`UiRightPanel`): scene/camera/model/mesh/material properties
- **Status bar** (`UiStatusBar`): render state, samples, timing, triangle count
- Settings/About/Save prompts (`UiSettingsWindow`, `UiAboutWindow`, `UiSaveDialog`)

Keyboard/mouse highlights:
- **Right mouse drag** (in preview mode) rotates the camera; **W/A/S/D** moves while in camera-control mode.
- **Shift+1 / Shift+2** switches display mode (Preview / Path tracer output).
- **F5/F6/F7/F8** map to render Start/Pause/Stop/Restart (see `onKeyboardAction`).

Strings are localized via JSON in `resources/strings/en_US.json` and `resources/strings/zh_CN.json`.

---

## Scene files: `*.pls`

A `*.pls` file is a binary dump of the custom DB:
- Magic: `P L S`
- Version: derived from app version
- Root object: a `PtScene`

The stored scene graph is made of typed objects:
- `PtScene`: resolution, trace depth, camera settings, list of models
- `PtModel`: name, file path (OBJ), transforms, list of meshes
- `PtMesh`: name, links to model + material
- `PtMaterial`: material model + parameters + texture references used by the renderer/UI

Opening a scene resets the DB, loads it from file, validates the root `PtScene`, then populates the previewer and UI.

---

## Polarization theory

Polarization is implemented inside the path tracing compute shader (`resources/shaders/pathTracer.comp`), using:

- Per-bounce stored data `PolarInfo` (normal, outgoing dir, Fresnel reflection coefficients `rs/rp`, and a Stokes-like vector)
- Fresnel s/p reflection coefficient computation (`getRsRp`)
- Rotation of the Stokes-like vector between interaction frames via a 3×3 rotation matrix (`stokesRotation(phi)`)
- Backward accumulation across the bounce history (`calcPolarResult(rayDir, depth)`), rotating into each bounce frame and applying Fresnel-derived updates

This is intended to model how polarization state changes through specular-like interactions and basis changes along a path.

---

## Module map

This repo is organized as several mostly-independent modules. Here’s a high-level map of what each does:

### App module (`src/app`, `inc/app`)
- Owns the main loop, window/view composition, input routing, and orchestration of previewer/path tracer/post-process.
- Scene loading/saving, model import, image export, and user workflows live here.

### DB module (`src/db`, `inc/db`)
- Generic typed object store (`DB`) with:
  - type registry (`DbTypeRegistry`)
  - binary serialization/deserialization (`DbSerializer`-based)
  - transaction log + undo/redo support
- Used as the authoritative backing store for `.pls` scenes and editor operations.

### Data types / scene graph (`src/app/data`, `inc/app/data`)
- App-specific DB-registered types: `PtScene`, `PtModel`, `PtMesh`, `PtMaterial`.
- These define exactly what is persisted in `.pls`.

### GFX module (`src/gfx`, `inc/gfx`)
- Rendering abstraction layer with multiple backends:
  - **OpenGL** backend
  - **Vulkan** backend
- Provides cross-backend concepts like renderer, buffers, images, pipelines, descriptor binding, and dispatch/copy operations.
- Used by:
  - the previewer (raster pipeline)
  - the path tracer (compute pipeline)
  - post-processing

### GUI module (`src/gui`, `inc/gui`, plus `inc/app/ui`)
- Windowing/input layer + view/event system on top of GLFW/ImGui-style widgets.
- UI widgets/views emit `GuiEvent`s which `PathTracerApp` consumes to update DB state and rendering state.
- Localization uses `GuiText` with JSON string tables.

### Resources module (`resources/*`, `inc/res/*`)
- Shaders under `resources/shaders/*`
- UI language strings under `resources/strings/*`
- Shaders are also embedded into a generated header (`inc/res/ShaderStrings.hpp`) via CMake (`cmake/GenerateShaderStrings.cmake`).

### Utils module (`src/utils`, `inc/utils`)
- Shared helpers (math types, logging, image IO via stb, OBJ parsing/mesh loading, timers, etc.)

---

## License

Source files include Apache-2.0 headers. See repository contents for the full license text if provided.

## Acknowledgements

- stb (image loading/writing)
- GLFW
- ImGui
- Vulkan / OpenGL ecosystem
- nlohmann/json
- tinyfiledialogs
