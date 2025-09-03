# 3d-printing-tools
3d-printing tools for 3d fractals

## GenMesh
- Purpose: generate meshes from fractal scenes (OpenVDB/TBB).
- Default output: `.glb` (binary glTF). Use flags to also export `.obj`/`.vdb`.

Quick start:

```
cd GenMesh
premake5 gmake
make -C build config=release
./build/Release/bin/GenMesh -i scene.json [-o basename] [--obj] [--vdb]
```

Common options:
- `-s, --sliceStep <float>`: voxel/slice step (default: 0.01)
- `--smooth <int>`: smoothing iterations (default: 0)
- `-f, --finite`: finite limit set

## GenTex
