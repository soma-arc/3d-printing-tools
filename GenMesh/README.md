# GenMesh
Generate mesh

## Build

```
premake5 --openvdb-path=/usr/local/OpenVDB gmake
```

## Usage

Default output is a binary glTF file (`.glb`). Use flags to additionally export `.obj` or `.vdb`.

Examples:

```
# Generate GLB from scene.json (default basename = scene)
./build/Release/bin/GenMesh -i scene.json

# Specify basename and also export OBJ
./build/Release/bin/GenMesh -i scene.json -o output_name --obj

# Export VDB together with GLB
./build/Release/bin/GenMesh -i scene.json --vdb

# Common options
#  -s, --sliceStep <float>   voxel/slice step (default: 0.01)
#  --smooth <int>            smoothing iterations (default: 0)
#  -f, --finite              generate finite limit set
```

## Third party libraries
- [OpenVDB](http://www.openvdb.org/)
- [args](https://github.com/Taywee/args)
  A simple header-only C++ argument parser library. Supposed to be flexible and powerful, and attempts to be compatible with the functionality of the Python standard argparse library (though not necessarily the API).
 - [tinygltf](https://github.com/syoyo/tinygltf)
   Header-only glTF 2.0 loader/serializer used for `.glb` export.
