#include <openvdb/openvdb.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/LevelSetFilter.h>
#include <openvdb/tools/ParticlesToLevelSet.h>
#include "args.hxx"
#include "nlohmann/json.hpp"

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

template <typename T>
class Vec3 {
public:
    Vec3(){}
    Vec3(T _v) {
        v[0] = _v;
        v[1] = _v;
        v[2] = _v;
    }
    Vec3(T x, T y, T z) {
        v[0] = x;
        v[1] = y;
        v[2] = z;
    }

    inline T x() const { return v[0]; }
    inline T y() const { return v[1]; }
    inline T z() const { return v[2]; }

    Vec3 operator*(T nv) const {
        return Vec3(x() * nv, y() * nv, z() * nv);
    }
    Vec3 operator-(const Vec3 &nv) const {
        return Vec3(x() - nv.x(), y() - nv.y(), z() - nv.z());
    }
    Vec3 operator*(const Vec3 &nv) const {
        return Vec3(x() * nv.x(), y() * nv.y(), z() * nv.z());
    }
    Vec3 operator+(const Vec3 &nv) const {
        return Vec3(x() + nv.x(), y() + nv.y(), z() + nv.z());
    }
    Vec3 &operator+=(const Vec3 &nv) {
        v[0] += nv.x();
        v[1] += nv.y();
        v[2] += nv.z();
        return (*this);
    }
    Vec3 &operator-=(const Vec3 &nv) {
        v[0] -= nv.x();
        v[1] -= nv.y();
        v[2] -= nv.z();
        return (*this);
    }
    Vec3 operator/(const Vec3 &nv) const {
        return Vec3(x() / nv.x(),
                    y() / nv.y(),
                    z() / nv.z());
    }
    Vec3 operator-() const { return Vec3(-x(), -y(), -z()); }
    T operator[](int i) const { return v[i]; }
    T &operator[](int i) { return v[i]; }

    T v[3];
};

template <typename T>
inline Vec3<T> operator*(T nv, const Vec3<T> &v) {
    return Vec3<T>(v.x() * nv, v.y() * nv, v.z() * nv);
}

template <typename T>
inline Vec3<T> vneg(const Vec3<T> &vec) {
    return Vec3<T>(-vec.x(), -vec.y(), -vec.z());
}

template <typename T>
inline T vlength(const Vec3<T> &vec) {
    return std::sqrt(vec.x() * vec.x() + vec.y() * vec.y() + vec.z() * vec.z());
}

template <typename T>
inline Vec3<T> vnormalize(const Vec3<T> &vec) {
    Vec3<T> v = vec;
    T len = vlength(vec);
    if (std::fabs(len) > 1.0e-6f) {
        float inv_len = 1.0f / len;
        v.v[0] *= inv_len;
        v.v[1] *= inv_len;
        v.v[2] *= inv_len;
    }
    return v;
}

template <typename T>
inline Vec3<T> vcross(Vec3<T> a, Vec3<T> b) {
    Vec3<T> c;
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    return c;
}

template <typename T>
inline T vdot(Vec3<T> a, Vec3<T> b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template<typename T>
inline T vdistance(Vec3<T> a, Vec3<T> b) {
    Vec3<T> d = a - b;
    return std::sqrt(vdot(d, d));
}

template<typename T>
inline Vec3<T> vclamp(Vec3<T> a, T minVal, T maxVal) {
    return Vec3<T>(std::min(std::max(a[0], minVal), maxVal),
                   std::min(std::max(a[1], minVal), maxVal),
                   std::min(std::max(a[2], minVal), maxVal));
}

template<typename T>
inline Vec3<T> vmin(Vec3<T> a, Vec3<T> b) {
    return Vec3<T>(std::min(a[0], b[0]),
                   std::min(a[1], b[1]),
                   std::min(a[2], b[2]));
}

template<typename T>
inline Vec3<T> vmax(Vec3<T> a, Vec3<T> b) {
    return Vec3<T>(std::max(a[0], b[0]),
                   std::max(a[1], b[1]),
                   std::max(a[2], b[2]));
}

template<typename T>
inline T clamp(T a, T minVal, T maxVal) {
    return std::min(std::max(a, minVal), maxVal);
}

template<typename T>
inline Vec3<T> vmix(Vec3<T> x, Vec3<T> y, T a) {
    return x * (1 - a) + y * a;
}

template<typename T>
inline Vec3<T> vfract(Vec3<T> x) {
    return Vec3<T>(x[0] - std::floor(x[0]),
                   x[1] - std::floor(x[1]),
                   x[2] - std::floor(x[2])) ;
}

template<typename T>
inline Vec3<T> vabs(Vec3<T> x) {
    return Vec3<T>(std::abs(x[0]),
                   std::abs(x[1]),
                   std::abs(x[2]));
}

typedef Vec3<float> Vec3f;
typedef Vec3<int> Vec3i;

class Sphere {
public:
    Sphere() {}
    Sphere(Vec3f _center, float _r) {
        center = _center;
        r = _r;
    }

    void invert(Vec3f &pos, float &dr) {
        Vec3f diff = pos - center;
        float lenSq = vdot(diff, diff);
        float k = (r * r) / lenSq;
        dr *= k; // (r * r) / lenSq
        pos = diff * k + center;
    }

    float distance(const Vec3f pos) {
        return vdistance(pos, center) - r;
    }

    Vec3f center;
    float r;
};

class Plane {
public:
    Plane(Vec3f _origin, Vec3f _normal) {
        origin = _origin;
        normal = _normal;
    }

    float distance(const Vec3f pos) {
        return vdot(pos - origin, normal);
    }

    Vec3f origin;
    Vec3f normal;
};

class Sphairahedron {
public:
    Sphairahedron(Vec3f _bboxMin, Vec3f _bboxMax,
                  std::vector<Sphere> _spheres,
                  std::vector<Plane> _planes,
                  std::vector<Plane> _boundingPlanes,
                  std::vector<Plane> _dividePlanes,
                  std::vector<Sphere> _finiteSpheres,
                  std::vector<Sphere> _convexSpheres,
                  Sphere _boundingSphere) {
        bboxMin = _bboxMin;
        bboxMax = _bboxMax;
        spheres = _spheres;
        planes = _planes;
        boundingPlanes = _boundingPlanes;
        dividePlanes = _dividePlanes;
        finiteSpheres = _finiteSpheres;
        convexSpheres = _convexSpheres;
        boundingSphere = _boundingSphere;
    }

    float iisInfSphairahedron(Vec3f pos) {
        bool outside = false;
        bool inside = false;
        float minDist = std::numeric_limits<float>::max();
        float maxDist = -std::numeric_limits<float>::max();

        // bottom surface
        pos = pos - Vec3f(0, bboxMin.y(), 0);
        const Vec3f bottomNormal(0, -1, 0);
        const float d = vdot(pos, bottomNormal);
        if (d > 0.0) {
            minDist = std::min(d, minDist);
            outside = true;
        }
        pos = pos + Vec3f(0, bboxMin.y(), 0);

        std::for_each(boundingPlanes.begin(), boundingPlanes.end(),
                      [&](Plane p){
                          pos = pos - p.origin;
                          const float d = vdot(pos, p.normal);
                          if (abs(d) < 0.0001) {
                              if (d > 0.) {
                                  outside = true;
                                  minDist = std::min(d, minDist);
                              } else {
                                  inside = true;
                                  maxDist = std::max(d, maxDist);
                              }
                          }
                          pos = pos + p.origin;
                      });

        if(outside) return minDist;

        int invNum = 0;
        float dr = 1.0;
        for(int n = 0; n < MAX_IIS_ITER_COUNT; n++) {
            bool inFund = true;
            std::for_each(spheres.begin(), spheres.end(),
                          [&](Sphere s){
                              if (vdistance(pos, s.center) < s.r) {
                                  s.invert(pos, dr);
                                  invNum++;
                                  inFund = false;
                              }
                          });
            std::for_each(planes.begin(), planes.end(),
                          [&](Plane p){
                              pos = pos - p.origin;
                              float d = vdot(pos, p.normal);
                              if (d > 0.) {
                                  pos = pos - 2.f * d * p.normal;
                                  invNum++;
                                  inFund = false;
                              }
                              pos = pos + p.origin;
                          });
            if (inFund) break;
        }

        if (inside) {
            float finalDist = distInfSphairahedron(pos) / dr;
            if (finalDist < 0.0) {
                return maxDist;
            }
        }
        return distInfSphairahedron(pos) / dr;
    }

    float distInfSphairahedron(const Vec3f pos) {
        float d = -1;
        std::for_each(planes.begin(), planes.end(),
                      [&](Plane p){
                          d = std::max(p.distance(pos), d);
                      });
        std::for_each(dividePlanes.begin(), dividePlanes.end(),
                      [&](Plane p){
                          d = std::max(p.distance(pos), d);
                      });
        std::for_each(spheres.begin(), spheres.end(),
                      [&](Sphere s){
                          d = std::max(-s.distance(pos), d);
                      });
        return d;
    }

    float iisFiniteSphairahedron(Vec3f pos) {
        int invNum = 0;
        float dr = 1.0;
        for(int n = 0; n < MAX_IIS_ITER_COUNT; n++) {
            bool inFund = true;
            std::for_each(finiteSpheres.begin(), finiteSpheres.end(),
                          [&](Sphere s){
                              if (vdistance(pos, s.center) < s.r) {
                                  s.invert(pos, dr);
                                  invNum++;
                                  inFund = false;
                              }
                          });
            if (inFund) break;
        }

        return distFiniteSphairahedron(pos) / dr;
    }

    float distFiniteSphairahedron(const Vec3f pos) {
        float d = -1.;
        std::for_each(convexSpheres.begin(), convexSpheres.end(),
                      [&](Sphere s){
                          d = std::max(s.distance(pos), d);
                      });
        std::for_each(finiteSpheres.begin(), finiteSpheres.end(),
                      [&](Sphere s){
                          d = std::max(-s.distance(pos), d);
                      });
        return d;
    }

    Vec3f bboxMin;
    Vec3f bboxMax;
    std::vector<Plane> boundingPlanes;
    std::vector<Sphere> spheres;
    std::vector<Plane> planes;
    std::vector<Plane> dividePlanes;

    std::vector<Sphere> finiteSpheres;
    std::vector<Sphere> convexSpheres;
    Sphere boundingSphere;
    const int MAX_IIS_ITER_COUNT = 10000;
};

void writeObj(
    std::string fileName, std::vector<openvdb::Vec3s>& points,
    std::vector<openvdb::Vec3s>& normals,
    std::vector<openvdb::Vec4I>& quads,
    std::vector<openvdb::Vec3I>& triangles) {
    std::cout << "\nWriting \"" << fileName << ".obj\" to file\n";
    std::ofstream objfile(fileName + ".obj");

    for (auto n : points) {
        objfile << "v " << n.x() << " " << n.y() << " " << n.z() << std::endl;
    }
    for (auto n : normals) {
        objfile << "vn " << n.x() << " " << n.y() << " " << n.z() << std::endl;
    }
    objfile << std::endl;

    if (normals.size() > 0) {
        assert(points.size() == normals.size());
        for (auto q : quads) {
            objfile <<
                "f " << q.x() + 1 << "//" << q.x() + 1 <<
                " " << q.y() + 1 << "//" << q.y() + 1 <<
                " " << q.z() + 1 << "//" << q.z() + 1 <<
                " " << q.w() + 1 << "//" << q.w() + 1 << std::endl;
        }
        for (auto q : triangles) {
            objfile << "f " << q.x() + 1 << "//" << q.x() + 1 <<
                " " << q.y() + 1 << "//" << q.y() + 1 <<
                " " << q.z() + 1 << "//" << q.z() + 1<< std::endl;
        }
    } else {
        for (auto q : quads) {
            objfile <<
                "f " << q.x() + 1 <<
                " " << q.y() + 1 <<
                " " << q.z() + 1 <<
                " " << q.w() + 1 << std::endl;
        }
        for (auto q : triangles) {
            objfile << "f " << q.x() + 1 << " " << q.y() + 1 << " " << q.z() + 1 << std::endl;
        }

    }
    objfile.close();
}

void writeGrid(openvdb::GridBase::Ptr grid, std::string fileName) {
    std::cout << "\nWriting \"" << fileName << ".vdb\" to file\n";
    grid->setName("Limit set of the sphairahedron group.");
    openvdb::GridPtrVec grids;
    grids.push_back(grid);
    openvdb::io::File file(fileName + ".vdb");
    file.write(grids);
    file.close();
}

std::vector<openvdb::Vec3s> computeNormals(openvdb::FloatGrid::Ptr grid,
                                           std::vector<openvdb::Vec3s> points,
                                           bool flipNormal) {
    openvdb::Coord ikj;
    openvdb::Vec3d pos, tmpNormal, normal(0.0, -1.0, 0.0);

    //openvdb::math::Gradient<openvdb::math::GenericMap, openvdb::math::CD_2ND> Gradient;
    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
    openvdb::Coord ijk;

    std::vector<openvdb::Vec3s> normals;
    const size_t numVerts = points.size();
    normals.resize(numVerts);

    for (size_t v = 0; v < numVerts; v++) {
        pos[0] = points[v][0];
        pos[1] = points[v][1];
        pos[2] = points[v][2];

        pos = grid->worldToIndex(pos);
        ijk[0] = int(pos[0]);
        ijk[1] = int(pos[1]);
        ijk[2] = int(pos[2]);
        accessor.getValue(ijk);
        tmpNormal = openvdb::math::ISGradient<openvdb::math::CD_2ND>::result(accessor, ijk);
        double length = tmpNormal.length();
        if( length > 1.0e-7 ) {
            tmpNormal *= 1.0/length;
            normal = tmpNormal;
        }

        if (flipNormal) {
            normals[v][0] = -(float)normal[0];
            normals[v][1] = -(float)normal[1];
            normals[v][2] = -(float)normal[2];
        } else {
            normals[v][0] = (float)normal[0];
            normals[v][1] = (float)normal[1];
            normals[v][2] = (float)normal[2];
        }
    }

    return normals;
}

openvdb::FloatGrid::Ptr computeVolumeGrid(Sphairahedron sphairahedron, Vec3f sliceStep) {
    Vec3f bboxMin = sphairahedron.bboxMin;
    Vec3f bboxMax = sphairahedron.bboxMax;
    Vec3i dim (int((bboxMax.x() - bboxMin.x()) / sliceStep.x()),
               int((bboxMax.y() - bboxMin.y()) / sliceStep.y()),
               int((bboxMax.z() - bboxMin.z()) / sliceStep.z()));
    int numPoints = dim.x() * dim.y() * dim.z();

    std::cout << "dim "
              << dim.x() << " "
              << dim.y() << " "
              << dim.z() << std::endl;
    std::cout << "numPoints " << numPoints << std::endl;

    std::cout << "start" << std::endl;

//    const float voxelSize = 1.0f, halfWidth = 2.0f;
    const float voxelSize = sliceStep.x(), halfWidth = 2.0 * sliceStep.x();

     openvdb::FloatGrid::Ptr grid =
         openvdb::createLevelSet<openvdb::FloatGrid>(voxelSize, halfWidth);
     openvdb::FloatGrid::Accessor accessor = grid->getAccessor();

     std::cout << bboxMin.y() + sliceStep.y() << std::endl;
     for(int zi = 0; zi < dim.z(); zi++) {
         std::cout << zi << std::endl;
         for(int yi = 0; yi < dim.y(); yi++) {
             for(int xi = 0; xi < dim.x(); xi++) {

                 const Vec3f p (bboxMin.x() + sliceStep.x() * float(xi),
                                bboxMin.y() + sliceStep.y() * float(yi),
                                bboxMin.z() + sliceStep.z() * float(zi));
                 const float dist = sphairahedron.iisInfSphairahedron(p);
                 if(abs(dist) < 0.1 ) {
                     accessor.setValue(openvdb::Coord(xi, yi, zi), dist);
                 }
             }
         }
     }
     openvdb::tools::signedFloodFill(grid->tree());

     return grid;
}

openvdb::FloatGrid::Ptr computeVolumeGridFinite(Sphairahedron sphairahedron, Vec3f sliceStep) {
    Vec3f bboxMin = sphairahedron.boundingSphere.center - Vec3f(sphairahedron.boundingSphere.r) * 1.1;
    Vec3f bboxMax = sphairahedron.boundingSphere.center + Vec3f(sphairahedron.boundingSphere.r) * 1.1;

    Vec3i dim (int((bboxMax.x() - bboxMin.x()) / sliceStep.x()),
               int((bboxMax.y() - bboxMin.y()) / sliceStep.y()),
               int((bboxMax.z() - bboxMin.z()) / sliceStep.z()));
    int numPoints = dim.x() * dim.y() * dim.z();

    std::cout << "dim "
              << dim.x() << " "
              << dim.y() << " "
              << dim.z() << std::endl;
    std::cout << "numPoints " << numPoints << std::endl;

    std::cout << "start" << std::endl;

//    const float voxelSize = 1.0f, halfWidth = 2.0f;
    const float voxelSize = sliceStep.x(), halfWidth = 2.0 * sliceStep.x();

     openvdb::FloatGrid::Ptr grid =
         openvdb::createLevelSet<openvdb::FloatGrid>(voxelSize, halfWidth);
     openvdb::FloatGrid::Accessor accessor = grid->getAccessor();

     std::cout << bboxMin.y() + sliceStep.y() << std::endl;
     for(int zi = 0; zi < dim.z(); zi++) {
         std::cout << zi << std::endl;
         for(int yi = 0; yi < dim.y(); yi++) {
             for(int xi = 0; xi < dim.x(); xi++) {

                 const Vec3f p (bboxMin.x() + sliceStep.x() * float(xi),
                                bboxMin.y() + sliceStep.y() * float(yi),
                                bboxMin.z() + sliceStep.z() * float(zi));
                 const float dist = sphairahedron.iisFiniteSphairahedron(p);
                 if(abs(dist) < 0.1 ) {
                     accessor.setValue(openvdb::Coord(xi, yi, zi), dist);
                 }
             }
         }
     }
     openvdb::tools::signedFloodFill(grid->tree());

     return grid;
}


template<typename TreeType>
struct OffsetAndMinComp
{
    typedef typename TreeType::LeafNodeType     LeafNodeType;
    typedef typename TreeType::ValueType        ValueType;

    OffsetAndMinComp(std::vector<LeafNodeType*>& lhsNodes, const TreeType& rhsTree, ValueType offset)
        : mLhsNodes(lhsNodes.empty() ? NULL : &lhsNodes[0]), mRhsTree(&rhsTree), mOffset(offset)
    {
    }

    void operator()(const tbb::blocked_range<size_t>& range) const
    {
        typedef typename LeafNodeType::ValueOnIter Iterator;

        openvdb::tree::ValueAccessor<const TreeType> rhsAcc(*mRhsTree);
        const ValueType offset = mOffset;

        for (size_t n = range.begin(), N = range.end(); n < N; ++n) {

            LeafNodeType& lhsNode = *mLhsNodes[n];
            const LeafNodeType * rhsNodePt = rhsAcc.probeConstLeaf(lhsNode.origin());
            if (!rhsNodePt) continue;

            for (Iterator it = lhsNode.beginValueOn(); it; ++it) {
                ValueType& val = const_cast<ValueType&>(it.getValue());
                val = std::min(val, offset + rhsNodePt->getValue(it.pos()));
            }
        }
    }

private:
    LeafNodeType    *       * const mLhsNodes;
    TreeType          const * const mRhsTree;
    ValueType                 const mOffset;
}; // struct OffsetAndMinComp


template<typename GridType, typename InterrupterType>
inline void
normalizeLevelSet(GridType& grid, const int halfWidthInVoxels, InterrupterType* interrupt = NULL)
{
    openvdb::tools::LevelSetFilter<GridType, GridType, InterrupterType> filter(grid, interrupt);
    filter.setSpatialScheme(openvdb::math::FIRST_BIAS);
    filter.setNormCount(halfWidthInVoxels);
    filter.normalize();
    filter.prune();
}

// https://github.com/dreamworksanimation/openvdb/blob/a7a1abde7c955d5a83acac2c82341497d73576c2/openvdb/tools/TopologyToLevelSet.h
template<typename GridType, typename InterrupterType>
inline void
smoothLevelSet(GridType& grid, int iterations, int halfBandWidthInVoxels, InterrupterType* interrupt = NULL)
{
    typedef typename GridType::ValueType        ValueType;
    typedef typename GridType::TreeType         TreeType;
    typedef typename TreeType::LeafNodeType     LeafNodeType;

    GridType filterGrid(grid);

    openvdb::tools::LevelSetFilter<GridType, GridType, InterrupterType> filter(filterGrid, interrupt);
    filter.setSpatialScheme(openvdb::math::FIRST_BIAS);

    for (int n = 0; n < iterations; ++n) {
        if (interrupt && interrupt->wasInterrupted()) break;
        filter.mean(1);
    }

    std::vector<LeafNodeType*> nodes;
    grid.tree().getNodes(nodes);

    const ValueType offset = ValueType(double(0.01) * grid.transform().voxelSize()[0]);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, nodes.size()),
        OffsetAndMinComp<TreeType>(nodes, filterGrid.tree(), -offset));

    (void) halfBandWidthInVoxels;
    // Clean up any damanage that was done by the min operation
    normalizeLevelSet(grid, halfBandWidthInVoxels, interrupt);
}

void makeMesh(Sphairahedron sphairahedron, Vec3f sliceStep, int smoothIterations, bool isFinite) {
    bool flipNormal = true;
    bool enableAdaptiveMeshing = false;
    float adaptivity = 0.0f;
    float isovalue = 0.f;

    openvdb::FloatGrid::Ptr grid;
    if (isFinite) {
        grid = computeVolumeGridFinite(sphairahedron, sliceStep);
    } else {
        grid = computeVolumeGrid(sphairahedron, sliceStep);
    }

    grid->tree().print(std::cout, 4);
    writeGrid(grid, "IISVolume");

    std::vector<openvdb::Vec3s> points;
    std::vector<openvdb::Vec4I> quads;
    std::vector<openvdb::Vec3I> triangles;

    if (smoothIterations > 0) {
        openvdb::util::NullInterrupter interrupt;
        const float halfWidth = 2.0 * sliceStep.x();
        smoothLevelSet(*grid, smoothIterations, halfWidth, &interrupt);
    }

    if (enableAdaptiveMeshing) {
        openvdb::tools::volumeToMesh(*grid, points, triangles, quads, isovalue, adaptivity);
    } else {
        // Uniformly mesh any scalar grid that has a continuous isosurface.
        openvdb::tools::volumeToMesh(*grid, points, quads, isovalue);
    }

    std::vector<openvdb::Vec3s> normals = computeNormals(grid, points, flipNormal);

    writeGrid(grid, "IISVolumeProcessed");
    writeObj("IIS", points, normals, quads, triangles);
}

Sphairahedron createSphairahedronFromJson(nlohmann::json jsonObj) {
    Vec3f bboxMin(jsonObj["bboxMin"][0], jsonObj["bboxMin"][1], jsonObj["bboxMin"][2]);
    Vec3f bboxMax(jsonObj["bboxMax"][0], jsonObj["bboxMax"][1], jsonObj["bboxMax"][2]);

    std::vector<Sphere> spheres;
    for (auto data: jsonObj["prismSpheres"]) {
        spheres.push_back(Sphere(Vec3f(data["center"][0],
                                       data["center"][1],
                                       data["center"][2]),
                                 data["r"].get<float>()));
    }

    std::vector<Sphere> finiteSpheres;
    for (auto data: jsonObj["finiteSpheres"]) {
        finiteSpheres.push_back(Sphere(Vec3f(data["center"][0],
                                             data["center"][1],
                                             data["center"][2]),
                                       data["r"].get<float>()));
    }

    Sphere boundingSphere(Vec3f(jsonObj["boundingSphere"]["center"][0],
                                jsonObj["boundingSphere"]["center"][1],
                                jsonObj["boundingSphere"]["center"][2]),
                          jsonObj["boundingSphere"]["r"].get<float>());

    std::vector<Sphere> convexSpheres;
    for (auto data: jsonObj["convexSpheres"]) {
        convexSpheres.push_back(Sphere(Vec3f(data["center"][0],
                                             data["center"][1],
                                             data["center"][2]),
                                       data["r"].get<float>()));
    }

    std::vector<Plane> planes;
    for (auto data: jsonObj["prismPlanes"]) {
        planes.push_back(Plane(Vec3f(data["p1"][0],
                                     data["p1"][1],
                                     data["p1"][2]),
                               Vec3f(data["normal"][0],
                                     data["normal"][1],
                                     data["normal"][2])));
    }

    std::vector<Plane> boundingPlanes;
    for (auto data: jsonObj["boundingPlanes"]) {
        Vec3f origin(data["p1"][0],
                     data["p1"][1],
                     data["p1"][2]);
        Vec3f normal(data["normal"][0],
                     data["normal"][1],
                     data["normal"][2]);
        boundingPlanes.push_back(Plane(origin, normal));
    }

    std::vector<Plane> dividePlanes;
    for (auto data: jsonObj["dividePlanes"]) {
        dividePlanes.push_back(Plane(Vec3f(data["p1"][0],
                                           data["p1"][1],
                                           data["p1"][2]),
                                     Vec3f(data["normal"][0],
                                           data["normal"][1],
                                           data["normal"][2])));
    }

    std::cout << "bbox min (" << bboxMin.x() << ", "
              << bboxMin.y() << ", "
              << bboxMin.z()  << ") "<< std::endl;
    std::cout << "bbox max (" << bboxMax.x() << ", "
              << bboxMax.y() << ", "
              << bboxMax.z()  << ") "<< std::endl;
    std::cout << "number of prism spheres " << spheres.size() << std::endl;
    std::cout << "number of prism planes " << planes.size() << std::endl;

    return Sphairahedron(bboxMin, bboxMax,
                         spheres, planes,
                         boundingPlanes,
                         dividePlanes,
                         finiteSpheres,
                         convexSpheres,
                         boundingSphere);
}

// Hello World for OpenVDB
// http://www.openvdb.org/documentation/doxygen/codeExamples.html
int main(int argc, char** argv) {
    args::ArgumentParser parser("Generate mesh.");
    args::ValueFlag<float> sliceStep(parser, "step", "slice step", {'s', "sliceStep"}, 0.01f);
    args::ValueFlag<std::string> inputJson(parser, "input", "Input json file", {'i', "input"}, "scene.json");
    args::ValueFlag<int> smoothIterations(parser, "smoothIterations", "Iterations of smoothing", {"smooth"}, 0);
    args::Flag isFinite(parser, "isFinite", "Generate finite limit set", {'f', "finite"});
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    } catch (args::ValidationError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    std::string inputJsonFileName;
    inputJsonFileName = args::get(inputJson);

    std::ifstream ifs(inputJsonFileName);
    nlohmann::json jsonObj;
    if (!ifs) {
        std::cout << "Can't open " << inputJsonFileName << std::endl;
        return 1;
    }
    ifs >> jsonObj;
    ifs.close();

    std::cout << "slice step: " << args::get(sliceStep) << std::endl;
    std::cout << "smooth iterations: " << args::get(smoothIterations) << std::endl;
    Sphairahedron s = createSphairahedronFromJson(jsonObj);
    makeMesh(s, Vec3f(args::get(sliceStep)), args::get(smoothIterations), isFinite);
}
