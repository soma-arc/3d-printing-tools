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

using json = nlohmann::json;

template <typename T>
class Vec3 {
public:
    Vec3(){}
    Vec3(T v) {
        v[0] = v;
        v[1] = v;
        v[2] = v;
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

    Vec3f center;
    float r;
};

class Plane {
public:
    Plane(Vec3f _origin, Vec3f _normal) {
        origin = _origin;
        normal = _normal;
    }

    Vec3f origin;
    Vec3f normal;
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

openvdb::FloatGrid::Ptr computeVolumeGrid(Vec3f bboxMin, Vec3f bboxMax, Vec3f sliceStep) {
    Vec3i dim (int((bboxMax.x() - bboxMin.x()) / sliceStep.x()),
               int((bboxMax.y() - bboxMin.y()) / sliceStep.y()),
               int((bboxMax.z() - bboxMin.z()) / sliceStep.z()));
    int numPoints = dim.x() * dim.y() * dim.z();

    std::vector<int> invNums;
    invNums.resize(numPoints);
    std::vector<float> distances;
    distances.resize(numPoints);
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

     int i = 0;
     for(int zi = 0; zi < dim.z(); zi++) {
         std::cout << zi << std::endl;
         for(int yi = 0; yi < dim.y(); yi++) {
             for(int xi = 0; xi < dim.x(); xi++) {
                 if (yi == 0) {
                     distances[i] = 0.;
                     i++;
                     continue;
                 }
                 Vec3f p (bboxMin.x() + sliceStep.x() * float(xi),
                          bboxMin.y() + sliceStep.y() * float(yi),
                          bboxMin.z() + sliceStep.z() * float(zi));
                 (void) p;
//                int invNum = IIS(p, prismSpheres, prismPlanes, divPlane);
                 float dd = 1.f;//distIIS(p, prismSpheres, prismPlanes, divPlane);
//                distances[i] = invNum < 0 ? 0 : 10;
                 if(abs(dd) < 0.01 ) {
                     accessor.setValue(openvdb::Coord(xi, yi, zi), dd);
                 }
                 i++;
             }
         }
     }
     openvdb::tools::signedFloodFill(grid->tree());

     return grid;
}

void makeMesh(Vec3f bboxMin, Vec3f bboxMax, Vec3f sliceStep) {
    bool flipNormal = false;
    int smoothIterations = 0;
    bool enableAdaptiveMeshing = false;
    float adaptivity = 0.0f;
    float isovalue = 0.f;

    openvdb::FloatGrid::Ptr grid = computeVolumeGrid(bboxMin, bboxMax, sliceStep);

    grid->tree().print(std::cout, 4);
    writeGrid(grid, "IISVolume");

    std::vector<openvdb::Vec3s> points;
    std::vector<openvdb::Vec4I> quads;
    std::vector<openvdb::Vec3I> triangles;

    if (smoothIterations > 0) {
        openvdb::util::NullInterrupter interrupt;
        //smoothLevelSet(*grid, smoothIterations, halfWidth, &interrupt);
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

// Hello World for OpenVDB
// http://www.openvdb.org/documentation/doxygen/codeExamples.html
int main(int argc, char** argv) {
    args::ArgumentParser parser("This is a test program.");
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
}
