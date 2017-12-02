#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "args.hxx"

#define TINYOBJLOADER_IMPLEMENTATION
#include "./tiny_obj_loader.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

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

inline Vec3f hsv2rgb(float h, float s, float v){
    Vec3f c = Vec3f(h, s, v);
    const Vec3f hue(h, h, h);
    const Vec3f K = Vec3f(1.0, 2.0 / 3.0, 1.0 / 3.0);
    const Vec3f Kx_ONE = Vec3f(1, 1, 1);
    const Vec3f Kw = Vec3f(3, 3, 3);
    Vec3f p = vabs(vfract(hue + K) * 6.0 - Kw);
    return c.z() * vmix(Kx_ONE, vclamp(p - Kx_ONE, 0.0f, 1.0f), c.y());
}

class Plane {
public:
    Plane(Vec3f _origin, Vec3f _normal) {
        origin = _origin;
        normal = _normal;
    }
    Vec3f origin;
    Vec3f normal;
};

inline float distSphere(const Vec3f pos, const Sphere s) {
    return vdistance(pos, s.center) - s.r;
}

inline float distPlane(const Vec3f pos, const Plane p) {
    return vdot(pos - p.origin, p.normal);
}

inline float distInfSphairahedron(const Vec3f pos,
                                  const std::vector<Sphere> spheres,
                                  const std::vector<Plane> planes,
                                  const Plane divPlane) {
    float d = -1;
    d = std::max(distPlane(pos, planes[0]), d);
    d = std::max(distPlane(pos, planes[1]), d);
    d = std::max(distPlane(pos, planes[2]), d);
    d = std::max(distPlane(pos, divPlane), d);

    d = std::max(-distSphere(pos, spheres[0]), d);
    d = std::max(-distSphere(pos, spheres[1]), d);
    d = std::max(-distSphere(pos, spheres[2]), d);
    return d;
}

const int MAX_ITER_COUNT = 1000;
int IIS(Vec3f pos,
        std::vector<Sphere> spheres,
        std::vector<Plane> planes,
        const Plane divPlane) {
    int invNum = 0;
    float dr = 1.0;
    for(int n = 0; n < MAX_ITER_COUNT; n++) {
        bool inFund = true;
        std::for_each(spheres.begin(), spheres.end(),
                      [&](Sphere s){
                          if (vdistance(pos, s.center) < s.r) {
                              //std::cout << pos.x() << std::endl;
                              s.invert(pos, dr);
                              //std::cout << pos.x() << std::endl;
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

    return invNum;
    // if(distInfSphairahedron(pos,
    //                         spheres,
    //                         planes,
    //                         divPlane) / dr <= 0.01) {
    //     return invNum;
    // }
    // return -1;
}

const Vec3f offset(0.5, 0, std::sqrt(3.) * 0.5);
const std::vector<Sphere> spheres = {
    Sphere(Vec3f(0.25450630091522675, 0, 0)  + offset,
           0.7454936990847733),
    Sphere(Vec3f(-0.3116633792528053,
                 0.6053931133878944,
                 0.5398168077244666) + offset,
           0.37667324149438935),
    Sphere(Vec3f(-0.12608782704164367,
                 1.2165336555165491,
                 -0.21839052265208383) + offset,
           0.7478243459167127)
};

const std::vector<Plane> planes = {
    Plane(Vec3f(0,
                5,
                0.5773502691896258) + offset,
          Vec3f(0.5,
                0,
                0.8660254037844388)),
    Plane(Vec3f(0,
                3,
                -0.5773502691896258) + offset,
          Vec3f(0.5,
                0,
                -0.8660254037844388)),
    Plane(Vec3f(-0.5, 0, 1) + offset,
          Vec3f(-1, 0, 0))
};

const Plane divPlane(Vec3f(0.9999999999999948,
                           -1.3100631690576847e-14,
                           7.91033905045424e-15) + offset,
                     Vec3f(0.4969732028017572,
                           0.8183202716219945,
                           0.2887378893554992));

inline Vec3f computeColor (Vec3f pos) {
    pos += Vec3f(2, -2, 2);
    float sink = std::sin(M_PI / 3.);
    float cosk = std::cos(M_PI / 3.);
    pos = Vec3f(pos.x() * cosk - pos.z() * sink,
                pos.y(),
                pos.x() * sink + pos.z() * cosk);
    int invNum = IIS(pos, spheres, planes, divPlane);
    if (invNum == -1) {
        return Vec3f(0, 0, 0);
    }
    return hsv2rgb((float(invNum) * 0.01), 1., 1.);
//    return hsv2rgb((-0.13 + pos.y()), 1., 1.);
}

int main(int argc, char** argv) {
    args::ArgumentParser parser("This is a test program.", "This goes after the options.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::Positional<std::string> inputObjArg(parser, "input",
                                              "input .obj file");
    args::ValueFlag<int> texSizeArg(parser, "integer",
                                         "The size of the generated texture. (default: 2048)", {'s'});
    args::ValueFlag<std::string> outTexArg(parser, "string",
                                        "The name of the generated texture. (default: texture.png)", {'o'});

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

    const std::string inputObjName = args::get(inputObjArg);
    const int texSize = texSizeArg ? args::get(texSizeArg) : 2048;
    const std::string outTexName = outTexArg ? args::get(outTexArg) : "texture.png";

    const float uvStep = 1.0f / texSize / 2.f;
    printf("texture size ... %d x %d\n", texSize, texSize);
    printf("uv step ... %f\n", uvStep);

    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;

    std::string err;

    std::cout << "Load " << inputObjName << std::endl;
    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, inputObjName.c_str());

    if (!err.empty()) { // `err` may contain warning message.
        std::cerr << err << std::endl;
    }

    if (!ret) {
        std::cerr << "Failed to load " << inputObjName << std::endl;
        exit(1);
    }

    if (attrib.texcoords.size() == 0) {
        std::cerr << "This object has no texture coordinates."  << std::endl;
        exit(1);
    }

    printf("# of vertices  = %d\n", (int) (attrib.vertices.size()) / 3);
    printf("# of normals   = %d\n", (int) (attrib.normals.size()) / 3);
    printf("# of texcoords = %d\n", (int) (attrib.texcoords.size()) / 2);
    printf("# of materials = %d\n", (int) materials.size());
    printf("# of shapes    = %d\n", (int) shapes.size());
    printf("\n");

    unsigned char *textureData = new unsigned char[texSize * texSize * 3];
    std::fill(textureData, textureData + texSize * texSize * 3, 0);

    for (size_t s = 0; s < shapes.size(); s++) {
        printf("Shape %d\n", (int) shapes.size());
        printf("# of faces %d\n", (int) shapes[s].mesh.num_face_vertices.size());
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
            // The number of vertexes in the face. 3, 4 or 5?
            int fv = shapes[s].mesh.num_face_vertices[f];
            if(fv > 3) {
                printf("quad!!");
                continue;
            }

            std::vector<Vec3f> faceVert;
            std::vector<Vec3f> faceUv;
            Vec3f uvMin(2, 0, 2);
            Vec3f uvMax(-1, 0, -1);
            // Loop over vertices in the face.
            for (int v = 0; v < fv; v++) {
                // access to vertex
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                tinyobj::real_t vx = attrib.vertices[3*idx.vertex_index+0];
                tinyobj::real_t vy = attrib.vertices[3*idx.vertex_index+1];
                tinyobj::real_t vz = attrib.vertices[3*idx.vertex_index+2];
                tinyobj::real_t nx = attrib.normals[3*idx.normal_index+0];
                tinyobj::real_t ny = attrib.normals[3*idx.normal_index+1];
                tinyobj::real_t nz = attrib.normals[3*idx.normal_index+2];
                tinyobj::real_t tx = attrib.texcoords[2*idx.texcoord_index+0];
                tinyobj::real_t ty = attrib.texcoords[2*idx.texcoord_index+1];
                faceVert.push_back(Vec3f(vx, vy, vz));
                faceUv.push_back(Vec3f(tx, 0, ty));
                uvMin[0] = std::min(tx, uvMin[0]);
                uvMin[2] = std::min(ty, uvMin[2]);
                uvMax[0] = std::max(tx, uvMax[0]);
                uvMax[2] = std::max(ty, uvMax[2]);
            }
            index_offset += fv;

            // printf("(%f, %f), (%f, %f), (%f, %f)\n",
            //        faceUv[0].x(), faceUv[0].z(),
            //        faceUv[1].x(), faceUv[1].z(),
            //        faceUv[2].x(), faceUv[2].z());
            // printf("BBox (%f, %f) ~ (%f, %f)\n",
            //        uvMin.x(), uvMin.z(),
            //        uvMax.x(), uvMax.z());
            float area2 = vlength(vcross(faceUv[1] - faceUv[0],
                                         faceUv[2] - faceUv[0]));
            // printf("Squared area ... %f\n", area2);
            for(float u = uvMin[0] - uvStep * 2; u < uvMax[0] + uvStep * 2; u += uvStep) {
                for (float v = uvMin[2] - uvStep * 2; v < uvMax[2] + uvStep * 2; v += uvStep) {
                    Vec3f uv(u, 0, v);
                    Vec3f e1 = faceUv[2] - faceUv[1];
                    float pu = vlength(vcross(e1, uv - faceUv[1])) / area2;

                    Vec3f e2 = faceUv[0] - faceUv[2];
                    float pv = vlength(vcross(e2, uv - faceUv[2])) / area2;
                    float pw = 1.f - pu - pv;
                    if(pu > 1.01f || pv > 1.01f || pw > 1.01f ||
                       pu < -0.01f || pv < -0.01f || pw < -0.01f) {
                        // Outside of the face
                        // printf("barycentric coordinates (%f, %f, %f)\n", pu, pv, pw);
                    } else
                    {
                        const int x = clamp(float(round(u * (texSize -1))),
                                            0.f, float(texSize - 1));
                        const int y = clamp(float(round((1.f - v) * (texSize - 1))),
                                            0.f, float(texSize - 1));
                        const int index = y * texSize + x;
                        Vec3f coord = faceVert[0] * pu + faceVert[1] * pv + faceVert[2] * pw;
                        Vec3f rgb = computeColor(coord);
                        textureData[index * 3 + 0] =
                            (unsigned char) std::max(0.0f,
                                                     std::min(rgb[0] * 255.f, 255.0f));
                        textureData[index * 3 + 1] =
                            (unsigned char) std::max(0.0f,
                                                     std::min(rgb[1] * 255.f, 255.0f));
                        textureData[index * 3 + 2] =
                            (unsigned char) std::max(0.0f,
                                                     std::min(rgb[2] * 255.f, 255.0f));
                        //printf("barycentric coordinates (%f, %f, %f)\n", pu, pv, pw);
                    }
                }
            }
            // printf("\n");
        }
    }

    int n = stbi_write_png(outTexName.c_str(), texSize, texSize, 3, textureData, texSize * 3);
    delete[] textureData;

    if (n == 0) {
        fprintf(stderr, "Error to save PNG image:\n");
    } else {
        printf("Saved image to [ %s ]\n", outTexName.c_str());
    }

    return 0;
}
