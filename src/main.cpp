#include <iostream>
#include <vector>

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
        return Vec3(x() / nv.x(), y() / nv.y(), z() / nv.z());
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

typedef Vec3<float> Vec3f;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Needs input.obj\n" << std::endl;
        return 0;
    }

    const int texSize = 1024;
    const float uvStep = 1.0f / texSize;
    printf("texture size ... %d x %d\n", texSize, texSize);
    printf("uv step ... %f\n", uvStep);

    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;

    std::string err;

    std::cout << "Load " << argv[1] << std::endl;
    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, argv[1]);

    if (!err.empty()) { // `err` may contain warning message.
        std::cerr << err << std::endl;
    }

    if (!ret) {
        std::cerr << "Failed to load " << argv[1] << std::endl;
        exit(1);
    }

    printf("# of vertices  = %d\n", (int) (attrib.vertices.size()) / 3);
    printf("# of normals   = %d\n", (int) (attrib.normals.size()) / 3);
    printf("# of texcoords = %d\n", (int) (attrib.texcoords.size()) / 2);
    printf("# of materials = %d\n", (int) materials.size());
    printf("# of shapes    = %d\n", (int) shapes.size());
    printf("\n");

    unsigned char *textureData = new unsigned char[texSize * texSize * 3];
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

            printf("(%f, %f), (%f, %f), (%f, %f)\n",
                   faceUv[0].x(), faceUv[0].z(),
                   faceUv[1].x(), faceUv[1].z(),
                   faceUv[2].x(), faceUv[2].z());
            printf("BBox (%f, %f) ~ (%f, %f)\n",
                   uvMin.x(), uvMin.z(),
                   uvMax.x(), uvMax.z());
            float area2 = vlength(vcross(faceUv[1] - faceUv[0],
                                         faceUv[2] - faceUv[0]));
            printf("Squared area ... %f\n", area2);
            for(float u = uvMin[0]; u < uvMax[0]; u += uvStep) {
                for (float v = uvMin[2]; v < uvMax[2]; v += uvStep) {
                    Vec3f uv(u, 0, v);
                    Vec3f e1 = faceUv[2] - faceUv[1];
                    float pu = vlength(vcross(e1, uv - faceUv[1])) / area2;

                    Vec3f e2 = faceUv[0] - faceUv[2];
                    float pv = vlength(vcross(e2, uv - faceUv[2])) / area2;
                    float pw = 1.f - pu - pv;
                    if(pu > 1.f || pv > 1.f || pw > 1.f ||
                       pu < 0.f || pv < 0.f || pw < 0.f) {
                        // Outside of the face
                        // printf("barycentric coordinates (%f, %f, %f)\n", pu, pv, pw);
                    } else {
                        const int x = u * texSize;
                        const int y = (1.f - v) * texSize;
                        const int index = y * texSize + x;
                        Vec3f coord = faceVert[0] * pu + faceVert[1] * pv + faceVert[2] * pw;
                        textureData[index * 3 + 0] = (unsigned char)std::max(0.0,
                                                                             std::min((coord.y() + 1.) / 2.f * 255.f,
                                                                                      255.0));
                        textureData[index * 3 + 1] = (unsigned char) 0.f;
                        textureData[index * 3 + 2] = (unsigned char) 0.f;
                        //printf("barycentric coordinates (%f, %f, %f)\n", pu, pv, pw);
                    }
                }
            }
            printf("\n");
        }
    }

    int n = stbi_write_png("texture.png", texSize, texSize, 3, textureData, texSize * 3);
    delete[] textureData;

    if (n == 0) {
        fprintf(stderr, "Error to save PNG image:\n");
    } else {
        printf("Saved image to [ %s ]\n", "texture.png");
    }

    return 0;
}
