#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <GLFW/glfw3.h>

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

const float GAMMA_COEFF = 1. / 2.2;
inline Vec3f gammaCorrect(Vec3f &rgb) {
    rgb[0] = std::min(std::pow(rgb[0], GAMMA_COEFF), 1.f);
    rgb[1] = std::min(std::pow(rgb[1], GAMMA_COEFF), 1.f);
    rgb[2] = std::min(std::pow(rgb[2], GAMMA_COEFF), 1.f);
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
    return hsv2rgb(((invNum) * 0.01), 1., 1.);
//    return hsv2rgb((-0.13 + pos.y()), 1., 1.);
}

int main(int argc, char** argv) {
    GLFWwindow* window;

    /* Initialize the library */
    if (!glfwInit())
        return -1;

    /* Create a windowed mode window and its OpenGL context */
    ::glfwWindowHint( GLFW_VISIBLE, 0 );
    window = glfwCreateWindow(640, 480, "Hello World", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    /* Make the window's context current */
    glfwMakeContextCurrent(window);

    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window))
    {
        /* Render here */
        glClear(GL_COLOR_BUFFER_BIT);

        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        /* Poll for and process events */
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}
