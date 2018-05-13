#pragma once
#include <cmath>

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
