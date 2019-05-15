#pragma once

#include <cmath>

class Vec3 {
public:
    double e[3] = {0.0,0.0,0.0};
    Vec3() { }
    Vec3(double x) { e[0] = x; e[1] = x; e[2] = x; }
    Vec3(double x, double y) { e[0] = x; e[1] = y; e[2] = 0.; }
    Vec3(double x, double y, double z) { e[0] = x; e[1] = y; e[2] = z; }
    inline double x() const { return e[0]; }
    inline double y() const { return e[1]; }
    inline double z() const { return e[2]; }

    inline Vec3 operator-() const { return Vec3(-e[0], -e[1], -e[2]); }
    inline double operator[](int i) const { return e[i]; }
    inline double& operator[](int i) { return e[i]; }

    inline Vec3 operator+(const Vec3 &v) const; //Element Addition
    inline Vec3 operator-(const Vec3 &v) const; //Element Subtraction
    inline Vec3 operator*(const Vec3 &v) const; //Element Multiplication
    inline Vec3 operator/(const Vec3 &v) const; //Element Division

    inline Vec3 operator+(double t) const; //float Addition
    inline Vec3 operator-(double t) const; //float Subtraction
    inline Vec3 operator*(double t) const; //float Multiplication
    inline Vec3 operator/(double t) const; //float Division
    
    inline Vec3& operator+=(const Vec3 &v);
    inline Vec3& operator-=(const Vec3 &v);
    inline Vec3& operator*=(const Vec3 &v);
    inline Vec3& operator/=(const Vec3 &v);
    inline Vec3& operator*=(const double t);
    inline Vec3& operator/=(const double t);

    inline double Length() const { 
        return std::sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
    }
    inline double LengthSquared() const {
        return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
    }
    inline void Make_Unit();
};

inline void Vec3::Make_Unit() {
    double inv = 1.0 / Length();
    e[0] *= inv; e[1] *= inv; e[2] *= inv;
}

inline Vec3 Vec3::operator+(const Vec3 &v) const {
    return Vec3(e[0] + v.e[0], e[1] + v.e[1], e[2] + v.e[2]);
}

inline Vec3 Vec3::operator-(const Vec3 &v) const {
    return Vec3(e[0] - v.e[0], e[1] - v.e[1], e[2] - v.e[2]);
}

inline Vec3 Vec3::operator*(const Vec3 &v) const {
    return Vec3(e[0] * v.e[0], e[1] * v.e[1], e[2] * v.e[2]);
}

inline Vec3 Vec3::operator/(const Vec3 &v) const {
    return Vec3(e[0] / v.e[0], e[1] / v.e[1], e[2] / v.e[2]);
}

inline Vec3 Vec3::operator+(double t) const {
    return Vec3(e[0] + t, e[1] + t, e[2] + t);
}

inline Vec3 Vec3::operator-(double t) const {
    return Vec3(e[0] - t, e[1] - t, e[2] - t);
}

inline Vec3 Vec3::operator*(double t) const {
    return Vec3(e[0] * t, e[1] * t, e[2] * t);
}

inline Vec3 Vec3::operator/(double t) const {
    double inv = 1.0 / t;
    return Vec3(e[0] * inv, e[1] * inv, e[2] * inv);
}

inline Vec3& Vec3::operator+=(const Vec3 &v) {
    e[0] += v.e[0]; e[1] += v.e[1]; e[2] += v.e[2];
    return *this;
}

inline Vec3& Vec3::operator-=(const Vec3 &v) {
    e[0] -= v.e[0]; e[1] -= v.e[1]; e[2] -= v.e[2];
    return *this;
}

inline Vec3& Vec3::operator*=(const Vec3 &v) {
    e[0] *= v.e[0]; e[1] *= v.e[1]; e[2] *= v.e[2];
    return *this;
}

inline Vec3& Vec3::operator/=(const Vec3 &v) {
    e[0] /= v.e[0]; e[1] /= v.e[1]; e[2] /= v.e[2];
    return *this;
}

inline Vec3& Vec3::operator*=(const double t) {
    e[0] *= t; e[1] *= t; e[2] *= t;
    return *this;
}

inline Vec3& Vec3::operator/=(const double t) {
    double inv = 1.0 / t;
    e[0] *= inv; e[1] *= inv; e[2] *= inv;
    return *this;
}

inline Vec3 operator*(double t, const Vec3 &v) {
    return v * t;
}

inline double Dot(const Vec3 &v1, const Vec3 &v2) {
    return v1.e[0]*v2.e[0] + v1.e[1]*v2.e[1] + v1.e[2]*v2.e[2];
}

inline Vec3 Cross(const Vec3 &v1, const Vec3 &v2) {
    double v1x = v1.e[0], v1y = v1.e[1], v1z = v1.e[2];
    double v2x = v2.e[0], v2y = v2.e[1], v2z = v2.e[2];
    return Vec3((v1y * v2z) - (v1z * v2y), 
                (v1z * v2x) - (v1x * v2z),
                (v1x * v2y) - (v1y * v2x));
}

inline Vec3 Unit_Vector(Vec3 v) {
    return v / v.Length();
}