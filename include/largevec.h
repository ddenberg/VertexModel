#pragma once

#include <cmath>

class LargeVec {
private:
    double *data;
    int n;
public:
    LargeVec() {}
    LargeVec(const int _n) { n = _n; data = new double[_n]; }
    inline void Init(const int _n) { n = _n; data = new double[_n]; }
    ~LargeVec() { delete[] data; }

    inline double operator[](int i) const { return data[i]; }
    inline double& operator[](int i) { return data[i]; }

    inline LargeVec operator+(const LargeVec &v) const; //Element Addition
    inline LargeVec operator-(const LargeVec &v) const; //Element Subtraction
    inline LargeVec operator*(const LargeVec &v) const; //Element Multiplication
    inline LargeVec operator/(const LargeVec &v) const; //Element Division

    inline LargeVec operator+(double t) const; //float Addition
    inline LargeVec operator-(double t) const; //float Subtraction
    inline LargeVec operator*(double t) const; //float Multiplication
    inline LargeVec operator/(double t) const; //float Division
    
    inline LargeVec& operator+=(const LargeVec &v);
    inline LargeVec& operator-=(const LargeVec &v);
    inline LargeVec& operator*=(const LargeVec &v);
    inline LargeVec& operator/=(const LargeVec &v);
    inline LargeVec& operator*=(const double t);
    inline LargeVec& operator/=(const double t);

    inline LargeVec& operator=(const LargeVec &v);

    inline void setZero();
    inline double norm();
    inline double normSquared();
    inline double dot(const LargeVec &v);
};

inline void LargeVec::setZero() {
    for (int i = 0; i < n; i++) {
        data[i] = 0;
    }
}

inline double LargeVec::norm() {
    return std::sqrt(normSquared());
}

inline double LargeVec::normSquared() {
    double value = 0;
    for (int i = 0; i < n; i++) {
        value += data[i] * data[i];
    }
    return value;
}

inline double LargeVec::dot(const LargeVec &v) {
    double value = 0;
    for (int i = 0; i < n; i++) {
        value += data[i] * v[i];
    }
    return value;
}

inline LargeVec LargeVec::operator+(const LargeVec &v) const {
    LargeVec temp = LargeVec(n);
    for (int i = 0; i < n; i++) {
        temp[i] = data[i] + v[i];
    }
    return temp;
}

inline LargeVec LargeVec::operator-(const LargeVec &v) const {
    LargeVec temp = LargeVec(n);
    for (int i = 0; i < n; i++) {
        temp[i] = data[i] - v[i];
    }
    return temp;
}

inline LargeVec LargeVec::operator*(const LargeVec &v) const {
    LargeVec temp = LargeVec(n);
    for (int i = 0; i < n; i++) {
        temp[i] = data[i] * v[i];
    }
    return temp;
}

inline LargeVec LargeVec::operator/(const LargeVec &v) const {
    LargeVec temp = LargeVec(n);
    for (int i = 0; i < n; i++) {
        temp[i] = data[i] / v[i];
    }
    return temp;
}

inline LargeVec LargeVec::operator+(double t) const {
    LargeVec temp = LargeVec(n);
    for (int i = 0; i < n; i++) {
        temp[i] = data[i] + t;
    }
    return temp;
}

inline LargeVec LargeVec::operator-(double t) const {
    LargeVec temp = LargeVec(n);
    for (int i = 0; i < n; i++) {
        temp[i] = data[i] - t;
    }
    return temp;
}

inline LargeVec LargeVec::operator*(double t) const {
    LargeVec temp = LargeVec(n);
    for (int i = 0; i < n; i++) {
        temp[i] = data[i] * t;
    }
    return temp;
}

inline LargeVec LargeVec::operator/(double t) const {
    LargeVec temp = LargeVec(n);
    double inv = 1.0 / t;
    for (int i = 0; i < n; i++) {
        temp[i] = data[i] * inv;
    }
    return temp;
}

inline LargeVec& LargeVec::operator+=(const LargeVec &v) {
    for (int i = 0; i < n; i++) {
        data[i] += v[i];
    }
    return *this;
}

inline LargeVec& LargeVec::operator-=(const LargeVec &v) {
    for (int i = 0; i < n; i++) {
        data[i] -= v[i];
    }
    return *this;
}

inline LargeVec& LargeVec::operator*=(const LargeVec &v) {
    for (int i = 0; i < n; i++) {
        data[i] *= v[i];
    }
    return *this;
}

inline LargeVec& LargeVec::operator/=(const LargeVec &v) {
    for (int i = 0; i < n; i++) {
        data[i] /= v[i];
    }
    return *this;
}

inline LargeVec& LargeVec::operator*=(const double t) {
    for (int i = 0; i < n; i++) {
        data[i] *= t;
    }
    return *this;
}

inline LargeVec& LargeVec::operator/=(const double t) {
    double inv = 1.0 / t;
    for (int i = 0; i < n; i++) {
        data[i] *= inv;
    }
    return *this;
}

inline LargeVec& LargeVec::operator=(const LargeVec &v) {
    for (int i = 0; i < n; i++) {
        data[i] = v[i];
    }
    return *this;
}

inline LargeVec operator*(double t, const LargeVec &v) {
    return v * t;
}