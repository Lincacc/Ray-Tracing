#pragma once 

#include "core.h"

#include <cmath>
#include <iostream>

namespace LinRT {
    class vec3 {
    public:
        vec3 () : x(0), y(0), z(0) {}
        vec3 (global_precision value) : x(value), y(value), z(value) {}
        vec3 (global_precision x, global_precision y, global_precision z) : x(x), y(y), z(z) {}

        vec3 operator-() const { return vec3(-x, -y, -z); }
        global_precision operator[](int i) const { return data[i]; }
        global_precision& operator[](int i) { return data[i]; }

        vec3 operator*=(global_precision scalar) {
            x *= scalar;
            y *= scalar;
            z *= scalar;
            return *this;
        }

        vec3 operator*=(const vec3& other) {
            x *= other.x;
            y *= other.y;
            z *= other.z;
            return *this;
        }

        vec3 operator+=(const vec3& other) {
            x += other.x;
            y += other.y;
            z += other.z;
            return *this;
        }

        vec3 operator/=(global_precision scalar) {
            x /= scalar;
            y /= scalar;
            z /= scalar;
            return *this;
        }

        global_precision length() const {
            return std::sqrt(x * x + y * y + z * z);
        }

        global_precision length_squared() const {
            return x * x + y * y + z * z;
        }

        bool is_zero() const {
            return x == 0.0f && y == 0.0f && z == 0.0f;
        }

    public:
        union {
            global_precision data[3];
            struct { global_precision x, y, z; };
            struct { global_precision r, g, b; };
        };
    };

    inline vec3 operator+(const vec3& a, const vec3& b) {
        return vec3(a.x + b.x, a.y + b.y, a.z + b.z);
    }

    inline vec3 operator-(const vec3& a, const vec3& b) {
        return vec3(a.x - b.x, a.y - b.y, a.z - b.z);
    }

    inline vec3 operator*(const vec3& a, global_precision scalar) {
        return vec3(a.x * scalar, a.y * scalar, a.z * scalar);
    }

    inline vec3 operator*(global_precision scalar, const vec3& a) {
        return a * scalar;
    }

    inline vec3 operator*(const vec3& a, const vec3& b) {
        return vec3(a.x * b.x, a.y * b.y, a.z * b.z);
    }

    inline vec3 operator/(const vec3& a, global_precision scalar) {
        return vec3(a.x / scalar, a.y / scalar, a.z / scalar);
    }

    inline global_precision dot(const vec3& a, const vec3& b) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    inline vec3 cross(const vec3& a, const vec3& b) {
        return vec3(
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x
        );
    }

    inline vec3 normalize(const vec3& v) {
        return v / v.length();
    }

    inline std::ostream& operator<<(std::ostream& os, const vec3& v) {
        os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
        return os;
    }
}