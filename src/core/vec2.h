#pragma once 

#include "core.h"

#include <cmath>
#include <iostream>

namespace LinRT {
    class vec2 {
    public:
        vec2 () : x(0), y(0) {}
        vec2 (global_precision value) : x(value), y(value) {}
        vec2 (global_precision x, global_precision y) : x(x), y(y) {}

        vec2 operator-() const { return vec2(-x, -y); }
        global_precision operator[](int i) const { return data[i]; }
        global_precision& operator[](int i) { return data[i]; }

        vec2 operator*=(global_precision scalar) {
            x *= scalar;
            y *= scalar;
            return *this;
        }

        vec2 operator+=(const vec2& other) {
            x += other.x;
            y += other.y;
            return *this;
        }

        vec2 operator/=(global_precision scalar) {
            x /= scalar;
            y /= scalar;
            return *this;
        }

    public:
        union {
            global_precision data[2];
            struct { global_precision x, y; };
            struct { global_precision r, g; };
        };
    };

    inline vec2 operator+(const vec2& a, const vec2& b) {
        return vec2(a.x + b.x, a.y + b.y);
    }

    inline vec2 operator-(const vec2& a, const vec2& b) {
        return vec2(a.x - b.x, a.y - b.y);
    }

    inline vec2 operator*(const vec2& a, global_precision scalar) {
        return vec2(a.x * scalar, a.y * scalar);
    }

    inline vec2 operator*(global_precision scalar, const vec2& a) {
        return a * scalar;
    }

    inline vec2 operator/(const vec2& a, global_precision scalar) {
        return vec2(a.x / scalar, a.y / scalar);
    }

    inline std::ostream& operator<<(std::ostream& os, const vec2& v) {
        os << "(" << v.x << ", " << v.y << ")";
        return os;
    }
}