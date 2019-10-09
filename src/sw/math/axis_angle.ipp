//
// Created by selya on 06.10.2019.
//

#include "axis_angle.hpp"

#include <cmath>

#if SANDWICH_MATH_INLINE
    #define INLINE inline
#else
    #define INLINE
#endif

namespace sw {
namespace math {

INLINE AxisAngle::AxisAngle() : axis(0.0, 0.0, 1.0), angle(0.0) {}
INLINE AxisAngle::AxisAngle(double angle, double x, double y, double z)
    : angle(angle), axis(x, y, z) {}
INLINE AxisAngle::AxisAngle(double angle, const Vector3 &axis)
    : angle(angle), axis(axis) {}
INLINE AxisAngle::AxisAngle(const AxisAngle &other)
    : angle(other.angle), axis(other.axis) {}

INLINE AxisAngle& AxisAngle::Rotate(double angle_) {
    const double pi = 3.14159265358979323846;
    angle += angle_;
    angle = std::fmod((angle < 0.0 ? pi + pi + std::fmod(angle, pi + pi) : angle), pi + pi);
    return *this;
}

INLINE AxisAngle& AxisAngle::Set(double angle_, double x, double y, double z) {
    angle = angle_;
    axis.Set(x, y, z);
    return *this;
}

INLINE AxisAngle& AxisAngle::Set(double angle_, const Vector3 &axis_) {
    angle = angle_;
    axis = axis_;
    return *this;
}

INLINE AxisAngle& AxisAngle::Set(const AxisAngle &other) {
    angle = other.angle;
    axis = other.axis;
    return *this;
}

INLINE AxisAngle& AxisAngle::operator=(const AxisAngle &other) {
    angle = other.angle;
    axis = other.axis;
    return *this;
}

} //namespace math
} //namespace sw