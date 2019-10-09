//
// Created by selya on 05.10.2019.
//

#ifndef SANDWICH_MATH_AXIS_ANGLE_HPP
#define SANDWICH_MATH_AXIS_ANGLE_HPP

#include "vector.hpp"

namespace sw {
namespace math {

class AxisAngle {
public:
    Vector3 axis;
    double angle;

    AxisAngle();
    AxisAngle(double angle, double x, double y, double z);
    AxisAngle(double angle, const Vector3 &axis);
    AxisAngle(const AxisAngle &other);

    AxisAngle &Rotate(double angle);

    AxisAngle &Set(double angle, double x, double y, double z);
    AxisAngle &Set(double angle, const Vector3 &axis);
    AxisAngle &Set(const AxisAngle &other);

    AxisAngle &operator=(const AxisAngle &other);
};

} //namespace math
} //namespace sw

#if SANDWICH_MATH_INLINE
    #include "axis_angle.ipp"
#endif

#endif //SANDWICH_MATH_AXIS_ANGLE_HPP
