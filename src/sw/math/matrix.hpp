//
// Created by selya on 24.07.2019.
//

#ifndef SANDWICH_MATH_MATRIX_HPP
#define SANDWICH_MATH_MATRIX_HPP

#include "vector.hpp"
#include "axis_angle.hpp"

namespace sw {
namespace math {

class Matrix4 {
    // Matrix is column-major
    // 0 4 8  12
    // 1 5 9  13
    // 2 6 10 14
    // 3 7 11 15
public:
    enum Element {
        M00, M01, M02, M03,
        M10, M11, M12, M13,
        M20, M21, M22, M23,
        M30, M31, M32, M33
    };

    double m[16];

    Matrix4();
    explicit Matrix4(const double data[16]);
    Matrix4(double m00, double m01, double m02, double m03,
            double m10, double m11, double m12, double m13,
            double m20, double m21, double m22, double m23,
            double m30, double m31, double m32, double m33);
    Matrix4(const Matrix4 &other);

    Matrix4 &Identity();
    Matrix4 &Invert();
    Matrix4 &Rotate(double angle, double x, double y, double z);
    Matrix4 &Rotate(double angle, const Vector3 &axis);
    Matrix4 &Rotate(const AxisAngle &axis_angle);
    Matrix4 &Scale(double x, double y, double z);
    Matrix4 &Scale(const Vector3 &scale);
    Matrix4 &Translate(double x, double y, double z);
    Matrix4 &Translate(const Vector3 &translation);
    Matrix4 &Transpose();

    Matrix4 &Set(const double data[16]);
    Matrix4 &Set(double m00, double m01, double m02, double m03,
                 double m10, double m11, double m12, double m13,
                 double m20, double m21, double m22, double m23,
                 double m30, double m31, double m32, double m33);
    Matrix4 &Set(const Matrix4 &other);
    Matrix4 &SetOrtho(double left, double right, double bottom, double top,
            double z_near, double z_far, bool z_zero_to_one = false);
    Matrix4 &SetOrtho2D(double left, double right, double bottom, double top);
    Matrix4 &SetPerspective(double fov_y, double aspect,
            double z_near, double z_far, bool z_zero_to_one = false);

    AxisAngle GetRotation() const;
    Vector3 GetScale() const;
    Vector3 GetTranslation() const;

    double Determinant() const;

    Matrix4 &operator=(const Matrix4 &other);
    Matrix4 &operator*=(const Matrix4 &other);

    double *operator[](int index);

    Matrix4 operator*(const Matrix4 &other) const;
    Vector4 operator*(const Vector4 &other) const;

    const double *operator[](int index) const;
};

} //namespace math
} //namespace sw

#if SANDWICH_MATH_INLINE
    #include "matrix.ipp"
#endif

#endif //SANDWICH_MATH_MATRIX_HPP
