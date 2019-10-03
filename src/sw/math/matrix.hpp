//
// Created by selya on 24.07.2019.
//

#ifndef SANDWICH_MATH_MATRIX_HPP
#define SANDWICH_MATH_MATRIX_HPP

#include "vector.hpp"

namespace sw {
namespace math {

// Matrix is column-major
// 0 4 8  12
// 1 5 9  13
// 2 6 10 14
// 3 7 11 15

class Matrix4 {
private:
    double data[16];

public:
    Matrix4();

    Matrix4(const Matrix4 &other);

    Matrix4(const double data[16]);

    //
    // Non const methods
    //

    Matrix4 &Identity();

    Matrix4 &Inverse();

    Matrix4 &Transpose();

    Matrix4 &Set(const Matrix4 &other);

    Matrix4 &Set(const double data[16]);

    Matrix4 &Scale(double x, double y, double z);

    Matrix4 &Scale(const Vector3 &scale);

    Matrix4 &Rotate(double ang, double x, double y, double z);

    Matrix4 &Rotate(double ang, const Vector3 &axis);

    Matrix4 &Translate(double x, double y, double z);

    Matrix4 &Translate(const Vector3 &translate);

    Matrix4 &SetOrtho2D(double left, double right, double bottom, double top);

    Matrix4 &operator=(const Matrix4 &other);

    double *operator[](int index);

    //
    // Const methods
    //

    Matrix4 operator*(const Matrix4 &other) const;

    Vector4 operator*(const Vector4 &v) const;

    const double *operator[](int index) const;

    /*Vector4 Unproject(double x, double y, double z) const;
    Vector4 Unproject(const Vector3 &v) const;*/
    Vector3 GetTranslation() const;

    double GetRotationZ() const;
};

} //namespace math
} //namespace sw

#endif //SANDWICH_MATH_MATRIX_HPP
