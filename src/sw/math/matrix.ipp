//
// Created by selya on 24.07.2019.
//

#include "matrix.hpp"

#include "utils.hpp"
#include "simd.hpp"

#include <cmath>
#include <cstring>

#if SANDWICH_MATH_INLINE
    #define INLINE inline
#else
    #define INLINE
#endif

namespace sw {
namespace math {

INLINE Matrix4::Matrix4() : m {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
} {}

INLINE Matrix4::Matrix4(const double *data) {
    std::memcpy(m, data, sizeof(m));
}

INLINE Matrix4::Matrix4(
        double m00, double m01, double m02, double m03,
        double m10, double m11, double m12, double m13,
        double m20, double m21, double m22, double m23,
        double m30, double m31, double m32, double m33) {
    m[M00] = m00;
    m[M01] = m01;
    m[M02] = m02;
    m[M03] = m03;

    m[M10] = m10;
    m[M11] = m11;
    m[M12] = m12;
    m[M13] = m13;

    m[M20] = m20;
    m[M21] = m21;
    m[M22] = m22;
    m[M23] = m23;

    m[M30] = m30;
    m[M31] = m31;
    m[M32] = m32;
    m[M33] = m33;
}

INLINE Matrix4::Matrix4(const Matrix4 &other) {
    std::memcpy(m, other.m, sizeof(m));
}

INLINE Matrix4 &Matrix4::Identity() {
    m[M00] = 1;
    m[M01] = 0;
    m[M02] = 0;
    m[M03] = 0;

    m[M10] = 0;
    m[M11] = 1;
    m[M12] = 0;
    m[M13] = 0;

    m[M20] = 0;
    m[M21] = 0;
    m[M22] = 1;
    m[M23] = 0;

    m[M30] = 0;
    m[M31] = 0;
    m[M32] = 0;
    m[M33] = 1;

    return *this;
}

INLINE Matrix4 &Matrix4::Invert() {
    double a = m[M00] * m[M11] - m[M01] * m[M10];
    double b = m[M00] * m[M12] - m[M02] * m[M10];
    double c = m[M00] * m[M13] - m[M03] * m[M10];
    double d = m[M01] * m[M12] - m[M02] * m[M11];
    double e = m[M01] * m[M13] - m[M03] * m[M11];
    double f = m[M02] * m[M13] - m[M03] * m[M12];
    double g = m[M20] * m[M31] - m[M21] * m[M30];
    double h = m[M20] * m[M32] - m[M22] * m[M30];
    double i = m[M20] * m[M33] - m[M23] * m[M30];
    double j = m[M21] * m[M32] - m[M22] * m[M31];
    double k = m[M21] * m[M33] - m[M23] * m[M31];
    double l = m[M22] * m[M33] - m[M23] * m[M32];

    double det_inv = 1.0 / (a * l - b * k + c * j + d * i - e * h + f * g);
    double nm[16];

    nm[M00] = ( m[M11] * l - m[M12] * k + m[M13] * j) * det_inv;
    nm[M01] = (-m[M01] * l + m[M02] * k - m[M03] * j) * det_inv;
    nm[M02] = ( m[M31] * f - m[M32] * e + m[M33] * d) * det_inv;
    nm[M03] = (-m[M21] * f + m[M22] * e - m[M23] * d) * det_inv;
    nm[M10] = (-m[M10] * l + m[M12] * i - m[M13] * h) * det_inv;
    nm[M11] = ( m[M00] * l - m[M02] * i + m[M03] * h) * det_inv;
    nm[M12] = (-m[M30] * f + m[M32] * c - m[M33] * b) * det_inv;
    nm[M13] = ( m[M20] * f - m[M22] * c + m[M23] * b) * det_inv;
    nm[M20] = ( m[M10] * k - m[M11] * i + m[M13] * g) * det_inv;
    nm[M21] = (-m[M00] * k + m[M01] * i - m[M03] * g) * det_inv;
    nm[M22] = ( m[M30] * e - m[M31] * c + m[M33] * a) * det_inv;
    nm[M23] = (-m[M20] * e + m[M21] * c - m[M23] * a) * det_inv;
    nm[M30] = (-m[M10] * j + m[M11] * h - m[M12] * g) * det_inv;
    nm[M31] = ( m[M00] * j - m[M01] * h + m[M02] * g) * det_inv;
    nm[M32] = (-m[M30] * d + m[M31] * b - m[M32] * a) * det_inv;
    nm[M33] = ( m[M20] * d - m[M21] * b + m[M22] * a) * det_inv;

    std::memcpy(m, nm, sizeof(m));

    return *this;
}

INLINE Matrix4 &Matrix4::Rotate(double angle, double x, double y, double z) {
    double sin = std::sin(angle);
    double cos = std::cos(angle);

    double C = 1.0f - cos;

    double xx = x * x;
    double xy = x * y;
    double xz = x * z;
    double yy = y * y;
    double yz = y * z;
    double zz = z * z;

    double rm00 = xx * C + cos;
    double rm01 = xy * C + z * sin;
    double rm02 = xz * C - y * sin;
    double rm10 = xy * C - z * sin;
    double rm11 = yy * C + cos;
    double rm12 = yz * C + x * sin;
    double rm20 = xz * C + y * sin;
    double rm21 = yz * C - x * sin;
    double rm22 = zz * C + cos;
    double nm00 = m[M00] * rm00 + m[M10] * rm01 + m[M20] * rm02;
    double nm01 = m[M01] * rm00 + m[M11] * rm01 + m[M21] * rm02;
    double nm02 = m[M02] * rm00 + m[M12] * rm01 + m[M22] * rm02;
    double nm03 = m[M03] * rm00 + m[M13] * rm01 + m[M23] * rm02;
    double nm10 = m[M00] * rm10 + m[M10] * rm11 + m[M20] * rm12;
    double nm11 = m[M01] * rm10 + m[M11] * rm11 + m[M21] * rm12;
    double nm12 = m[M02] * rm10 + m[M12] * rm11 + m[M22] * rm12;
    double nm13 = m[M03] * rm10 + m[M13] * rm11 + m[M23] * rm12;
    m[M20] = m[M00] * rm20 + m[M10] * rm21 + m[M20] * rm22;
    m[M21] = m[M01] * rm20 + m[M11] * rm21 + m[M21] * rm22;
    m[M22] = m[M02] * rm20 + m[M12] * rm21 + m[M22] * rm22;
    m[M23] = m[M03] * rm20 + m[M13] * rm21 + m[M23] * rm22;
    m[M00] = nm00;
    m[M01] = nm01;
    m[M02] = nm02;
    m[M03] = nm03;
    m[M10] = nm10;
    m[M11] = nm11;
    m[M12] = nm12;
    m[M13] = nm13;

    return *this;
}

INLINE Matrix4 &Matrix4::Rotate(double angle, const Vector3 &axis) {
    return Rotate(angle, axis.x, axis.y, axis.z);
}

INLINE Matrix4 &Matrix4::Rotate(const AxisAngle &axis_angle) {
    return Rotate(axis_angle.angle, axis_angle.axis.x, axis_angle.axis.y, axis_angle.axis.z);
}

INLINE Matrix4 &Matrix4::Scale(double x, double y, double z) {
    m[M00] *= x;
    m[M01] *= x;
    m[M02] *= x;
    m[M03] *= x;

    m[M10] *= y;
    m[M11] *= y;
    m[M12] *= y;
    m[M13] *= y;

    m[M20] *= z;
    m[M21] *= z;
    m[M22] *= z;
    m[M23] *= z;

    return *this;
}

INLINE Matrix4 &Matrix4::Scale(const Vector3 &scale) {
    return Scale(scale.x, scale.y, scale.z);
}

INLINE Matrix4 &Matrix4::Translate(double x, double y, double z) {
    m[M30] = m[M00] * x + m[M10] * y + m[M20] * z + m[M30];
    m[M31] = m[M01] * x + m[M11] * y + m[M21] * z + m[M31];
    m[M32] = m[M02] * x + m[M12] * y + m[M22] * z + m[M32];
    m[M33] = m[M03] * x + m[M13] * y + m[M23] * z + m[M33];
    return *this;
}

INLINE Matrix4 &Matrix4::Translate(const Vector3 &translation) {
    return Translate(translation.x, translation.y, translation.z);
}

INLINE Matrix4 &Matrix4::Transpose() {
    double temp;
    temp = m[M01]; m[M01] = m[M10]; m[M10] = temp;
    temp = m[M02]; m[M02] = m[M20]; m[M20] = temp;
    temp = m[M03]; m[M03] = m[M30]; m[M30] = temp;
    temp = m[M12]; m[M12] = m[M21]; m[M21] = temp;
    temp = m[M13]; m[M13] = m[M31]; m[M31] = temp;
    temp = m[M23]; m[M23] = m[M32]; m[M32] = temp;
    return *this;
}

INLINE Matrix4 &Matrix4::Set(const double *data) {
    std::memcpy(m, data, sizeof(m));
    return *this;
}

INLINE Matrix4 &Matrix4::Set(
        double m00, double m01, double m02, double m03,
        double m10, double m11, double m12, double m13,
        double m20, double m21, double m22, double m23,
        double m30, double m31, double m32, double m33) {
    m[M00] = m00;
    m[M01] = m01;
    m[M02] = m02;
    m[M03] = m03;

    m[M10] = m10;
    m[M11] = m11;
    m[M12] = m12;
    m[M13] = m13;

    m[M20] = m20;
    m[M21] = m21;
    m[M22] = m22;
    m[M23] = m23;

    m[M30] = m30;
    m[M31] = m31;
    m[M32] = m32;
    m[M33] = m33;

    return *this;
}

INLINE Matrix4 &Matrix4::Set(const Matrix4 &other) {
    std::memcpy(m, other.m, sizeof(m));
    return *this;
}

INLINE Matrix4 &Matrix4::SetOrtho(double left, double right, double bottom, double top,
        double z_near, double z_far, bool z_zero_to_one) {
    m[M00] = 2.0 / (right - left);
    m[M01] = 0.0;
    m[M02] = 0.0;
    m[M03] = 0.0;

    m[M10] = 0.0;
    m[M11] = 2.0 / (top - bottom);
    m[M12] = 0.0;
    m[M13] = 0.0;

    m[M20] = 0.0;
    m[M21] = 0.0;
    m[M22] = (z_zero_to_one ? 1.0 : 2.0) / (z_near - z_far);
    m[M23] = 0.0;

    m[M30] = (right + left) / (left - right);
    m[M31] = (top + bottom) / (bottom - top);
    m[M32] = (z_zero_to_one ? z_near : (z_far + z_near)) / (z_near - z_far);
    m[M33] = 1.0;

    return *this;
}

INLINE Matrix4 &Matrix4::SetOrtho2D(double left, double right, double bottom, double top) {
    m[M00] = 2.0 / (right - left);
    m[M01] = 0.0;
    m[M02] = 0.0;
    m[M03] = 0.0;

    m[M10] = 0.0;
    m[M11] = 2.0 / (top - bottom);
    m[M12] = 0.0;
    m[M13] = 0.0;

    m[M20] = 0.0;
    m[M21] = 0.0;
    m[M22] = -1.0;
    m[M23] = 0.0;

    m[M30] = (right + left) / (left - right);
    m[M31] = (top + bottom) / (bottom - top);
    m[M32] = 0.0;
    m[M33] = 1.0;

    return *this;
}

INLINE Matrix4 &Matrix4::SetPerspective(double fov_y, double aspect,
        double z_near, double z_far, bool z_zero_to_one) {
    double h = std::tan(fov_y * 0.5);

    m[M00] = 1.0 / (h * aspect);
    m[M01] = 0.0;
    m[M02] = 0.0;
    m[M03] = 0.0;

    m[M10] = 0.0;
    m[M11] = 1.0 / h;
    m[M12] = 0.0;
    m[M13] = 0.0;

    m[M20] = 0.0;
    m[M21] = 0.0;

    bool far_inf = z_far > 0.0 && std::isinf(z_far);
    bool near_inf = z_near > 0.0 && std::isinf(z_near);
    if (far_inf) {
        // See: "Infinite Projection Matrix" (http://www.terathon.com/gdc07_lengyel.pdf)
        double e = 1E-6;
        m[M22] = e - 1.0;
        m[M32] = (e - (z_zero_to_one ? 1.0 : 2.0)) * z_near;
    } else if (near_inf) {
        double e = 1E-6;
        m[M22] = (z_zero_to_one ? 0.0 : 1.0) - e;
        m[M32] = ((z_zero_to_one ? 1.0 : 2.0) - e) * z_far;
    } else {
        m[M22] = (z_zero_to_one ? z_far : z_far + z_near) / (z_near - z_far);
        m[M32] = (z_zero_to_one ? z_far : z_far + z_far) * z_near / (z_near - z_far);
    }

    m[M23] = -1.0;
    m[M30] = 0.0;
    m[M31] = 0.0;
    m[M33] = 0.0;

    return *this;
}

INLINE AxisAngle Matrix4::GetRotation() const {
    double epsilon = 1E-4;
    if (std::abs(m[M10] - m[M01]) < epsilon &&
        std::abs(m[M20] - m[M02]) < epsilon &&
        std::abs(m[M21] - m[M12]) < epsilon) {
        AxisAngle a(3.14159265358979323846, 0.0, 0.0, 0.0);
        double xx = (m[M00] + 1.0) * 0.5;
        double yy = (m[M11] + 1.0) * 0.5;
        double zz = (m[M22] + 1.0) * 0.5;
        double xy = (m[M10] + m[M01]) * 0.25;
        double xz = (m[M20] + m[M02]) * 0.25;
        double yz = (m[M21] + m[M12]) * 0.25;
        if (xx > yy && xx > zz) {
            a.axis.x = std::sqrt(xx);
            a.axis.y = xy / a.axis.x;
            a.axis.z = xz / a.axis.x;
        } else if (yy > zz) {
            a.axis.y = std::sqrt(yy);
            a.axis.x = xy / a.axis.y;
            a.axis.z = yz / a.axis.y;
        } else {
            a.axis.z = std::sqrt(zz);
            a.axis.x = xz / a.axis.z;
            a.axis.y = yz / a.axis.z;
        }
        return a;
    }
    double s = 1.0 / std::sqrt(
            (m[M12] - m[M21]) * (m[M12] - m[M21]) +
            (m[M20] - m[M02]) * (m[M20] - m[M02]) +
            (m[M01] - m[M10]) * (m[M01] - m[M10])
    );
    double cos = (m[M00] + m[M11] + m[M22] - 1.0) * 0.5;
    return AxisAngle(
            std::acos(cos < -1.0 ? -1.0 : (cos > 1.0 ? 1.0 : cos)),
            (m[M12] - m[M21]) * s,
            (m[M20] - m[M02]) * s,
            (m[M01] - m[M10]) * s
    );
}

INLINE Vector3 Matrix4::GetScale() const {
    return Vector3(
        std::sqrt(m[M00] * m[M00] + m[M01] * m[M01] + m[M02] * m[M02]),
        std::sqrt(m[M10] * m[M10] + m[M11] * m[M11] + m[M12] * m[M12]),
        std::sqrt(m[M20] * m[M20] + m[M21] * m[M21] + m[M22] * m[M22])
    );
}

INLINE Vector3 Matrix4::GetTranslation() const {
    return Vector3(m[M30], m[M31], m[M32]);
}

INLINE double Matrix4::Determinant() const {
    return (m[M00] * m[M11] - m[M01] * m[M10]) * (m[M22] * m[M33] - m[M23] * m[M32]) +
           (m[M02] * m[M10] - m[M00] * m[M12]) * (m[M21] * m[M33] - m[M23] * m[M31]) +
           (m[M00] * m[M13] - m[M03] * m[M10]) * (m[M21] * m[M32] - m[M22] * m[M31]) +
           (m[M01] * m[M12] - m[M02] * m[M11]) * (m[M20] * m[M33] - m[M23] * m[M30]) +
           (m[M03] * m[M11] - m[M01] * m[M13]) * (m[M20] * m[M32] - m[M22] * m[M30]) +
           (m[M02] * m[M13] - m[M03] * m[M12]) * (m[M20] * m[M31] - m[M21] * m[M30]);
}

INLINE Matrix4 &Matrix4::operator=(const Matrix4 &other) {
    std::memcpy(m, other.m, sizeof(m));
    return *this;
}

INLINE Matrix4 &Matrix4::operator*=(const Matrix4 &other) {
    if (utils::CPUFeatures::instance.HW_AVX) {
        simd::mat4d_mul(m, m, other.m);
        return *this;
    }

    double nm[16];

    nm[M00] = m[M00] * other.m[M00] +
              m[M10] * other.m[M01] +
              m[M20] * other.m[M02] +
              m[M30] * other.m[M03];

    nm[M01] = m[M01] * other.m[M00] +
              m[M11] * other.m[M01] +
              m[M21] * other.m[M02] +
              m[M31] * other.m[M03];

    nm[M02] = m[M02] * other.m[M00] +
              m[M12] * other.m[M01] +
              m[M22] * other.m[M02] +
              m[M32] * other.m[M03];

    nm[M03] = m[M03] * other.m[M00] +
              m[M13] * other.m[M01] +
              m[M23] * other.m[M02] +
              m[M33] * other.m[M03];

    nm[M10] = m[M00] * other.m[M10] +
              m[M10] * other.m[M11] +
              m[M20] * other.m[M12] +
              m[M30] * other.m[M13];

    nm[M11] = m[M01] * other.m[M10] +
              m[M11] * other.m[M11] +
              m[M21] * other.m[M12] +
              m[M31] * other.m[M13];

    nm[M12] = m[M02] * other.m[M10] +
              m[M12] * other.m[M11] +
              m[M22] * other.m[M12] +
              m[M32] * other.m[M13];

    nm[M13] = m[M03] * other.m[M10] +
              m[M13] * other.m[M11] +
              m[M23] * other.m[M12] +
              m[M33] * other.m[M13];

    nm[M20] = m[M00] * other.m[M20] +
              m[M10] * other.m[M21] +
              m[M20] * other.m[M22] +
              m[M30] * other.m[M23];

    nm[M21] = m[M01] * other.m[M20] +
              m[M11] * other.m[M21] +
              m[M21] * other.m[M22] +
              m[M31] * other.m[M23];

    nm[M22] = m[M02] * other.m[M20] +
              m[M12] * other.m[M21] +
              m[M22] * other.m[M22] +
              m[M32] * other.m[M23];

    nm[M23] = m[M03] * other.m[M20] +
              m[M13] * other.m[M21] +
              m[M23] * other.m[M22] +
              m[M33] * other.m[M23];

    nm[M30] = m[M00] * other.m[M30] +
              m[M10] * other.m[M31] +
              m[M20] * other.m[M32] +
              m[M30] * other.m[M33];

    nm[M31] = m[M01] * other.m[M30] +
              m[M11] * other.m[M31] +
              m[M21] * other.m[M32] +
              m[M31] * other.m[M33];

    nm[M32] = m[M02] * other.m[M30] +
              m[M12] * other.m[M31] +
              m[M22] * other.m[M32] +
              m[M32] * other.m[M33];

    nm[M33] = m[M03] * other.m[M30] +
              m[M13] * other.m[M31] +
              m[M23] * other.m[M32] +
              m[M33] * other.m[M33];

    std::memcpy(m, nm, sizeof(m));

    return *this;
}

INLINE double *Matrix4::operator[](int index) {
    return m + index * 4;
}

INLINE Matrix4 Matrix4::operator*(const Matrix4 &other) const {
    Matrix4 nm(m);

    if (utils::CPUFeatures::instance.HW_AVX) {
        simd::mat4d_mul(nm.m, m, other.m);
        return nm;
    }

    nm.m[M00] = m[M00] * other.m[M00] +
                m[M10] * other.m[M01] +
                m[M20] * other.m[M02] +
                m[M30] * other.m[M03];

    nm.m[M01] = m[M01] * other.m[M00] +
                m[M11] * other.m[M01] +
                m[M21] * other.m[M02] +
                m[M31] * other.m[M03];

    nm.m[M02] = m[M02] * other.m[M00] +
                m[M12] * other.m[M01] +
                m[M22] * other.m[M02] +
                m[M32] * other.m[M03];

    nm.m[M03] = m[M03] * other.m[M00] +
                m[M13] * other.m[M01] +
                m[M23] * other.m[M02] +
                m[M33] * other.m[M03];

    nm.m[M10] = m[M00] * other.m[M10] +
                m[M10] * other.m[M11] +
                m[M20] * other.m[M12] +
                m[M30] * other.m[M13];

    nm.m[M11] = m[M01] * other.m[M10] +
                m[M11] * other.m[M11] +
                m[M21] * other.m[M12] +
                m[M31] * other.m[M13];

    nm.m[M12] = m[M02] * other.m[M10] +
                m[M12] * other.m[M11] +
                m[M22] * other.m[M12] +
                m[M32] * other.m[M13];

    nm.m[M13] = m[M03] * other.m[M10] +
                m[M13] * other.m[M11] +
                m[M23] * other.m[M12] +
                m[M33] * other.m[M13];

    nm.m[M20] = m[M00] * other.m[M20] +
                m[M10] * other.m[M21] +
                m[M20] * other.m[M22] +
                m[M30] * other.m[M23];

    nm.m[M21] = m[M01] * other.m[M20] +
                m[M11] * other.m[M21] +
                m[M21] * other.m[M22] +
                m[M31] * other.m[M23];

    nm.m[M22] = m[M02] * other.m[M20] +
                m[M12] * other.m[M21] +
                m[M22] * other.m[M22] +
                m[M32] * other.m[M23];

    nm.m[M23] = m[M03] * other.m[M20] +
                m[M13] * other.m[M21] +
                m[M23] * other.m[M22] +
                m[M33] * other.m[M23];

    nm.m[M30] = m[M00] * other.m[M30] +
                m[M10] * other.m[M31] +
                m[M20] * other.m[M32] +
                m[M30] * other.m[M33];

    nm.m[M31] = m[M01] * other.m[M30] +
                m[M11] * other.m[M31] +
                m[M21] * other.m[M32] +
                m[M31] * other.m[M33];

    nm.m[M32] = m[M02] * other.m[M30] +
                m[M12] * other.m[M31] +
                m[M22] * other.m[M32] +
                m[M32] * other.m[M33];

    nm.m[M33] = m[M03] * other.m[M30] +
                m[M13] * other.m[M31] +
                m[M23] * other.m[M32] +
                m[M33] * other.m[M33];

    return nm;
}

INLINE Vector4 Matrix4::operator*(const Vector4 &other) const {
    return Vector4(
            m[M00] * other.x + m[M10] * other.y + m[M20] * other.z + m[M30] * other.w,
            m[M01] * other.x + m[M11] * other.y + m[M21] * other.z + m[M31] * other.w,
            m[M02] * other.x + m[M12] * other.y + m[M22] * other.z + m[M32] * other.w,
            m[M03] * other.x + m[M13] * other.y + m[M23] * other.z + m[M33] * other.w
    );
}

INLINE const double *Matrix4::operator[](int index) const {
    return m + index * 4;
}

} //namespace math
} //namespace sw

#undef INLINE