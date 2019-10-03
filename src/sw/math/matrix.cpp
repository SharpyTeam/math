//
// Created by selya on 24.07.2019.
//

#include "matrix.hpp"

#include <cmath>

namespace sw {
namespace math {

Matrix4::Matrix4() : data{
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
} {}

Matrix4::Matrix4(const Matrix4 &other) {
    data[0] = other.data[0];
    data[1] = other.data[1];
    data[2] = other.data[2];
    data[3] = other.data[3];

    data[4] = other.data[4];
    data[5] = other.data[5];
    data[6] = other.data[6];
    data[7] = other.data[7];

    data[8] = other.data[8];
    data[9] = other.data[9];
    data[10] = other.data[10];
    data[11] = other.data[11];

    data[12] = other.data[12];
    data[13] = other.data[13];
    data[14] = other.data[14];
    data[15] = other.data[15];
}

Matrix4::Matrix4(const double *data) {
    this->data[0] = data[0];
    this->data[1] = data[1];
    this->data[2] = data[2];
    this->data[3] = data[3];

    this->data[4] = data[4];
    this->data[5] = data[5];
    this->data[6] = data[6];
    this->data[7] = data[7];

    this->data[8] = data[8];
    this->data[9] = data[9];
    this->data[10] = data[10];
    this->data[11] = data[11];

    this->data[12] = data[12];
    this->data[13] = data[13];
    this->data[14] = data[14];
    this->data[15] = data[15];
}

Matrix4 &Matrix4::Identity() {
    data[0] = 1;
    data[1] = 0;
    data[2] = 0;
    data[3] = 0;

    data[4] = 0;
    data[5] = 1;
    data[6] = 0;
    data[7] = 0;

    data[8] = 0;
    data[9] = 0;
    data[10] = 1;
    data[11] = 0;

    data[12] = 0;
    data[13] = 0;
    data[14] = 0;
    data[15] = 1;

    return *this;
}

Matrix4 &Matrix4::Inverse() {
    double a = data[0] * data[5] - data[1] * data[4];
    double b = data[0] * data[6] - data[2] * data[4];
    double c = data[0] * data[7] - data[3] * data[4];
    double d = data[1] * data[6] - data[2] * data[5];
    double e = data[1] * data[7] - data[3] * data[5];
    double f = data[2] * data[7] - data[3] * data[6];
    double g = data[8] * data[13] - data[9] * data[12];
    double h = data[8] * data[14] - data[10] * data[12];
    double i = data[8] * data[15] - data[11] * data[12];
    double j = data[9] * data[14] - data[10] * data[13];
    double k = data[9] * data[15] - data[11] * data[13];
    double l = data[10] * data[15] - data[11] * data[14];

    double det = a * l - b * k + c * j + d * i - e * h + f * g;
    double nm00, nm01, nm02, nm03, nm10, nm11, nm12, nm13, nm20, nm21, nm22, nm23, nm30, nm31, nm32, nm33;
    det = 1.0f / det;

    nm00 = (data[5] * l - data[6] * k + data[7] * j) * det;
    nm01 = (-data[1] * l + data[2] * k - data[3] * j) * det;
    nm02 = (data[13] * f - data[14] * e + data[15] * d) * det;
    nm03 = (-data[9] * f + data[10] * e - data[11] * d) * det;
    nm10 = (-data[4] * l + data[6] * i - data[7] * h) * det;
    nm11 = (data[0] * l - data[2] * i + data[3] * h) * det;
    nm12 = (-data[12] * f + data[14] * c - data[15] * b) * det;
    nm13 = (data[8] * f - data[10] * c + data[11] * b) * det;
    nm20 = (data[4] * k - data[5] * i + data[7] * g) * det;
    nm21 = (-data[0] * k + data[1] * i - data[3] * g) * det;
    nm22 = (data[12] * e - data[13] * c + data[15] * a) * det;
    nm23 = (-data[8] * e + data[9] * c - data[11] * a) * det;
    nm30 = (-data[4] * j + data[5] * h - data[6] * g) * det;
    nm31 = (data[0] * j - data[1] * h + data[2] * g) * det;
    nm32 = (-data[12] * d + data[13] * b - data[14] * a) * det;
    nm33 = (data[8] * d - data[9] * b + data[10] * a) * det;

    data[0] = nm00;
    data[1] = nm01;
    data[2] = nm02;
    data[3] = nm03;

    data[4] = nm10;
    data[5] = nm11;
    data[6] = nm12;
    data[7] = nm13;

    data[8] = nm20;
    data[9] = nm21;
    data[10] = nm22;
    data[11] = nm23;

    data[12] = nm30;
    data[13] = nm31;
    data[14] = nm32;
    data[15] = nm33;

    return *this;
}

Matrix4 &Matrix4::Transpose() {
    double nm00 = data[0];
    double nm10 = data[1];
    double nm20 = data[2];
    double nm30 = data[3];

    double nm01 = data[4];
    double nm11 = data[5];
    double nm21 = data[6];
    double nm31 = data[7];

    double nm02 = data[8];
    double nm12 = data[9];
    double nm22 = data[10];
    double nm32 = data[11];

    double nm03 = data[12];
    double nm13 = data[13];
    double nm23 = data[14];
    double nm33 = data[15];

    data[0] = nm00;
    data[1] = nm01;
    data[2] = nm02;
    data[3] = nm03;

    data[4] = nm10;
    data[5] = nm11;
    data[6] = nm12;
    data[7] = nm13;

    data[8] = nm20;
    data[9] = nm21;
    data[10] = nm22;
    data[11] = nm23;

    data[12] = nm30;
    data[13] = nm31;
    data[14] = nm32;
    data[15] = nm33;

    return *this;
}

Matrix4 &Matrix4::Set(const Matrix4 &other) {
    data[0] = other.data[0];
    data[1] = other.data[1];
    data[2] = other.data[2];
    data[3] = other.data[3];

    data[4] = other.data[4];
    data[5] = other.data[5];
    data[6] = other.data[6];
    data[7] = other.data[7];

    data[8] = other.data[8];
    data[9] = other.data[9];
    data[10] = other.data[10];
    data[11] = other.data[11];

    data[12] = other.data[12];
    data[13] = other.data[13];
    data[14] = other.data[14];
    data[15] = other.data[15];

    return *this;
}

Matrix4 &Matrix4::Set(const double *data) {
    this->data[0] = data[0];
    this->data[1] = data[1];
    this->data[2] = data[2];
    this->data[3] = data[3];

    this->data[4] = data[4];
    this->data[5] = data[5];
    this->data[6] = data[6];
    this->data[7] = data[7];

    this->data[8] = data[8];
    this->data[9] = data[9];
    this->data[10] = data[10];
    this->data[11] = data[11];

    this->data[12] = data[12];
    this->data[13] = data[13];
    this->data[14] = data[14];
    this->data[15] = data[15];

    return *this;
}

Matrix4 &Matrix4::Scale(double x, double y, double z) {
    this->data[0] *= x;
    this->data[1] *= x;
    this->data[2] *= x;
    this->data[3] *= x;

    this->data[4] *= y;
    this->data[5] *= y;
    this->data[6] *= y;
    this->data[7] *= y;

    this->data[8] *= z;
    this->data[9] *= z;
    this->data[10] *= z;
    this->data[11] *= z;

    return *this;
}

Matrix4 &Matrix4::Scale(const Vector3 &scale) {
    return Matrix4::Scale(scale.x, scale.y, scale.z);
}

Matrix4 &Matrix4::Rotate(double ang, double x, double y, double z) {
    double sin = std::sin(ang);
    double cos;

    // TODO move to function
    static constexpr double PI = 3.14159265358979323846;
    static constexpr double PI2 = PI * 2.0;
    static constexpr double PIHalf = PI * 0.5;

    double cos_v = std::sqrt(1.0 - sin * sin),
            a = ang + PIHalf,
            b = a - int(a / PI2) * PI2;

    if (b < 0.0)
        b = PI2 + b;
    cos = (double) (b >= PI ? -cos_v : cos_v);
    //

    double C = 1.0f - cos;

    double xx = x * x, xy = x * y, xz = x * z;
    double yy = y * y, yz = y * z;
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

    double nm00 = this->data[0] * rm00 + this->data[4] * rm01 + this->data[8] * rm02;
    double nm01 = this->data[1] * rm00 + this->data[5] * rm01 + this->data[9] * rm02;
    double nm02 = this->data[2] * rm00 + this->data[6] * rm01 + this->data[10] * rm02;
    double nm03 = this->data[3] * rm00 + this->data[7] * rm01 + this->data[11] * rm02;
    double nm10 = this->data[0] * rm10 + this->data[4] * rm11 + this->data[8] * rm12;
    double nm11 = this->data[1] * rm10 + this->data[5] * rm11 + this->data[9] * rm12;
    double nm12 = this->data[2] * rm10 + this->data[6] * rm11 + this->data[10] * rm12;
    double nm13 = this->data[3] * rm10 + this->data[7] * rm11 + this->data[11] * rm12;

    this->data[8] = this->data[0] * rm20 + this->data[4] * rm21 + this->data[8] * rm22;
    this->data[9] = this->data[1] * rm20 + this->data[5] * rm21 + this->data[9] * rm22;
    this->data[10] = this->data[2] * rm20 + this->data[6] * rm21 + this->data[10] * rm22;
    this->data[11] = this->data[3] * rm20 + this->data[7] * rm21 + this->data[11] * rm22;

    this->data[0] = nm00;
    this->data[1] = nm01;
    this->data[2] = nm02;
    this->data[3] = nm03;
    this->data[4] = nm10;
    this->data[5] = nm11;
    this->data[6] = nm12;
    this->data[7] = nm13;

    return *this;
}

Matrix4 &Matrix4::Rotate(double ang, const Vector3 &axis) {
    return Matrix4::Rotate(ang, axis.x, axis.y, axis.z);
}

Matrix4 &Matrix4::Translate(double x, double y, double z) {
    this->data[12] = data[0] * x + data[4] * y + data[8] * z + data[12];
    this->data[13] = data[1] * x + data[5] * y + data[9] * z + data[13];
    this->data[14] = data[2] * x + data[6] * y + data[10] * z + data[14];
    this->data[15] = data[3] * x + data[7] * y + data[11] * z + data[15];

    return *this;
}

Matrix4 &Matrix4::Translate(const Vector3 &translate) {
    return Matrix4::Translate(translate.x, translate.y, translate.z);
}

Matrix4 &Matrix4::SetOrtho2D(double left, double right, double bottom, double top) {
    data[0] = 2.0f / (right - left);
    data[1] = 0;
    data[2] = 0;
    data[3] = 0;

    data[4] = 0;
    data[5] = 2.0f / (top - bottom);
    data[6] = 0;
    data[7] = 0;

    data[8] = 0;
    data[9] = 0;
    data[10] = -1.0f;
    data[11] = 0;

    data[12] = (right + left) / (left - right);
    data[13] = (top + bottom) / (bottom - top);
    data[14] = 0;
    data[15] = 1;

    return *this;
}

Matrix4 &Matrix4::operator=(const Matrix4 &other) {
    Set(other);
    return *this;
}

double *Matrix4::operator[](int index) {
    return &data[index * 4];
}

Matrix4 Matrix4::operator*(const Matrix4 &other) const {
    Matrix4 nm;
    nm.data[0] = data[0] * other.data[0] +
                 data[4] * other.data[1] +
                 data[8] * other.data[2] +
                 data[12] * other.data[3];

    nm.data[1] = data[1] * other.data[0] +
                 data[5] * other.data[1] +
                 data[9] * other.data[2] +
                 data[13] * other.data[3];

    nm.data[2] = data[2] * other.data[0] +
                 data[6] * other.data[1] +
                 data[10] * other.data[2] +
                 data[14] * other.data[3];

    nm.data[3] = data[3] * other.data[0] +
                 data[7] * other.data[1] +
                 data[11] * other.data[2] +
                 data[15] * other.data[3];

    nm.data[4] = data[0] * other.data[4] +
                 data[4] * other.data[5] +
                 data[8] * other.data[6] +
                 data[12] * other.data[7];

    nm.data[5] = data[1] * other.data[4] +
                 data[5] * other.data[5] +
                 data[9] * other.data[6] +
                 data[13] * other.data[7];

    nm.data[6] = data[2] * other.data[4] +
                 data[6] * other.data[5] +
                 data[10] * other.data[6] +
                 data[14] * other.data[7];

    nm.data[7] = data[3] * other.data[4] +
                 data[7] * other.data[5] +
                 data[11] * other.data[6] +
                 data[15] * other.data[7];

    nm.data[8] = data[0] * other.data[8] +
                 data[4] * other.data[9] +
                 data[8] * other.data[10] +
                 data[12] * other.data[11];

    nm.data[9] = data[1] * other.data[8] +
                 data[5] * other.data[9] +
                 data[9] * other.data[10] +
                 data[13] * other.data[11];

    nm.data[10] = data[2] * other.data[8] +
                  data[6] * other.data[9] +
                  data[10] * other.data[10] +
                  data[14] * other.data[11];

    nm.data[11] = data[3] * other.data[8] +
                  data[7] * other.data[9] +
                  data[11] * other.data[10] +
                  data[15] * other.data[11];

    nm.data[12] = data[0] * other.data[12] +
                  data[4] * other.data[13] +
                  data[8] * other.data[14] +
                  data[12] * other.data[15];

    nm.data[13] = data[1] * other.data[12] +
                  data[5] * other.data[13] +
                  data[9] * other.data[14] +
                  data[13] * other.data[15];

    nm.data[14] = data[2] * other.data[12] +
                  data[6] * other.data[13] +
                  data[10] * other.data[14] +
                  data[14] * other.data[15];

    nm.data[15] = data[3] * other.data[12] +
                  data[7] * other.data[13] +
                  data[11] * other.data[14] +
                  data[15] * other.data[15];

    return nm;
}

Vector4 Matrix4::operator*(const Vector4 &v) const {
    return Vector4(
        data[0] * v.x + data[4] * v.y + data[8] * v.z + data[12] * v.w,
        data[1] * v.x + data[5] * v.y + data[9] * v.z + data[13] * v.w,
        data[2] * v.x + data[6] * v.y + data[10] * v.z + data[14] * v.w,
        data[3] * v.x + data[7] * v.y + data[11] * v.z + data[15] * v.w
    );
}

const double *Matrix4::operator[](int index) const {
    return &data[index * 4];
}

/*template<class double>
Vector<4, double> Matrix<4, 4, double>::Unproject(double x, double y, double z) const {
    Vector<4, double> v;

    double a = data[0]  * data[5]  - data[1]  * data[4];
    double b = data[0]  * data[6]  - data[2]  * data[4];
    double c = data[0]  * data[7]  - data[3]  * data[4];
    double d = data[1]  * data[6]  - data[2]  * data[5];
    double e = data[1]  * data[7]  - data[3]  * data[5];
    double f = data[2]  * data[7]  - data[3]  * data[6];
    double g = data[8]  * data[13] - data[9]  * data[12];
    double h = data[8]  * data[14] - data[10] * data[12];
    double i = data[8]  * data[15] - data[11] * data[12];
    double j = data[9]  * data[14] - data[10] * data[13];
    double k = data[9]  * data[15] - data[11] * data[13];
    double l = data[10] * data[15] - data[11] * data[14];

    double det = a * l - b * k + c * j + d * i - e * h + f * g;
    det = 1.0f / det;

    double im00 = ( data[5]  * l - data[6]  * k + data[7]  * j) * det;
    double im01 = (-data[1]  * l + data[2]  * k - data[3]  * j) * det;
    double im02 = ( data[13] * f - data[14] * e + data[15] * d) * det;
    double im03 = (-data[9]  * f + data[10] * e - data[11] * d) * det;
    double im10 = (-data[4]  * l + data[6]  * i - data[7]  * h) * det;
    double im11 = ( data[0]  * l - data[2]  * i + data[3]  * h) * det;
    double im12 = (-data[12] * f + data[14] * c - data[15] * b) * det;
    double im13 = ( data[8]  * f - data[10] * c + data[11] * b) * det;
    double im20 = ( data[4]  * k - data[5]  * i + data[7]  * g) * det;
    double im21 = (-data[0]  * k + data[1]  * i - data[3]  * g) * det;
    double im22 = ( data[12] * e - data[13] * c + data[15] * a) * det;
    double im23 = (-data[8]  * e + data[9]  * c - data[11] * a) * det;
    double im30 = (-data[4]  * j + data[5]  * h - data[6]  * g) * det;
    double im31 = ( data[0]  * j - data[1]  * h + data[2]  * g) * det;
    double im32 = (-data[12] * d + data[13] * b - data[14] * a) * det;
    double im33 = ( data[8]  * d - data[9]  * b + data[10] * a) * det;

    // Assume that it is NDC
    double ndcX = x;
    double ndcY = y;
    double ndcZ = z;

    double invW = 1.0f / (im03 * ndcX + im13 * ndcY + im23 * ndcZ + im33);

    v.x = (im00 * ndcX + im10 * ndcY + im20 * ndcZ + im30) * invW;
    v.y = (im01 * ndcX + im11 * ndcY + im21 * ndcZ + im31) * invW;
    v.z = (im02 * ndcX + im12 * ndcY + im22 * ndcZ + im32) * invW;
    v.w = 1.0f;

    return v;
}

template<class double>
Vector<4, double> Matrix<4, 4, double>::Unproject(const Vector<3, double> &v) const {
    return Unproject(v.x, v.y, v.z);
}*/

Vector3 Matrix4::GetTranslation() const {
    return Vector3(data[12], data[13], data[14]);
}

double Matrix4::GetRotationZ() const {
    Vector4 rot = *this * Vector4(1, 0, 0, 0);
    return std::atan2(rot.y, rot.x);
}

} //namespace math
} //namespace sw
