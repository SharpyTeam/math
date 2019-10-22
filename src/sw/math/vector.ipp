//
// Created by Ilya on 30.06.2019.
//

#include "vector.hpp"

#include <cmath>
#include <cassert>

#if SANDWICH_MATH_INLINE
    #define INLINE inline
#else
    #define INLINE
#endif

namespace sw {
namespace math {

INLINE Vector2::Vector2() : x(0.0), y(0.0) {}
INLINE Vector2::Vector2(double scalar) : x(scalar), y(scalar) {}
INLINE Vector2::Vector2(double x, double y) : x(x), y(y) {}
INLINE Vector2::Vector2(const Vector2 &other) : x(other.x), y(other.y) {}

INLINE double Vector2::Angle(double x_, double y_) const {
    return std::abs(std::atan2(Determinant(x_, y_), Dot(x_, y_)));
}

INLINE double Vector2::Angle(const Vector2 &other) const {
    return Angle(other.x, other.y);
}

INLINE double Vector2::Determinant(double x_, double y_) const {
    return x * y_ - y * x_;
}

INLINE double Vector2::Determinant(const Vector2 &other) const {
    return Determinant(other.x, other.y);
}

INLINE double Vector2::Distance(double x_, double y_) const {
    return std::sqrt(DistanceSquared(x_, y_));
}

INLINE double Vector2::Distance(const Vector2 &other) const {
    return Distance(other.x, other.y);
}

INLINE double Vector2::DistanceSquared(double x_, double y_) const {
    double x_d = x - x_;
    double y_d = y - y_;
    return x_d * x_d + y_d * y_d;
}

INLINE double Vector2::DistanceSquared(const Vector2 &other) const {
    return DistanceSquared(other.x, other.y);
}

INLINE double Vector2::Dot(double x_, double y_) const {
    return x * x_ + y * y_;
}

INLINE double Vector2::Dot(const Vector2 &other) const {
    return Dot(other.x, other.y);
}

INLINE double Vector2::Length() const {
    return std::sqrt(LengthSquared());
}

INLINE double Vector2::LengthSquared() const {
    return x * x + y * y;
}

INLINE Vector2 Vector2::Lerp(double x_, double y_, double t) const {
    return Vector2(x + (x_ - x) * t, y + (y_ - y) * t);
}

INLINE Vector2 Vector2::Lerp(const Vector2 &other, double t) const {
    return Lerp(other.x, other.y, t);
}

INLINE Vector2 Vector2::Max(double x_, double y_) const {
    return Vector2(x > x_ ? x : x_, y > y_ ? y : y_);
}

INLINE Vector2 Vector2::Max(const Vector2 &other) const {
    return Max(other.x, other.y);
}

INLINE Vector2 Vector2::Min(double x_, double y_) const {
    return Vector2(x < x_ ? x : x_, y < y_ ? y : y_);
}

INLINE Vector2 Vector2::Min(const Vector2 &other) const {
    return Min(other.x, other.y);
}

INLINE Vector2 Vector2::Normalized() const {
    double l = Length();
    assert(l != 0.0);
    return Vector2(x / l, y / l);
}

INLINE Vector2 Vector2::Perpendicular() const {
    return Vector2(y, -x);
}

INLINE Vector2 Vector2::operator-() const {
    return Vector2(-x, -y);
}

INLINE Vector2 Vector2::operator+(const Vector2 &other) const {
    return Vector2(x + other.x, y + other.y);
}

INLINE Vector2 Vector2::operator-(const Vector2 &other) const {
    return Vector2(x - other.x, y - other.y);
}

INLINE Vector2 Vector2::operator*(const Vector2 &other) const {
    return Vector2(x * other.x, y * other.y);
}

INLINE Vector2 Vector2::operator/(const Vector2 &other) const {
    return Vector2(x / other.x, y / other.y);
}

INLINE Vector2 Vector2::operator+(double scalar) const {
    return Vector2(x + scalar, y + scalar);
}

INLINE Vector2 Vector2::operator-(double scalar) const {
    return Vector2(x - scalar, y - scalar);
}

INLINE Vector2 Vector2::operator*(double scalar) const {
    return Vector2(x * scalar, y * scalar);
}

INLINE Vector2 Vector2::operator/(double scalar) const {
    return Vector2(x / scalar, y / scalar);
}

INLINE bool Vector2::operator==(double scalar) const {
    // use this operator only with integers
    assert(std::trunc(scalar) == scalar);
    assert(std::trunc(x) == x);
    assert(std::trunc(y) == y);
    return x == scalar && y == scalar;
}

INLINE bool Vector2::operator==(const Vector2 &other) const {
    // use this operator only with integers
    assert(std::trunc(other.x) == other.x);
    assert(std::trunc(other.y) == other.y);
    assert(std::trunc(x) == x);
    assert(std::trunc(y) == y);
    return x == other.x && y == other.y;
}

INLINE Vector2 &Vector2::Ceil() {
    x = std::ceil(x);
    y = std::ceil(y);
    return *this;
}

INLINE Vector2 &Vector2::Floor() {
    x = std::floor(x);
    y = std::floor(y);
    return *this;
}

INLINE Vector2 &Vector2::Normalize() {
    double l = Length();
    assert(l != 0.0);
    x /= l;
    y /= l;
    return *this;
}

INLINE Vector2 &Vector2::Round() {
    x = std::round(x);
    y = std::round(y);
    return *this;
}

INLINE Vector2 &Vector2::Set(double scalar) {
    x = scalar;
    y = scalar;
    return *this;
}

INLINE Vector2 &Vector2::Set(double x_, double y_) {
    x = x_;
    y = y_;
    return *this;
}

INLINE Vector2 &Vector2::Set(const Vector2 &other) {
    return Set(other.x, other.y);
}

INLINE Vector2 &Vector2::Zero() {
    x = 0.0;
    y = 0.0;
    return *this;
}

INLINE Vector2 &Vector2::operator=(double scalar) {
    x = scalar;
    y = scalar;
    return *this;
}

INLINE Vector2 &Vector2::operator+=(double scalar) {
    x += scalar;
    y += scalar;
    return *this;
}

INLINE Vector2 &Vector2::operator-=(double scalar) {
    x -= scalar;
    y -= scalar;
    return *this;
}

INLINE Vector2 &Vector2::operator*=(double scalar) {
    x *= scalar;
    y *= scalar;
    return *this;
}

INLINE Vector2 &Vector2::operator/=(double scalar) {
    x /= scalar;
    y /= scalar;
    return *this;
}

INLINE Vector2 &Vector2::operator=(const Vector2 &other) {
    x = other.x;
    y = other.y;
    return *this;
}

INLINE Vector2 &Vector2::operator+=(const Vector2 &other) {
    x += other.x;
    y += other.y;
    return *this;
}

INLINE Vector2 &Vector2::operator-=(const Vector2 &other) {
    x -= other.x;
    y -= other.y;
    return *this;
}

INLINE Vector2 &Vector2::operator*=(const Vector2 &other) {
    x *= other.x;
    y *= other.y;
    return *this;
}

INLINE Vector2 &Vector2::operator/=(const Vector2 &other) {
    x /= other.x;
    y /= other.y;
    return *this;
}

INLINE Vector3::Vector3() : x(0.0), y(0.0), z(0.0) {}
INLINE Vector3::Vector3(double scalar) : x(scalar), y(scalar), z(scalar) {}
INLINE Vector3::Vector3(double x, double y, double z) : x(x), y(y), z(z) {}
INLINE Vector3::Vector3(const Vector3 &other) : x(other.x), y(other.y), z(other.z) {}
INLINE Vector3::Vector3(const Vector2 &other, double z) : x(other.x), y(other.y), z(z) {}

INLINE double Vector3::Angle(double x_, double y_, double z_) const {
    return std::acos(AngleCos(x_, y_, z_));
}

INLINE double Vector3::Angle(const Vector3 &other) const {
    return Angle(other.x, other.y, other.z);
}

INLINE double Vector3::AngleCos(double x_, double y_, double z_) const {
    double cos = (x * x_ + y * y_ + z * z_) /
                 std::sqrt((x * x + y * y + z * z) * (x_ * x_ + y_ * y_ * z_ * z_));
    // clamp it (can be out of [-1, 1] cause of rounding)
    cos = cos < 1 ? cos : 1;
    cos = cos > -1 ? cos : -1;
    return cos;
}

INLINE double Vector3::AngleCos(const Vector3 &other) const {
    return AngleCos(other.x, other.y, other.z);
}

INLINE double Vector3::Distance(double x_, double y_, double z_) const {
    return std::sqrt(DistanceSquared(x_, y_, z_));
}

INLINE double Vector3::Distance(const Vector3 &other) const {
    return Distance(other.x, other.y, other.z);
}

INLINE double Vector3::DistanceSquared(double x_, double y_, double z_) const {
    double x_d = x - x_;
    double y_d = y - y_;
    double z_d = z - z_;
    return x_d * x_d + y_d * y_d + z_d * z_d;
}

INLINE double Vector3::DistanceSquared(const Vector3 &other) const {
    return DistanceSquared(other.x, other.y, other.z);
}

INLINE double Vector3::Dot(double x_, double y_, double z_) const {
    return x * x_ + y * y_ + z * z_;
}

INLINE double Vector3::Dot(const Vector3 &other) const {
    return Dot(other.x, other.y, other.z);
}

INLINE double Vector3::Length() const {
    return std::sqrt(LengthSquared());
}

INLINE double Vector3::LengthSquared() const {
    return x * x + y * y + z * z;
}

INLINE Vector3 Vector3::Cross(double x_, double y_, double z_) const {
    return Vector3(y * z_ - z * y_, z * x_ - x * z_, x * y_ - y * x_);
}

INLINE Vector3 Vector3::Cross(const Vector3 &other) const {
    return Cross(other.x, other.y, other.z);
}

INLINE Vector3 Vector3::Lerp(double x_, double y_, double z_, double t) const {
    return Vector3(x + (x_ - x) * t, y + (y_ - y) * t,z + (z_ - z) * t);
}

INLINE Vector3 Vector3::Lerp(const Vector3 &other, double t) const {
    return Lerp(other.x, other.y, other.z, t);
}

INLINE Vector3 Vector3::Max(double x_, double y_, double z_) const {
    return Vector3(x > x_ ? x : x_, y > y_ ? y : y_, z > z_ ? z : z_);
}

INLINE Vector3 Vector3::Max(const Vector3 &other) const {
    return Max(other.x, other.y, other.z);
}

INLINE Vector3 Vector3::Min(double x_, double y_, double z_) const {
    return Vector3(x < x_ ? x : x_, y < y_ ? y : y_, z < z_ ? z : z_);
}

INLINE Vector3 Vector3::Min(const Vector3 &other) const {
    return Min(other.x, other.y, other.z);
}

INLINE Vector3 Vector3::Normalized() const {
    double l = Length();
    assert(l != 0.0);
    return Vector3(x / l, y / l, z / l);
}

INLINE Vector3 Vector3::operator-() const {
    return Vector3(-x, -y, -z);
}

INLINE Vector3 Vector3::operator+(const Vector3 &other) const {
    return Vector3(x + other.x, y + other.y, z + other.z);
}

INLINE Vector3 Vector3::operator-(const Vector3 &other) const {
    return Vector3(x - other.x, y - other.y, z - other.z);
}

INLINE Vector3 Vector3::operator*(const Vector3 &other) const {
    return Vector3(x * other.x, y * other.y, z * other.z);
}

INLINE Vector3 Vector3::operator/(const Vector3 &other) const {
    return Vector3(x / other.x, y / other.y, z / other.z);
}

INLINE Vector3 Vector3::operator+(double scalar) const {
    return Vector3(x + scalar, y + scalar, z + scalar);
}

INLINE Vector3 Vector3::operator-(double scalar) const {
    return Vector3(x - scalar, y - scalar, z - scalar);
}

INLINE Vector3 Vector3::operator*(double scalar) const {
    return Vector3(x * scalar, y * scalar, z * scalar);
}

INLINE Vector3 Vector3::operator/(double scalar) const {
    return Vector3(x / scalar, y / scalar, z / scalar);
}

INLINE bool Vector3::operator==(double scalar) const {
    // use this operator only with integers
    assert(std::trunc(scalar) == scalar);
    assert(std::trunc(x) == x);
    assert(std::trunc(y) == y);
    assert(std::trunc(z) == z);
    return x == scalar && y == scalar && z == scalar;
}

INLINE bool Vector3::operator==(const Vector3 &other) const {
    // use this operator only with integers
    assert(std::trunc(other.x) == other.x);
    assert(std::trunc(other.y) == other.y);
    assert(std::trunc(other.z) == other.z);
    assert(std::trunc(x) == x);
    assert(std::trunc(y) == y);
    assert(std::trunc(z) == z);
    return x == other.x && y == other.y && z == other.z;
}

INLINE Vector3 &Vector3::Ceil() {
    x = std::ceil(x);
    y = std::ceil(y);
    z = std::ceil(z);
    return *this;
}

INLINE Vector3 &Vector3::Floor() {
    x = std::floor(x);
    y = std::floor(y);
    z = std::floor(z);
    return *this;
}

INLINE Vector3 &Vector3::Normalize() {
    double l = Length();
    assert(l != 0.0);
    x /= l;
    y /= l;
    z /= l;
    return *this;
}

INLINE Vector3 &Vector3::Round() {
    x = std::round(x);
    y = std::round(y);
    z = std::round(z);
    return *this;
}

INLINE Vector3 &Vector3::Set(double scalar) {
    x = scalar;
    y = scalar;
    z = scalar;
    return *this;
}

INLINE Vector3 &Vector3::Set(double x_, double y_, double z_) {
    x = x_;
    y = y_;
    z = z_;
    return *this;
}

INLINE Vector3 &Vector3::Set(const Vector3 &other) {
    return Set(other.x, other.y, other.z);
}

INLINE Vector3 &Vector3::Set(const Vector2 &other, double z_) {
    return Set(other.x, other.y, z_);
}

INLINE Vector3 &Vector3::Zero() {
    x = 0;
    y = 0;
    z = 0;
    return *this;
}

INLINE Vector3 &Vector3::operator=(const Vector3 &other) {
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
}

INLINE Vector3 &Vector3::operator=(double scalar) {
    x = scalar;
    y = scalar;
    z = scalar;
    return *this;
}

INLINE Vector3 &Vector3::operator+=(double scalar) {
    x += scalar;
    y += scalar;
    z += scalar;
    return *this;
}

INLINE Vector3 &Vector3::operator-=(double scalar) {
    x -= scalar;
    y -= scalar;
    z -= scalar;
    return *this;
}

INLINE Vector3 &Vector3::operator*=(double scalar) {
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
}

INLINE Vector3 &Vector3::operator/=(double scalar) {
    x /= scalar;
    y /= scalar;
    z /= scalar;
    return *this;
}

INLINE Vector3 &Vector3::operator+=(const Vector3 &other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

INLINE Vector3 &Vector3::operator-=(const Vector3 &other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
}

INLINE Vector3 &Vector3::operator*=(const Vector3 &other) {
    x *= other.x;
    y *= other.y;
    z *= other.z;
    return *this;
}

INLINE Vector3 &Vector3::operator/=(const Vector3 &other) {
    x /= other.x;
    y /= other.y;
    z /= other.z;
    return *this;
}

INLINE Vector4::Vector4() : x(0.0), y(0.0), z(0.0), w(0.0) {}
INLINE Vector4::Vector4(double scalar) : x(scalar), y(scalar), z(scalar), w(scalar) {}
INLINE Vector4::Vector4(double x, double y, double z, double w) : x(x), y(y), z(z), w(w) {}
INLINE Vector4::Vector4(const Vector4 &other) : x(other.x), y(other.y), z(other.z), w(other.w) {}
INLINE Vector4::Vector4(const Vector3 &other, double w) : x(other.x), y(other.y), z(other.z), w(w) {}
INLINE Vector4::Vector4(const Vector2 &other, double z, double w) : x(other.x), y(other.y), z(z), w(w) {}

INLINE double Vector4::Angle(double x_, double y_, double z_, double w_) const {
    return std::acos(AngleCos(x_, y_, z_, w_));
}

INLINE double Vector4::Angle(const Vector4 &other) const {
    return Angle(other.x, other.y, other.z, other.w);
}

INLINE double Vector4::AngleCos(double x_, double y_, double z_, double w_) const {
    double cos = (x * x_ + y * y_ + z * z_ + w * w_) /
        std::sqrt((x * x + y * y + z * z + w * w) * (x_ * x_ + y_ * y_ * z_ * z_ + w_ * w_));
    // clamp it (can be out of [-1, 1] cause of rounding)
    cos = cos < 1 ? cos : 1;
    cos = cos > -1 ? cos : -1;
    return cos;
}

INLINE double Vector4::AngleCos(const Vector4 &other) const {
    return AngleCos(other.x, other.y, other.z, other.w);
}

INLINE double Vector4::Distance(double x_, double y_, double z_, double w_) const {
    return std::sqrt(DistanceSquared(x_, y_, z_, w_));
}

INLINE double Vector4::Distance(const Vector4 &other) const {
    return Distance(other.x, other.y, other.z, other.w);
}

INLINE double Vector4::DistanceSquared(double x_, double y_, double z_, double w_) const {
    double x_d = x - x_;
    double y_d = y - y_;
    double z_d = z - z_;
    double w_d = w - w_;
    return x_d * x_d + y_d * y_d + z_d * z_d + w_d * w_d;
}

INLINE double Vector4::DistanceSquared(const Vector4 &other) const {
    return DistanceSquared(other.x, other.y, other.z, other.w);
}

INLINE double Vector4::Dot(double x_, double y_, double z_, double w_) const {
    return x * x_ + y * y_ + z * z_ + w * w_;
}

INLINE double Vector4::Dot(const Vector4 &other) const {
    return Dot(other.x, other.y, other.z, other.w);
}

INLINE double Vector4::Length() const {
    return std::sqrt(LengthSquared());
}

INLINE double Vector4::LengthSquared() const {
    return x * x + y * y + z * z + w * w;
}

INLINE Vector4 Vector4::Lerp(double x_, double y_, double z_, double w_, double t) const {
    return Vector4(x + (x_ - x) * t, y + (y_ - y) * t, z + (z_ - z) * t, w + (w_ - w) * t);
}

INLINE Vector4 Vector4::Lerp(const Vector4 &other, double t) const {
    return Lerp(other.x, other.y, other.z, other.w, t);
}

INLINE Vector4 Vector4::Max(double x_, double y_, double z_, double w_) const {
    return Vector4(x > x_ ? x : x_, y > y_ ? y : y_, z > z_ ? z : z_, w > w_ ? w : w_);
}

INLINE Vector4 Vector4::Max(const Vector4 &other) const {
    return Max(other.x, other.y, other.z, other.w);
}

INLINE Vector4 Vector4::Min(double x_, double y_, double z_, double w_) const {
    return Vector4(x < x_ ? x : x_, y < y_ ? y : y_, z < z_ ? z : z_, w < w_ ? w : w_);
}

INLINE Vector4 Vector4::Min(const Vector4 &other) const {
    return Min(other.x, other.y, other.z, other.w);
}

INLINE Vector4 Vector4::Normalized() const {
    double l = Length();
    assert(l != 0.0);
    return Vector4(x / l, y / l, z / l, w / l);
}

INLINE Vector4 Vector4::operator-() const {
    return Vector4(-x, -y, -z, -w);
}

INLINE Vector4 Vector4::operator+(const Vector4 &other) const {
    return Vector4(x + other.x, y + other.y, z + other.z, w + other.w);
}

INLINE Vector4 Vector4::operator-(const Vector4 &other) const {
    return Vector4(x - other.x, y - other.y, z - other.z, w - other.w);
}

INLINE Vector4 Vector4::operator*(const Vector4 &other) const {
    return Vector4(x * other.x, y * other.y, z * other.z, w * other.w);
}

INLINE Vector4 Vector4::operator/(const Vector4 &other) const {
    return Vector4(x / other.x, y / other.y, z / other.z, w / other.w);
}

INLINE Vector4 Vector4::operator+(double scalar) const {
    return Vector4(x + scalar, y + scalar, z + scalar, w + scalar);
}

INLINE Vector4 Vector4::operator-(double scalar) const {
    return Vector4(x - scalar, y - scalar, z - scalar, w - scalar);
}

INLINE Vector4 Vector4::operator*(double scalar) const {
    return Vector4(x * scalar, y * scalar, z * scalar, w * scalar);
}

INLINE Vector4 Vector4::operator/(double scalar) const {
    return Vector4(x / scalar, y / scalar, z / scalar, w / scalar);
}

INLINE bool Vector4::operator==(double scalar) const {
    // use this operator only with integers
    assert(std::trunc(scalar) == scalar);
    assert(std::trunc(x) == x);
    assert(std::trunc(y) == y);
    assert(std::trunc(z) == z);
    assert(std::trunc(w) == w);
    return x == scalar && y == scalar && z == scalar && w == scalar;
}

INLINE bool Vector4::operator==(const Vector4 &other) const {
    // use this operator only with integers
    assert(std::trunc(other.x) == other.x);
    assert(std::trunc(other.y) == other.y);
    assert(std::trunc(other.z) == other.z);
    assert(std::trunc(other.w) == other.w);
    assert(std::trunc(x) == x);
    assert(std::trunc(y) == y);
    assert(std::trunc(z) == z);
    assert(std::trunc(w) == w);
    return x == other.x && y == other.y && z == other.z && w == other.w;
}

INLINE Vector4 &Vector4::Ceil() {
    x = std::ceil(x);
    y = std::ceil(y);
    z = std::ceil(z);
    w = std::ceil(w);
    return *this;
}

INLINE Vector4 &Vector4::Floor() {
    x = std::floor(x);
    y = std::floor(y);
    z = std::floor(z);
    w = std::floor(w);
    return *this;
}

INLINE Vector4 &Vector4::Normalize() {
    double l = Length();
    assert(l != 0.0);
    x /= l;
    y /= l;
    z /= l;
    w /= l;
    return *this;
}

INLINE Vector4 &Vector4::Round() {
    x = std::round(x);
    y = std::round(y);
    z = std::round(z);
    w = std::round(w);
    return *this;
}

INLINE Vector4 &Vector4::Set(double scalar) {
    x = scalar;
    y = scalar;
    z = scalar;
    w = scalar;
    return *this;
}

INLINE Vector4 &Vector4::Set(double x_, double y_, double z_, double w_) {
    x = x_;
    y = y_;
    z = z_;
    w = w_;
    return *this;
}

INLINE Vector4 &Vector4::Set(const Vector4 &other) {
    return Set(other.x, other.y, other.z, other.w);
}

INLINE Vector4 &Vector4::Set(const Vector3 &other, double w_) {
    return Set(other.x, other.y, other.z, w_);
}

INLINE Vector4 &Vector4::Set(const Vector2 &other, double z_, double w_) {
    return Set(other.x, other.y, z_, w_);
}

INLINE Vector4 &Vector4::Zero() {
    x = 0;
    y = 0;
    z = 0;
    w = 0;
    return *this;
}

INLINE Vector4 &Vector4::operator=(const Vector4 &other) {
    x = other.x;
    y = other.y;
    z = other.z;
    w = other.w;
    return *this;
}

INLINE Vector4 &Vector4::operator=(double scalar) {
    x = scalar;
    y = scalar;
    z = scalar;
    w = scalar;
    return *this;
}

INLINE Vector4 &Vector4::operator+=(double scalar) {
    x += scalar;
    y += scalar;
    z += scalar;
    w += scalar;
    return *this;
}

INLINE Vector4 &Vector4::operator-=(double scalar) {
    x -= scalar;
    y -= scalar;
    z -= scalar;
    w -= scalar;
    return *this;
}

INLINE Vector4 &Vector4::operator*=(double scalar) {
    x *= scalar;
    y *= scalar;
    z *= scalar;
    w *= scalar;
    return *this;
}

INLINE Vector4 &Vector4::operator/=(double scalar) {
    x /= scalar;
    y /= scalar;
    z /= scalar;
    w /= scalar;
    return *this;
}

INLINE Vector4 &Vector4::operator+=(const Vector4 &other) {
    x += other.x;
    y += other.y;
    z += other.z;
    w += other.w;
    return *this;
}

INLINE Vector4 &Vector4::operator-=(const Vector4 &other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    w -= other.w;
    return *this;
}

INLINE Vector4 &Vector4::operator*=(const Vector4 &other) {
    x *= other.x;
    y *= other.y;
    z *= other.z;
    w *= other.w;
    return *this;
}

INLINE Vector4 &Vector4::operator/=(const Vector4 &other) {
    x /= other.x;
    y /= other.y;
    z /= other.z;
    w /= other.w;
    return *this;
}

} //namespace math
} //namespace sw

#undef INLINE