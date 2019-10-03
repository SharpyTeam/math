//
// Created by Ilya on 30.06.2019.
//

#include "vector.hpp"

#include <cmath>
#include <string>

namespace sw {
namespace math {
// ----------------------------------------------- Vector2 ------------------------------------------------

Vector2::Vector2() : x(0.0), y(0.0) {}

Vector2::Vector2(double x, double y) : x(x), y(y) {}

Vector2::Vector2(double scalar) : x(scalar), y(scalar) {}

Vector2::Vector2(const Vector2 &other) : x(other.x), y(other.y) {}


double Vector2::Length() const {
    return std::sqrt(x * x + y * y);
}


double Vector2::LengthSquared() const {
    return x * x + y * y;
}


double Vector2::Distance(double tox, double toy) const {
    double x_diff = this->x - tox;
    double y_diff = this->y - toy;
    return std::sqrt(x_diff * x_diff + y_diff * y_diff);
}


double Vector2::Distance(const Vector2 &other) const {
    double x_diff = x - other.x;
    double y_diff = y - other.y;
    return std::sqrt(x_diff * x_diff + y_diff * y_diff);
}


double Vector2::DistanceSquared(double tox, double toy) const {
    double x_diff = this->x - tox;
    double y_diff = this->y - toy;
    return x_diff * x_diff + y_diff * y_diff;
}


double Vector2::DistanceSquared(const Vector2 &other) const {
    double x_diff = x - other.x;
    double y_diff = y - other.y;
    return x_diff * x_diff + y_diff * y_diff;
}


double Vector2::Dot(double tox, double toy) const {
    return this->x * tox + this->y * toy;
}


double Vector2::Dot(const Vector2 &other) const {
    return x * other.x + y * other.y;
}


double Vector2::Determinant(double ofx, double ofy) const {
    return this->x * ofy - this->y * ofx;
}


double Vector2::Determinant(const Vector2 &other) const {
    return x * other.y - y * other.x;
}


double Vector2::Angle(double tox, double toy) const {
    double dot = this->x * tox + this->y * toy;
    double det = this->x * toy - this->y * tox;
    return std::abs(std::atan2(det, dot));
}


double Vector2::Angle(const Vector2 &other) const {
    double dot = x * other.x + y * other.y;
    double det = x * other.y - y * other.x;
    return std::abs(std::atan2(det, dot));
}


Vector2 Vector2::Normalized() const {
    double length = std::sqrt(x * x + y * y);
    return Vector2(x / length, y / length);
}


Vector2 Vector2::Perpendicular() const {
    return Vector2(y, -x);
}


Vector2 Vector2::Min(double ofx, double ofy) const {
    return Vector2(this->x > ofx ? ofx : this->x, this->y > ofy ? ofy : this->y);
}


Vector2 Vector2::Min(const Vector2 &other) const {
    return Vector2(x > other.x ? other.x : x, y > other.y ? other.y : y);
}


Vector2 Vector2::Max(double ofx, double ofy) const {
    return Vector2(this->x > ofx ? this->x : ofx, this->y > ofy ? this->y : ofy);
}


Vector2 Vector2::Max(const Vector2 &other) const {
    return Vector2(x > other.x ? x : other.x, y > other.y ? y : other.y);
}


Vector2 Vector2::Lerp(double tox, double y, double t) const {
    return Vector2((double) (this->x + (tox - this->x) * t), (double) (this->y + (y - this->y) * t));
}


Vector2 Vector2::Lerp(const Vector2 &other, double t) const {
    return Vector2((double) (x + (other.x - x) * t), (double) (y + (other.y - y) * t));
}


Vector2 Vector2::operator+(const Vector2 &other) const {
    return Vector2(x + other.x, y + other.y);
}


Vector2 Vector2::operator-(const Vector2 &other) const {
    return Vector2(x - other.x, y - other.y);
}


Vector2 Vector2::operator*(const Vector2 &other) const {
    return Vector2(x * other.x, y * other.y);
}


Vector2 Vector2::operator/(const Vector2 &other) const {
    return Vector2(x / other.x, y / other.y);
}


Vector2 Vector2::operator+(double scalar) const {
    return Vector2(x + scalar, y + scalar);
}


Vector2 Vector2::operator-(double scalar) const {
    return Vector2(x - scalar, y - scalar);
}


Vector2 Vector2::operator*(double scalar) const {
    return Vector2(x * scalar, y * scalar);
}


Vector2 Vector2::operator/(double scalar) const {
    return Vector2(x / scalar, y / scalar);
}


bool Vector2::operator==(const Vector2 &other) const {
    return (x == other.x && y == other.y);
}


bool Vector2::operator==(double scalar) const {
    return (x == scalar && y == scalar);
}


Vector2 &Vector2::Normalize() {
    double length = std::sqrt(x * x + y * y);
    x /= length;
    y /= length;
    return *this;
}


Vector2 &Vector2::Zero() {
    x = 0;
    y = 0;
    return *this;
}


Vector2 &Vector2::Floor() {
    x = std::floor(x);
    y = std::floor(y);
    return *this;
}


Vector2 &Vector2::Ceil() {
    x = std::ceil(x);
    y = std::ceil(y);
    return *this;
}

Vector2 &Vector2::Round() {
    x = std::round(x);
    y = std::round(y);
    return *this;
}

Vector2 &Vector2::Set(double x, double y) {
    this->x = x;
    this->y = y;
    return *this;
}

Vector2 &Vector2::Set(const Vector2 &other) {
    x = (double) other.x;
    y = (double) other.y;
    return *this;
}

Vector2 &Vector2::Set(double scalar) {
    x = scalar;
    y = scalar;
    return *this;
}

Vector2 &Vector2::operator-() {
    x = -x;
    y = -y;
    return *this;
}

Vector2 &Vector2::operator=(const Vector2 &other) {
    x = other.x;
    y = other.y;
    return *this;
}

Vector2 &Vector2::operator+=(const Vector2 &other) {
    x += other.x;
    y += other.y;
    return *this;
}

Vector2 &Vector2::operator-=(const Vector2 &other) {
    x -= other.x;
    y -= other.y;
    return *this;
}

Vector2 &Vector2::operator*=(const Vector2 &other) {
    x *= other.x;
    y *= other.y;
    return *this;
}

Vector2 &Vector2::operator/=(const Vector2 &other) {
    x /= other.x;
    y /= other.y;
    return *this;
}

Vector2 &Vector2::operator=(double scalar) {
    x = scalar;
    y = scalar;
    return *this;
}

Vector2 &Vector2::operator+=(double scalar) {
    x += scalar;
    y += scalar;
    return *this;
}

Vector2 &Vector2::operator-=(double scalar) {
    x -= scalar;
    y -= scalar;
    return *this;
}

Vector2 &Vector2::operator*=(double scalar) {
    x *= scalar;
    y *= scalar;
    return *this;
}

Vector2 &Vector2::operator/=(double scalar) {
    x /= scalar;
    y /= scalar;
    return *this;
}


// ----------------------------------------------- Vector3 ------------------------------------------------





Vector3::Vector3() : x(0), y(0), z(0) {}


Vector3::Vector3(double scalar) : x(scalar), y(scalar), z(scalar) {}


Vector3::Vector3(double x, double y, double z) : x(x), y(y), z(z) {}


Vector3::Vector3(const Vector3 &other) : x((double) other.x), y((double) other.y), z((double) other.z) {}


Vector3::Vector3(const Vector2 &v2, double z) : x((double) v2.x), y((double) v2.y), z(z) {}


double Vector3::Length() const {
    return std::sqrt(x * x + y * y + z * z);
}


double Vector3::LengthSquared() const {
    return x * x + y * y + z * z;
}


double Vector3::Distance(const Vector3 &other) const {
    double x_diff = x - other.x;
    double y_diff = y - other.y;
    double z_diff = z - other.z;
    return std::sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
}


double Vector3::Distance(double x, double y, double z) const {
    double x_diff = this->x - x;
    double y_diff = this->y - y;
    double z_diff = this->z - z;
    return std::sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
}


double Vector3::DistanceSquared(const Vector3 &other) const {
    double x_diff = x - other.x;
    double y_diff = y - other.y;
    double z_diff = z - other.z;
    return x_diff * x_diff + y_diff * y_diff + z_diff * z_diff;
}


double Vector3::Dot(const Vector3 &other) const {
    return x * other.x + y * other.y + z * other.z;
}


double Vector3::AngleCos(const Vector3 &other) const {
    return (double) ((x * other.x + y * other.y + z * other.z) /
                     std::sqrt((x * x + y * y + z * z) * (other.x * other.x + other.y * other.y + other.z * other.z)));
}


double Vector3::Angle(const Vector3 &other) const {
    double cos = (x * other.x + y * other.y + z * other.z) /
                 std::sqrt((x * x + y * y + z * z) * (other.x * other.x + other.y * other.y + other.z * other.z));
    // doublehis is because sometimes cos goes above 1 or below -1 because of lost precision
    cos = cos < 1 ? cos : 1;
    cos = cos > -1 ? cos : -1;
    return std::acos(cos);
}


Vector3 Vector3::Normalized() const {
    double length = std::sqrt(x * x + y * y + z * z);
    return Vector3(x / length, y / length, z / length);
}


Vector3 Vector3::Cross(const Vector3 &other) const {
    return Vector3(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
}


Vector3 Vector3::Min(const Vector3 &other) const {
    return Vector3(x > other.x ? other.x : x, y > other.y ? other.y : y, z > other.z ? other.z : z);
}


Vector3 Vector3::Max(const Vector3 &other) const {
    return Vector3(x > other.x ? x : other.x, y > other.y ? y : other.y, z > other.z ? z : other.z);
}


Vector3 Vector3::Lerp(const Vector3 &other, double t) const {
    return Vector3((double) (x + (other.x - x) * t), (double) (y + (other.y - y) * t),
                   (double) (z + (other.z - z) * t));
}


Vector3 Vector3::operator+(const Vector3 &other) const {
    return Vector3(x + other.x, y + other.y, z + other.z);
}


Vector3 Vector3::operator-(const Vector3 &other) const {
    return Vector3(x - other.x, y - other.y, z - other.z);
}


Vector3 Vector3::operator*(const Vector3 &other) const {
    return Vector3(x * other.x, y * other.y, z * other.z);
}


Vector3 Vector3::operator/(const Vector3 &other) const {
    return Vector3(x / other.x, y / other.y, z / other.z);
}


Vector3 Vector3::operator+(double scalar) const {
    return Vector3(x + scalar, y + scalar, z + scalar);
}


Vector3 Vector3::operator-(double scalar) const {
    return Vector3(x - scalar, y - scalar, z - scalar);
}


Vector3 Vector3::operator*(double scalar) const {
    return Vector3(x * scalar, y * scalar, z * scalar);
}


Vector3 Vector3::operator/(double scalar) const {
    return Vector3(x / scalar, y / scalar, z / scalar);
}


bool Vector3::operator==(const Vector3 &other) const {
    return (x == other.x && y == other.y && z == other.z);
}


bool Vector3::operator==(double scalar) const {
    return (x == scalar && y == scalar && z == scalar);
}


Vector3 &Vector3::Normalize() {
    double length = std::sqrt(x * x + y * y + z * z);
    x /= length;
    y /= length;
    z /= length;
    return *this;
}


Vector3 &Vector3::Zero() {
    x = 0;
    y = 0;
    z = 0;
    return *this;
}


Vector3 &Vector3::Floor() {
    x = std::floor(x);
    y = std::floor(y);
    z = std::floor(z);
    return *this;
}


Vector3 &Vector3::Ceil() {
    x = std::ceil(x);
    y = std::ceil(y);
    z = std::ceil(z);
    return *this;
}


Vector3 &Vector3::Round() {
    x = std::round(x);
    y = std::round(y);
    z = std::round(z);
    return *this;
}


Vector3 &Vector3::Set(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
    return *this;
}


Vector3 &Vector3::Set(const Vector3 &other) {
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    return *this;
}


Vector3 &Vector3::Set(const Vector2 &other, double z) {
    this->x = other.x;
    this->y = other.y;
    this->z = z;
    return *this;
}


Vector3 &Vector3::Set(double scalar) {
    x = scalar;
    y = scalar;
    z = scalar;
    return *this;
}


Vector3 &Vector3::operator-() {
    x = -x;
    y = -y;
    z = -z;
    return *this;
}


Vector3 &Vector3::operator=(const Vector3 &other) {
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
}


Vector3 &Vector3::operator+=(const Vector3 &other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}


Vector3 &Vector3::operator-=(const Vector3 &other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
}


Vector3 &Vector3::operator*=(const Vector3 &other) {
    x *= other.x;
    y *= other.y;
    z *= other.z;
    return *this;
}


Vector3 &Vector3::operator/=(const Vector3 &other) {
    x /= other.x;
    y /= other.y;
    z /= other.z;
    return *this;
}


Vector3 &Vector3::operator=(double scalar) {
    x = scalar;
    y = scalar;
    z = scalar;
    return *this;
}


Vector3 &Vector3::operator+=(double scalar) {
    x += scalar;
    y += scalar;
    z += scalar;
    return *this;
}


Vector3 &Vector3::operator-=(double scalar) {
    x -= scalar;
    y -= scalar;
    z -= scalar;
    return *this;
}


Vector3 &Vector3::operator*=(double scalar) {
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
}


Vector3 &Vector3::operator/=(double scalar) {
    x /= scalar;
    y /= scalar;
    z /= scalar;
    return *this;
}

// ----------------------------------------------- Vector4 ------------------------------------------------





Vector4::Vector4() : x(0.0), y(0.0), z(0.0), w(0.0) {}


Vector4::Vector4(double scalar) : x(scalar), y(scalar), z(scalar), w(scalar) {}


Vector4::Vector4(double x, double y, double z, double w) : x(x), y(y), z(z), w(w) {}


Vector4::Vector4(const Vector4 &other) : x((double) other.x), y((double) other.y), z((double) other.z),
                                         w((double) other.w) {}


Vector4::Vector4(const Vector3 &v3, double w) : x((double) v3.x), y((double) v3.y), z((double) v3.z), w(w) {}


Vector4::Vector4(const Vector2 &v2, double z, double w) : x((double) v2.x), y((double) v2.y), z(z), w(w) {}


double Vector4::Length() const {
    return std::sqrt(x * x + y * y + z * z + w * w);
}


double Vector4::LengthSquared() const {
    return x * x + y * y + z * z + w * w;
}


double Vector4::Distance(const Vector4 &other) const {
    double x_diff = x - other.x;
    double y_diff = y - other.y;
    double z_diff = z - other.z;
    double w_diff = w - other.w;
    return std::sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff + w_diff * w_diff);
}


double Vector4::DistanceSquared(const Vector4 &other) const {
    double x_diff = x - other.x;
    double y_diff = y - other.y;
    double z_diff = z - other.z;
    double w_diff = w - other.w;
    return x_diff * x_diff + y_diff * y_diff + z_diff * z_diff + w_diff * w_diff;
}


double Vector4::Dot(const Vector4 &other) const {
    return x * other.x + y * other.y + z * other.z + w * other.w;
}


double Vector4::AngleCos(const Vector4 &other) const {
    return (double) ((x * other.x + y * other.y + z * other.z + w * other.w) /
                     std::sqrt((x * x + y * y + z * z + w * w) *
                               (other.x * other.x + other.y * other.y + other.z * other.z + other.w * other.w)));
}


double Vector4::Angle(const Vector4 &other) const {
    double cos = (x * other.x + y * other.y + z * other.z + w * other.w) /
                 std::sqrt((x * x + y * y + z * z + w * w) *
                           (other.x * other.x + other.y * other.y + other.z * other.z + other.w * other.w));
    // doublehis is because sometimes cos goes above 1 or below -1 because of lost precision
    cos = cos < 1 ? cos : 1;
    cos = cos > -1 ? cos : -1;
    return std::acos(cos);
}


Vector4 Vector4::Normalized() const {
    double length = std::sqrt(x * x + y * y + z * z + w * w);
    return Vector4(x / length, y / length, z / length, w / length);
}


Vector4 Vector4::Min(const Vector4 &other) const {
    return Vector4(x > other.x ? other.x : x, y > other.y ? other.y : y, z > other.z ? other.z : z,
                   w > other.w ? other.w : w);
}


Vector4 Vector4::Max(const Vector4 &other) const {
    return Vector4(x > other.x ? x : other.x, y > other.y ? y : other.y, z > other.z ? z : other.z,
                   w > other.w ? w : other.w);
}


Vector4 Vector4::Lerp(const Vector4 &other, double t) const {
    return Vector4((double) (x + (other.x - x) * t), (double) (y + (other.y - y) * t),
                   (double) (z + (other.z - z) * t),
                   (double) (w + (other.w - w) * t));
}


Vector4 Vector4::operator+(const Vector4 &other) const {
    return Vector4(x + other.x, y + other.y, z + other.z, w + other.w);
}


Vector4 Vector4::operator-(const Vector4 &other) const {
    return Vector4(x - other.x, y - other.y, z - other.z, w - other.w);
}


Vector4 Vector4::operator*(const Vector4 &other) const {
    return Vector4(x * other.x, y * other.y, z * other.z, w * other.w);
}


Vector4 Vector4::operator/(const Vector4 &other) const {
    return Vector4(x / other.x, y / other.y, z / other.z, w / other.w);
}


Vector4 Vector4::operator+(double scalar) const {
    return Vector4(x + scalar, y + scalar, z + scalar, w + scalar);
}


Vector4 Vector4::operator-(double scalar) const {
    return Vector4(x - scalar, y - scalar, z - scalar, w - scalar);
}


Vector4 Vector4::operator*(double scalar) const {
    return Vector4(x * scalar, y * scalar, z * scalar, w * scalar);
}


Vector4 Vector4::operator/(double scalar) const {
    return Vector4(x / scalar, y / scalar, z / scalar, w / scalar);
}


bool Vector4::operator==(const Vector4 &other) const {
    return (x == other.x && y == other.y && z == other.z && w == other.w);
}


bool Vector4::operator==(double scalar) const {
    return (x == scalar && y == scalar && z == scalar && w == scalar);
}


Vector4 &Vector4::Normalize() {
    double length = std::sqrt(x * x + y * y + z * z + w * w);
    x /= length;
    y /= length;
    z /= length;
    w /= length;
    return *this;
}


Vector4 &Vector4::Zero() {
    x = 0;
    y = 0;
    z = 0;
    w = 0;
    return *this;
}


Vector4 &Vector4::Floor() {
    x = std::floor(x);
    y = std::floor(y);
    z = std::floor(z);
    w = std::floor(w);
    return *this;
}


Vector4 &Vector4::Ceil() {
    x = std::ceil(x);
    y = std::ceil(y);
    z = std::ceil(z);
    w = std::ceil(w);
    return *this;
}


Vector4 &Vector4::Round() {
    x = std::round(x);
    y = std::round(y);
    z = std::round(z);
    w = std::round(w);
    return *this;
}


Vector4 &Vector4::Set(double x, double y, double z, double w) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = w;
    return *this;
}


Vector4 &Vector4::Set(const Vector4 &other) {
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    this->w = other.w;
    return *this;
}


Vector4 &Vector4::Set(const Vector3 &other, double w) {
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    this->w = w;
    return *this;
}


Vector4 &Vector4::Set(const Vector2 &other, double z, double w) {
    this->x = other.x;
    this->y = other.y;
    this->z = z;
    this->w = w;
    return *this;
}


Vector4 &Vector4::Set(double scalar) {
    x = scalar;
    y = scalar;
    z = scalar;
    w = scalar;
    return *this;
}


Vector4 &Vector4::operator-() {
    x = -x;
    y = -y;
    z = -z;
    w = -w;
    return *this;
}


Vector4 &Vector4::operator=(const Vector4 &other) {
    x = other.x;
    y = other.y;
    z = other.z;
    w = other.w;
    return *this;
}


Vector4 &Vector4::operator+=(const Vector4 &other) {
    x += other.x;
    y += other.y;
    z += other.z;
    w += other.w;
    return *this;
}


Vector4 &Vector4::operator-=(const Vector4 &other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    w -= other.w;
    return *this;
}


Vector4 &Vector4::operator*=(const Vector4 &other) {
    x *= other.x;
    y *= other.y;
    z *= other.z;
    w *= other.w;
    return *this;
}


Vector4 &Vector4::operator/=(const Vector4 &other) {
    x /= other.x;
    y /= other.y;
    z /= other.z;
    w /= other.w;
    return *this;
}


Vector4 &Vector4::operator=(double scalar) {
    x = scalar;
    y = scalar;
    z = scalar;
    w = scalar;
    return *this;
}


Vector4 &Vector4::operator+=(double scalar) {
    x += scalar;
    y += scalar;
    z += scalar;
    w += scalar;
    return *this;
}


Vector4 &Vector4::operator-=(double scalar) {
    x -= scalar;
    y -= scalar;
    z -= scalar;
    w -= scalar;
    return *this;
}


Vector4 &Vector4::operator*=(double scalar) {
    x *= scalar;
    y *= scalar;
    z *= scalar;
    w *= scalar;
    return *this;
}


Vector4 &Vector4::operator/=(double scalar) {
    x /= scalar;
    y /= scalar;
    z /= scalar;
    w /= scalar;
    return *this;
}

} //namespace math
} //namespace sw