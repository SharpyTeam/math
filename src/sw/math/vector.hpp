//
// Created by Ilya on 30.06.2019.
//

#ifndef SANDWICH_MATH_VECTOR_HPP
#define SANDWICH_MATH_VECTOR_HPP

namespace sw {
namespace math {

class Vector2 {
public:
    double x;
    double y;

    Vector2();
    explicit Vector2(double scalar);
    Vector2(double x, double y);
    Vector2(const Vector2 &other);

    double Angle(double x, double y) const;
    double Angle(const Vector2 &other) const;
    double Determinant(double x, double y) const;
    double Determinant(const Vector2 &other) const;
    double Distance(double x, double y) const;
    double Distance(const Vector2 &other) const;
    double DistanceSquared(double x, double y) const;
    double DistanceSquared(const Vector2 &other) const;
    double Dot(double x, double y) const;
    double Dot(const Vector2 &other) const;
    double Length() const;
    double LengthSquared() const;

    Vector2 Lerp(double x, double y, double t) const;
    Vector2 Lerp(const Vector2 &other, double t) const;
    Vector2 Max(double x, double y) const;
    Vector2 Max(const Vector2 &other) const;
    Vector2 Min(double x, double y) const;
    Vector2 Min(const Vector2 &other) const;
    Vector2 Normalized() const;
    Vector2 Perpendicular() const;

    Vector2 operator-() const;
    Vector2 operator+(double scalar) const;
    Vector2 operator-(double scalar) const;
    Vector2 operator*(double scalar) const;
    Vector2 operator/(double scalar) const;
    Vector2 operator+(const Vector2 &other) const;
    Vector2 operator-(const Vector2 &other) const;
    Vector2 operator*(const Vector2 &other) const;
    Vector2 operator/(const Vector2 &other) const;

    bool operator==(double scalar) const;
    bool operator==(const Vector2 &other) const;

    Vector2 &Ceil();
    Vector2 &Floor();
    Vector2 &Normalize();
    Vector2 &Round();
    Vector2 &Set(double scalar);
    Vector2 &Set(double x, double y);
    Vector2 &Set(const Vector2 &other);
    Vector2 &Zero();

    Vector2 &operator=(double scalar);
    Vector2 &operator+=(double scalar);
    Vector2 &operator-=(double scalar);
    Vector2 &operator*=(double scalar);
    Vector2 &operator/=(double scalar);
    Vector2 &operator=(const Vector2 &other);
    Vector2 &operator+=(const Vector2 &other);
    Vector2 &operator-=(const Vector2 &other);
    Vector2 &operator*=(const Vector2 &other);
    Vector2 &operator/=(const Vector2 &other);
};

class Vector3 {
public:
    double x;
    double y;
    double z;

    Vector3();
    explicit Vector3(double scalar);
    Vector3(double x, double y, double z);
    Vector3(const Vector3 &other);
    Vector3(const Vector2 &other, double z);

    double Angle(double x, double y, double z) const;
    double Angle(const Vector3 &other) const;
    double AngleCos(double x, double y, double z) const;
    double AngleCos(const Vector3 &other) const;
    double Distance(double x, double y, double z) const;
    double Distance(const Vector3 &other) const;
    double DistanceSquared(double x, double y, double z) const;
    double DistanceSquared(const Vector3 &other) const;
    double Dot(double x, double y, double z) const;
    double Dot(const Vector3 &other) const;
    double Length() const;
    double LengthSquared() const;

    Vector3 Cross(double x, double y, double z) const;
    Vector3 Cross(const Vector3 &other) const;
    Vector3 Lerp(double x, double y, double z, double t) const;
    Vector3 Lerp(const Vector3 &other, double t) const;
    Vector3 Max(double x, double y, double z) const;
    Vector3 Max(const Vector3 &other) const;
    Vector3 Min(double x, double y, double z) const;
    Vector3 Min(const Vector3 &other) const;
    Vector3 Normalized() const;

    Vector3 operator-() const;
    Vector3 operator+(double scalar) const;
    Vector3 operator-(double scalar) const;
    Vector3 operator*(double scalar) const;
    Vector3 operator/(double scalar) const;
    Vector3 operator+(const Vector3 &other) const;
    Vector3 operator-(const Vector3 &other) const;
    Vector3 operator*(const Vector3 &other) const;
    Vector3 operator/(const Vector3 &other) const;

    bool operator==(double scalar) const;
    bool operator==(const Vector3 &other) const;

    Vector3 &Ceil();
    Vector3 &Floor();
    Vector3 &Normalize();
    Vector3 &Round();
    Vector3 &Set(double scalar);
    Vector3 &Set(double x, double y, double z);
    Vector3 &Set(const Vector3 &other);
    Vector3 &Set(const Vector2 &other, double z);
    Vector3 &Zero();

    Vector3 &operator=(double scalar);
    Vector3 &operator+=(double scalar);
    Vector3 &operator-=(double scalar);
    Vector3 &operator*=(double scalar);
    Vector3 &operator/=(double scalar);
    Vector3 &operator=(const Vector3 &other);
    Vector3 &operator+=(const Vector3 &other);
    Vector3 &operator-=(const Vector3 &other);
    Vector3 &operator*=(const Vector3 &other);
    Vector3 &operator/=(const Vector3 &other);
};

class Vector4 {
public:
    double x;
    double y;
    double z;
    double w;

    Vector4();
    explicit Vector4(double scalar);
    Vector4(double x, double y, double z, double w);
    Vector4(const Vector4 &other);
    Vector4(const Vector3 &other, double w);
    Vector4(const Vector2 &other, double z, double w);

    double Angle(double x, double y, double z, double w) const;
    double Angle(const Vector4 &other) const;
    double AngleCos(double x, double y, double z, double w) const;
    double AngleCos(const Vector4 &other) const;
    double Distance(double x, double y, double z, double w) const;
    double Distance(const Vector4 &other) const;
    double DistanceSquared(double x, double y, double z, double w) const;
    double DistanceSquared(const Vector4 &other) const;
    double Dot(double x, double y, double z, double w) const;
    double Dot(const Vector4 &other) const;
    double Length() const;
    double LengthSquared() const;

    Vector4 Lerp(double x, double y, double z, double w, double t) const;
    Vector4 Lerp(const Vector4 &other, double t) const;
    Vector4 Max(double x, double y, double z, double w) const;
    Vector4 Max(const Vector4 &other) const;
    Vector4 Min(double x, double y, double z, double w) const;
    Vector4 Min(const Vector4 &other) const;
    Vector4 Normalized() const;

    Vector4 operator-() const;
    Vector4 operator+(double scalar) const;
    Vector4 operator-(double scalar) const;
    Vector4 operator*(double scalar) const;
    Vector4 operator/(double scalar) const;
    Vector4 operator+(const Vector4 &other) const;
    Vector4 operator-(const Vector4 &other) const;
    Vector4 operator*(const Vector4 &other) const;
    Vector4 operator/(const Vector4 &other) const;

    bool operator==(double scalar) const;
    bool operator==(const Vector4 &other) const;

    Vector4 &Ceil();
    Vector4 &Floor();
    Vector4 &Normalize();
    Vector4 &Round();
    Vector4 &Set(double scalar);
    Vector4 &Set(double x, double y, double z, double w);
    Vector4 &Set(const Vector4 &other);
    Vector4 &Set(const Vector3 &other, double w);
    Vector4 &Set(const Vector2 &other, double z, double w);
    Vector4 &Zero();

    Vector4 &operator=(double scalar);
    Vector4 &operator+=(double scalar);
    Vector4 &operator-=(double scalar);
    Vector4 &operator*=(double scalar);
    Vector4 &operator/=(double scalar);
    Vector4 &operator=(const Vector4 &other);
    Vector4 &operator+=(const Vector4 &other);
    Vector4 &operator-=(const Vector4 &other);
    Vector4 &operator*=(const Vector4 &other);
    Vector4 &operator/=(const Vector4 &other);
};

} //namespace math
} //namespace sw

#if SANDWICH_MATH_INLINE
    #include "vector.ipp"
#endif

#endif //SANDWICH_MATH_VECTOR_HPP
