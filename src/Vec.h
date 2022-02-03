#pragma once

struct Vec3
{
    double x, y, z;

    Vec3(double x, double y, double z);
    Vec3();

    Vec3 operator+(const Vec3 &v) const;
    Vec3 operator-(const Vec3 &v) const;
    Vec3 operator*(double f) const;
    Vec3 operator/(double f) const;
    Vec3 operator-() const;

    Vec3 &operator+=(const Vec3 &v);
    Vec3 &operator-=(const Vec3 &v);
    Vec3 &operator*=(double f);
    Vec3 &operator/=(double f);

    static double Dot(const Vec3 &v1, const Vec3 &v2);
    static Vec3 Cross(const Vec3 &v1, const Vec3 &v2);
    static double Length(const Vec3 &v);
    static Vec3 Normalize(const Vec3 &v);

    static Vec3 Rotate(const Vec3 &v, const Vec3 &axis, double angle);
};