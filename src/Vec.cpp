#include"Vec.h"

#include<cmath>

Vec3::Vec3(double x, double y, double z): x(x), y(y), z(z) {}
Vec3::Vec3(): x(0.0), y(0.0), z(0.0) {}

Vec3 Vec3::operator+(const Vec3 &v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
Vec3 Vec3::operator-(const Vec3 &v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
Vec3 Vec3::operator*(double f) const { return Vec3(x * f, y * f, z * f); }
Vec3 Vec3::operator/(double f) const { return Vec3(x / f, y / f, z / f); }
Vec3 Vec3::operator-() const { return Vec3(-x, -y, -z); }

Vec3 &Vec3::operator+=(const Vec3 &v) { return *this = *this + v; }
Vec3 &Vec3::operator-=(const Vec3 &v) { return *this = *this - v; }
Vec3 &Vec3::operator*=(double f) { return *this = *this * f; }
Vec3 &Vec3::operator/=(double f) { return *this = *this / f; }

double Vec3::Dot(const Vec3 &v1, const Vec3 &v2) { return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z; }
Vec3 Vec3::Cross(const Vec3 &v1, const Vec3 &v2) { return Vec3(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x); }
double Vec3::Length(const Vec3 &v) { return sqrt(Dot(v, v)); }
Vec3 Vec3::Normalize(const Vec3 &v) { return v / Length(v); }

Vec3 Vec3::Rotate(const Vec3 &v, const Vec3 &axis, double angle) { return axis * Vec3::Dot(axis, v) + Vec3::Cross(Vec3::Cross(axis, v), axis) * cos(angle) + Vec3::Cross(axis, v) * sin(angle); }