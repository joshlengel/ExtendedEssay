#pragma once

class Vec3;

void LambertDIzzo(const Vec3 &r1, const Vec3 &r2, double tof, double mu, Vec3 &v1, Vec3 &v2);
void LambertUV(const Vec3 &r1, const Vec3 &r2, double tof, double mu, Vec3 &v1, Vec3 &v2);