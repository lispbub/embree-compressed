#pragma once

#include <vector>
#include <cmath>

namespace embree
{
    namespace isa
    {
        namespace compressed
        {

            class leafStorageBase {
            public:
                inline void setZ(const float& z1, const float& z2, const float& z3, const float& z4, const float& zf) const {}
            };

            template<int bits>
                class leafStorage : leafStorageBase {};

            template<>
                struct leafStorage<4> : leafStorageBase {
                public:
                    unsigned char z12;
                    unsigned char z34;

                    inline void setZ1(const unsigned char z) { z12 = (z12 & 0x0f) | (z << 4); }
                    inline void setZ2(const unsigned char z) { z12 = (z12 & 0xf0) | (z & 0x0f); }
                    inline void setZ3(const unsigned char z) { z34 = (z34 & 0x0f) | (z << 4); }
                    inline void setZ4(const unsigned char z) { z34 = (z34 & 0xf0) | (z & 0x0f); }

                    inline void setZero() {
                        setZ1(0);
                        setZ2(0);
                        setZ3(0);
                        setZ4(0);
                    }


                    inline void setZ(const float& z1, const float& z2, const float& z3, const float& z4, const float& zf) {
                        const float rcpF = 16.f / zf;

                        setZ1(static_cast<unsigned char>(fmaxf(0.f, fminf(15.f, z1 * rcpF))));
                        setZ2(static_cast<unsigned char>(fmaxf(0.f, fminf(15.f, z2 * rcpF))));
                        setZ3(static_cast<unsigned char>(fmaxf(0.f, fminf(15.f, z3 * rcpF))));
                        setZ4(static_cast<unsigned char>(fmaxf(0.f, fminf(15.f, z4 * rcpF))));
                    }

                    inline void getZ(Vec3f& v1,
                                     Vec3f& v2,
                                     Vec3f& v3,
                                     Vec3f& v4,
                                     const BBox3f& pBox)  const  {
                        const float rcpF = rcp(16.f) * (pBox.upper.z - pBox.lower.z);

                        v1.x = v3.x = pBox.lower.x;
                        v2.x = v4.x = pBox.upper.x;

                        v1.y = v2.y = pBox.lower.y;
                        v3.y = v4.y = pBox.upper.y;

                        v1.z = pBox.lower.z + rcpF * static_cast<float>((z12 & 0xf0) >> 4);
                        v2.z = pBox.lower.z + rcpF * static_cast<float>((z12 & 0x0f));
                        v3.z = pBox.lower.z + rcpF * static_cast<float>((z34 & 0xf0) >> 4);
                        v4.z = pBox.lower.z + rcpF * static_cast<float>((z34 & 0x0f));


                    }

                    inline void getZ(Vec3f& v1,
                                     Vec3f& v2,
                                     Vec3f& v3,
                                     Vec3f& v4,
                                     const BBox3f& pBox,
                                     const float range,
                                     const float offset)  const  {
                        const float rcpF = rcp(16.f) * range;

                        v1.x = v3.x = pBox.lower.x;
                        v2.x = v4.x = pBox.upper.x;

                        v1.y = v2.y = pBox.lower.y;
                        v3.y = v4.y = pBox.upper.y;

                        v1.z = offset + rcpF * static_cast<float>((z12 & 0xf0) >> 4);
                        v2.z = offset + rcpF * static_cast<float>((z12 & 0x0f));
                        v3.z = offset + rcpF * static_cast<float>((z34 & 0xf0) >> 4);
                        v4.z = offset + rcpF * static_cast<float>((z34 & 0x0f));


                    }

                    inline void getZ(float& v1,
                                     float& v2,
                                     float& v3,
                                     float& v4,
                                     const float range,
                                     const float offset)  const  {
                        const float rcpF = rcp(16.f) * range;

                        v1 = offset + rcpF * static_cast<float>((z12 & 0xf0) >> 4);
                        v2 = offset + rcpF * static_cast<float>((z12 & 0x0f));
                        v3 = offset + rcpF * static_cast<float>((z34 & 0xf0) >> 4);
                        v4 = offset + rcpF * static_cast<float>((z34 & 0x0f));


                    }

                    inline float getDelta() const {
                        return rcp(16.f);
                    }

                };

            class baseLeaf {
            public:
                float refitTriangle(const Vec2f& p,
                                    const Vec3f& p1,
                                    const Vec3f& p2,
                                    const Vec3f& p3) {

                    Vec3f ta; ta.x = p1.x; ta.y = p1.y; ta.z = p1.z;
                    Vec3f tb; tb.x = p2.x; tb.y = p2.y; tb.z = p2.z;
                    Vec3f tc; tc.x = p3.x; tc.y = p3.y; tc.z = p3.z;

                    Vec3f o; o.x = p.x; o.y = p.y; o.z = 0.f;
                    Vec3f d; d.x = 0.f, d.y = 0.f; d.z = 1.f;

                    float tmp = ta.x, tt = ta.y, tmp2 = ta.z;
                    const float a = tmp - tb.x;
                    const float f = tmp - tc.x;
                    const float j = tmp - o.x;

                    const float b = tt - tb.y;
                    const float e = tt - tc.y;
                    const float k = tt - o.y;

                    const float c = tmp2 - tb.z;
                    tt = tmp2 - tc.z;
                    const float l = tmp2 - o.z;

                    const float g = d.x;
                    const float h = d.y;
                    float gamma = d.z;

                    tmp = e*gamma - h*tt;
                    tmp2 = g*tt - f*gamma;
                    float M = a * tmp;

                    tmp = f*h - e*g;
                    M += b*tmp2;

                    tmp2 = a*k - j*b;
                    M += c * tmp;

                    tmp = j*c - a*l;
                    tt *= tmp2;

                    M = 1.f / M;

                    tmp2 = b*l - k*c;
                    tt +=  e * tmp;

                    tt +=  f * tmp2;

                    tt *= -M;

                    return tt;
                };

            };

            template<int bits>
                class quantTris : baseLeaf {
                public:
                    typedef leafStorage<bits> leafStorage_t;
                    leafStorage<bits> z;

                    void setZ(const BBox3f& parentBox,
                              const float& z1,
                              const float& z2,
                              const float& z3,
                              const float& z4) {

                        z.setZ(z1 - parentBox.lower.z,
                               z2 - parentBox.lower.z,
                               z3 - parentBox.lower.z,
                               z4 - parentBox.lower.z,
                               parentBox.upper.z - parentBox.lower.z);

                    }

                    void setZ(const BBox3f& parentBox,
                              const Vec3f& v1,
                              const Vec3f& v2,
                              const Vec3f& v3,
                              const Vec3f& v4,
                              const float extent) {

                        const Vec2f lowerLeft  = Vec2f(parentBox.lower.x, parentBox.lower.y);
                        const Vec2f lowerRight = Vec2f(parentBox.upper.x, parentBox.lower.y);
                        const Vec2f upperLeft  = Vec2f(parentBox.lower.x, parentBox.upper.y);
                        const Vec2f upperRight = Vec2f(parentBox.upper.x, parentBox.upper.y);

                        const float z1 = refitTriangle(lowerLeft,  v1, v2, v3);
                        const float z2 = refitTriangle(lowerRight, v1, v2, v4);
                        const float z3 = refitTriangle(upperLeft,  v1, v3, v4);
                        const float z4 = refitTriangle(upperRight, v2, v3, v4);

                        const double zF = parentBox.upper.z - parentBox.lower.z;
                        if (zF != 0.0) {
                            z.setZ(z1 - (parentBox.lower.z - extent*zF),
                                   z2 - (parentBox.lower.z - extent*zF),
                                   z3 - (parentBox.lower.z - extent*zF),
                                   z4 - (parentBox.lower.z - extent*zF),
                                   (1.0f + 2.0f*extent)*zF);
                        }
                        else {
                            z.setZero();
                        }
                    }

                    float estimateExtent(const BBox3f& parentBox,
                                         const Vec3f& v1,
                                         const Vec3f& v2,
                                         const Vec3f& v3,
                                         const Vec3f& v4) {

                        const Vec2f lowerLeft  = Vec2f(parentBox.lower.x, parentBox.lower.y);
                        const Vec2f lowerRight = Vec2f(parentBox.upper.x, parentBox.lower.y);
                        const Vec2f upperLeft  = Vec2f(parentBox.lower.x, parentBox.upper.y);
                        const Vec2f upperRight = Vec2f(parentBox.upper.x, parentBox.upper.y);


                        const float z1 = refitTriangle(lowerLeft,  v1, v2, v3);
                        const float z2 = refitTriangle(lowerRight, v1, v2, v4);
                        const float z3 = refitTriangle(upperLeft,  v1, v3, v4);
                        const float z4 = refitTriangle(upperRight, v2, v3, v4);

                        const float t1 = fmaxf( fmaxf(z1 - parentBox.upper.z, 0.f), fabsf( fminf(z1 - parentBox.lower.z, 0.f)));
                        const float t2 = fmaxf( fmaxf(z2 - parentBox.upper.z, 0.f), fabsf( fminf(z2 - parentBox.lower.z, 0.f)));
                        const float t3 = fmaxf( fmaxf(z3 - parentBox.upper.z, 0.f), fabsf( fminf(z3 - parentBox.lower.z, 0.f)));
                        const float t4 = fmaxf( fmaxf(z4 - parentBox.upper.z, 0.f), fabsf( fminf(z4 - parentBox.lower.z, 0.f)));

                        const double zF = parentBox.upper.z - parentBox.lower.z;

                        if (zF == 0.0)
                            return 0.f;
                        else
                            return fmaxf(fmaxf(t1, t2), fmaxf(t3, t4)) / zF;
                    }




                    inline void getZ(Vec3f& v1,
                                     Vec3f& v2,
                                     Vec3f& v3,
                                     Vec3f& v4,
                                     const BBox3f& pBox)  const  {
                        z.getZ(v1, v2, v3, v4, pBox);
                    }

                    inline void getZ(Vec3f& v1,
                                     Vec3f& v2,
                                     Vec3f& v3,
                                     Vec3f& v4,
                                     const BBox3f& pBox,
                                     const float range,
                                     const float offset)  const  {
                        z.getZ(v1, v2, v3, v4, pBox, range, offset);
                    }

                    inline void getZ(float& v1,
                                     float& v2,
                                     float& v3,
                                     float& v4,
                                     const BBox3f& pBox)  const  {
                        z.getZ(v1, v2, v3, v4, pBox);
                    }

                    inline void getZ(float& v1,
                                     float& v2,
                                     float& v3,
                                     float& v4,
                                     const float range,
                                     const float offset)  const  {
                        z.getZ(v1, v2, v3, v4, range, offset);
                    }

                };
        }
    }
}

