#pragma once

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/Geometry>

namespace embree
{
    namespace isa
    {
        namespace compressed
        {


            // morton decoder code :
            //https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/
            // "Insert" a 0 bit after each of the 16 low bits of x
            //
            inline uint32_t Part1By1(uint32_t x) {
                x &= 0x0000ffff;                  // x = ---- ---- ---- ---- fedc ba98 7654 3210
                x = (x ^ (x <<  8)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
                x = (x ^ (x <<  4)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
                x = (x ^ (x <<  2)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
                x = (x ^ (x <<  1)) & 0x55555555; // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
                return x;
            }

            inline uint32_t EncodeMorton2(uint32_t x, uint32_t y) {
                return (Part1By1(y) << 1) + Part1By1(x);
            }

            inline uint32_t Compact1By1(uint32_t x) {
                x &= 0x55555555;                  // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
                x = (x ^ (x >>  1)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
                x = (x ^ (x >>  2)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
                x = (x ^ (x >>  4)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
                x = (x ^ (x >>  8)) & 0x0000ffff; // x = ---- ---- ---- ---- fedc ba98 7654 3210
                return x;

            }

            inline uint32_t DecodeMorton2X(uint32_t code)
            {
                return Compact1By1(code >> 0);
            }

            inline uint32_t DecodeMorton2Y(uint32_t code)
            {
                return Compact1By1(code >> 1);
            }


            //  XY projection
            inline Eigen::Matrix3f ComputeLinearEstimate(const std::vector<Eigen::Vector2f>& source, const std::vector<Eigen::Vector2f>& target) {
                Eigen::Matrix3f M;
                M.setZero();


                Eigen::MatrixXd A(8,8);

                for (size_t i = 0; i < 4; ++i) {

                    const Eigen::Vector2f& p = target[i];
                    const Eigen::Vector2f& q = source[i];

                    A.row(0+i) << q[0], q[1], 1.f, 0.f, 0.f, 0.f, -q[0]*p[0], -q[1]*p[0];
                    A.row(4+i) << 0.f, 0.f, 0.f, q[0], q[1], 1.f, -q[0]*p[1], -q[1]*p[1];
                }

                Eigen::VectorXd b(8);
                for (size_t i = 0; i < 2; ++i)
                    for (size_t k = 0; k < 4; ++k)
                        b.row(i*4+k) << target[k][i];

                Eigen::VectorXd x(8);

                x = A.fullPivLu().solve(b);

                M <<  x[0],  x[1], x[2],
                      x[3],  x[4], x[5],
                      x[6],  x[7],  1.f;

                return M;
            }

            inline Vec3f project(const Vec3f& a, const Eigen::Matrix3f& proj) {
                Eigen::Vector3f tmp = Eigen::Vector3f(a.x, a.y, 1.f);
                tmp = proj * tmp;
                return Vec3f(tmp.x() / tmp.z(), tmp.y() / tmp.z(), a.z);
            }

            // 2d line intersection
            inline  float intersect_line(const Vec2f& p2,
                                         const Vec2f& p3,
                                         const Vec3f& o,
                                         const Vec3f& d) {

                const Vec2f v = Vec2f(p2.x - o.x, p2.y - o.y);
                const Vec2f l = Vec2f(p3.x - p2.x, p3.y - p2.y);


                const float t1 = (l.y * v.x - l.x * v.y) / (l.y * d.x - l.x * d.y);
                const float t2 = (d.x * v.y - d.y * v.x) / (l.x * d.y - l.y * d.x);

                return t2 >= 0.f && t2 <= 1.f ? t1 : std::nanf("");
            }

            // 3d frustum box intersection
            inline bool intersect_frustum(const float* box, const Vec3f& ray_origin, const Vec3f& ray_dir, float &t_near, float &t_far) {

                const Vec3f ray_rdir = rcp_safe(ray_dir);
                const Vec3f ray_org_rdir = ray_origin*ray_rdir;

                const float t1z = box[0] * ray_rdir.z - ray_org_rdir.z;
                const float t2z = box[1] * ray_rdir.z - ray_org_rdir.z;

                const float t1x = intersect_line(*((Vec2f*)(&box[2])), *((Vec2f*)(&box[6])), ray_origin, ray_dir);
                const float t2x = intersect_line(*((Vec2f*)(&box[4])), *((Vec2f*)(&box[8])), ray_origin, ray_dir);
                const float t1y = intersect_line(*((Vec2f*)(&box[2])), *((Vec2f*)(&box[4])), ray_origin, ray_dir);
                const float t2y = intersect_line(*((Vec2f*)(&box[6])), *((Vec2f*)(&box[8])), ray_origin, ray_dir);

                // get valid frustum points
                const float near1 = fminf(fminf(t1x, t2x), fminf(t1y, t2y));
                const float far1  = fmaxf(fmaxf(t1x, t2x), fmaxf(t1y, t2y));


                // intersect frustum box
                t_near = fmaxf(fmaxf(fminf(t1z, t2z), near1), t_near);
                t_far = fminf(fminf(fmaxf(t1z, t2z), far1), t_far);

                return t_near <= t_far && near1 == near1 && far1 == far1;

            }

            inline bool intersect_patch(const unsigned idx,
                                        const float rcp_edges,
                                        const float dz,
                                        const float t1,
                                        const float t2,
                                        const float v0,
                                        const float v1,
                                        const float v2,
                                        const float v3,
                                        const BBox3fa& box,
                                        const Vec3f org,
                                        const Vec3f dir,
                                        float& u,
                                        float& v,
                                        float& tt) {



                const Vec3f p = t1 * dir + org;
                const Vec3f p2 = t2 * dir + org;

                const float lenX = rcp(box.upper.x - box.lower.x);
                const float lenY = rcp(box.upper.y - box.lower.y);

                const float fx1 = (p.x - box.lower.x) * lenX;
                const float fy1 = (p.y - box.lower.y) * lenY;
                

                // avoid too small patches
                if (t2 - t1 < 1.0E-6f) {
                    tt = t1;
                    u = (fx1 + static_cast<float>(DecodeMorton2X(idx))) * rcp_edges;
                    v = (fy1 + static_cast<float>(DecodeMorton2Y(idx))) * rcp_edges;
                    return true;
                }

                const float fx2 = (p2.x - box.lower.x) * lenX;
                const float fy2 = (p2.y - box.lower.y) * lenY;

                const float dx1 = 1.f - fx1;
                const float dy1 = 1.f - fy1;

                // entry point z-value on surface
                float z1 = v0 * dx1 * dy1 +
                           v1 * fx1 * dy1 +
                           v2 * dx1 * fy1 +
                           v3 * fx1 * fy1;

                const float dx2 = 1.f - fx2;
                const float dy2 = 1.f - fy2;

                // exit point z-value on surface
                float z2 = v0 * dx2 * dy2 +
                           v1 * fx2 * dy2 +
                           v2 * dx2 * fy2 +
                           v3 * fx2 * fy2;

                // entry point between lower and upper bound -> box hit
                if (p.z >= z1 && p.z <= z1 + dz) {
                    tt = t1;
                    u = (fx1 + static_cast<float>(DecodeMorton2X(idx))) * rcp_edges;
                    v = (fy1 + static_cast<float>(DecodeMorton2Y(idx))) * rcp_edges;
                    return true;
                }

                // entry point above upper bound -> shift z1 and z2
                if (p.z > z1 + dz) {
                    z1 += dz;
                    z2 += dz;
                }


                // 2d line intersection
                
                const float alpha = p2.z - z2;
                const float beta  = z1   - p.z;

                const float t = (t1 * alpha + t2 * beta) / (alpha + beta);

                const float d = (t - t1) / (t2 - t1);

                const float fx = fx2 - fx1;
                const float fy = fy2 - fy1;

                // test for valid hit (only if inside bbox)  and get uv-coords
                if (t < tt && t >= t1 && t <= t2) {
                    u = (fx*d + fx1 + static_cast<float>(DecodeMorton2X(idx))) * rcp_edges;
                    v = (fy*d + fy1 + static_cast<float>(DecodeMorton2Y(idx))) * rcp_edges;
                    tt = t;

                    return true;
                }

                return false;
            }

            // intersect single triangle
            inline bool intersect_triangle(const Vec3f& v0, const Vec3f& v1, const Vec3f& v2, RayHit &ray) {
                const float &a_x = v0.x;
                const float &a_y = v0.y;
                const float &a_z = v0.z;
                const float &a = a_x - v1.x;
                const float &b = a_y - v1.y;
                const float &c = a_z - v1.z;
                const float &d = a_x - v2.x;
                const float &e = a_y - v2.y;
                const float &f = a_z - v2.z;
                const float &g = ray.dir.x;
                const float &h = ray.dir.y;
                const float &i = ray.dir.z;
                const float &j = a_x - ray.org.x;
                const float &k = a_y - ray.org.y;
                const float &l = a_z - ray.org.z;

                float common1 = e*i - h*f;
                float common2 = g*f - d*i;
                float common3 = d*h - e*g;
                const float M   = rcp(a * common1  +  b * common2  +  c * common3);
                float beta  = j * common1  +  k * common2  +  l * common3;

                common1 = a*k - j*b;
                common2 = j*c - a*l;
                common3 = b*l - k*c;
                float gamma = i * common1  +  h * common2  +  g * common3;
                float tt    = -(f * common1  +  e * common2  +  d * common3);

                beta *= M;
                gamma *= M;
                tt *= M;

                if (tt > 0 && tt < ray.tfar && tt >= ray.tnear())
                    if (beta > 0 && gamma > 0 && beta + gamma <= 1)
                    {
                        ray.tfar = tt;
                        ray.u = beta;
                        ray.v = gamma;
                        return true;
                    }

                return false;
            }

            // intersects triangles built on-the-fly from grid vertices starting with idx
            inline bool intersect_triangles(const int idx, const float rcp_edges, const unsigned width, const Vec3f* grid, RayHit& ray)  {
                const unsigned x = DecodeMorton2X(idx);
                const unsigned y = DecodeMorton2Y(idx);

                const unsigned gridIdx0 = y*width+x;
                const unsigned gridIdx1 = gridIdx0+ 1;
                const unsigned gridIdx2 = gridIdx0+ width;
                const unsigned gridIdx3 = gridIdx2+ 1;

                const Vec3f& v0 = grid[gridIdx0];
                const Vec3f& v1 = grid[gridIdx1];
                const Vec3f& v2 = grid[gridIdx2];
                const Vec3f& v3 = grid[gridIdx3];

                const bool hit1 = intersect_triangle(v0, v1, v2, ray);
                const bool hit2 = intersect_triangle(v3, v2, v1, ray);

                if (hit2) {
                    ray.u = (static_cast<float>(x) + (1.f - ray.u)) * rcp_edges;
                    ray.v = (static_cast<float>(y) + (1.f - ray.v)) * rcp_edges;
                    return true;
                }

                if (hit1) {
                    ray.u = (static_cast<float>(x) + ray.u) * rcp_edges;
                    ray.v = (static_cast<float>(y) + ray.v) * rcp_edges;
                    return true;
                }

                return false;
            }
        }
    }
}


// vim: set foldmethod=marker :
