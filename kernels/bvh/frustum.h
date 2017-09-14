// ======================================================================== //
// Copyright 2009-2017 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#pragma once

#include "bvh.h"
#include "../common/ray.h"
#include "../common/stack_item.h"

namespace embree
{
  namespace isa
  {
    struct NearFarPreCompute
    {
      size_t nearX, nearY, nearZ;
      size_t farX, farY, farZ;

      __forceinline NearFarPreCompute() {}
      
      __forceinline NearFarPreCompute(const Vec3fa& dir, size_t N)
      {
        const size_t M = sizeof(float)*N;
        nearX = (dir.x < 0.0f) ? 1*M : 0*M;
        nearY = (dir.y < 0.0f) ? 3*M : 2*M;
        nearZ = (dir.z < 0.0f) ? 5*M : 4*M;
        farX  = nearX ^ M;
        farY  = nearY ^ M;
        farZ  = nearZ ^ M;
      }
    };

    template<int K>
    struct Packet
    {
      Vec3vf<K> rdir;
      Vec3vf<K> org_rdir;
      vfloat<K> min_dist;
      vfloat<K> max_dist;

      __forceinline Packet () {}
      
      __forceinline Packet (const Vec3vf<K>& org,
                            const Vec3vf<K>& dir,
                            const vfloat<K>& tnear,
                            const vfloat<K>& tfar,
                            bool robust)
      {
        rdir = rcp_safe(dir);
        if (robust) org_rdir = org;
        else        org_rdir = org * rdir;
        min_dist = tnear;
        max_dist = tfar;
      }

      __forceinline size_t intersectFast(const vfloat<K>& minX,
                                         const vfloat<K>& minY,
                                         const vfloat<K>& minZ,
                                         const vfloat<K>& maxX,
                                         const vfloat<K>& maxY,
                                         const vfloat<K>& maxZ) const
      {
        const vfloat<K> tminX = msub(minX, rdir.x, org_rdir.x);
        const vfloat<K> tminY = msub(minY, rdir.y, org_rdir.y);
        const vfloat<K> tminZ = msub(minZ, rdir.z, org_rdir.z);
        const vfloat<K> tmaxX = msub(maxX, rdir.x, org_rdir.x);
        const vfloat<K> tmaxY = msub(maxY, rdir.y, org_rdir.y);
        const vfloat<K> tmaxZ = msub(maxZ, rdir.z, org_rdir.z);
        const vfloat<K> tmin  = maxi(tminX, tminY, tminZ, min_dist);
        const vfloat<K> tmax  = mini(tmaxX, tmaxY, tmaxZ, max_dist);
        const vbool<K> vmask  = tmin <= tmax;
        return movemask(vmask);
      }

      __forceinline size_t intersectRobust(const vfloat<K>& minX,
                                           const vfloat<K>& minY,
                                           const vfloat<K>& minZ,
                                           const vfloat<K>& maxX,
                                           const vfloat<K>& maxY,
                                           const vfloat<K>& maxZ) const
      {
        const vfloat<K> tminX = (minX - org_rdir.x) * rdir.x;
        const vfloat<K> tminY = (minY - org_rdir.y) * rdir.y;
        const vfloat<K> tminZ = (minZ - org_rdir.z) * rdir.z;
        const vfloat<K> tmaxX = (maxX - org_rdir.x) * rdir.x;
        const vfloat<K> tmaxY = (maxY - org_rdir.y) * rdir.y;
        const vfloat<K> tmaxZ = (maxZ - org_rdir.z) * rdir.z;
        const float round_down = 1.0f-2.0f*float(ulp); // FIXME: use per instruction rounding for AVX512 
        const float round_up   = 1.0f+2.0f*float(ulp); 
        const vfloat<K> tmin  = round_down*max(tminX, tminY, tminZ, min_dist);
        const vfloat<K> tmax  = round_up  *min(tmaxX, tmaxY, tmaxZ, max_dist);
        const vbool<K> vmask  = tmin <= tmax;
        return movemask(vmask);
      }

      template<int N, int Nx>
        __forceinline size_t intersectRayFast(const typename BVHN<N>::AlignedNode* __restrict__ node, size_t rid, const NearFarPreCompute& nf) const
      {
        const vfloat<Nx> bminX = vfloat<Nx>(*(const vfloat<N>*)((const char*)&node->lower_x + nf.nearX));
        const vfloat<Nx> bminY = vfloat<Nx>(*(const vfloat<N>*)((const char*)&node->lower_x + nf.nearY));
        const vfloat<Nx> bminZ = vfloat<Nx>(*(const vfloat<N>*)((const char*)&node->lower_x + nf.nearZ));
        const vfloat<Nx> bmaxX = vfloat<Nx>(*(const vfloat<N>*)((const char*)&node->lower_x + nf.farX));
        const vfloat<Nx> bmaxY = vfloat<Nx>(*(const vfloat<N>*)((const char*)&node->lower_x + nf.farY));
        const vfloat<Nx> bmaxZ = vfloat<Nx>(*(const vfloat<N>*)((const char*)&node->lower_x + nf.farZ));
        
        const vfloat<Nx> rminX = msub(bminX, vfloat<Nx>(rdir.x[rid]), vfloat<Nx>(org_rdir.x[rid]));
        const vfloat<Nx> rminY = msub(bminY, vfloat<Nx>(rdir.y[rid]), vfloat<Nx>(org_rdir.y[rid]));
        const vfloat<Nx> rminZ = msub(bminZ, vfloat<Nx>(rdir.z[rid]), vfloat<Nx>(org_rdir.z[rid]));
        const vfloat<Nx> rmaxX = msub(bmaxX, vfloat<Nx>(rdir.x[rid]), vfloat<Nx>(org_rdir.x[rid]));
        const vfloat<Nx> rmaxY = msub(bmaxY, vfloat<Nx>(rdir.y[rid]), vfloat<Nx>(org_rdir.y[rid]));
        const vfloat<Nx> rmaxZ = msub(bmaxZ, vfloat<Nx>(rdir.z[rid]), vfloat<Nx>(org_rdir.z[rid]));
        const vfloat<Nx> rmin  = maxi(rminX, rminY, rminZ, vfloat<Nx>(min_dist[rid]));
        const vfloat<Nx> rmax  = mini(rmaxX, rmaxY, rmaxZ, vfloat<Nx>(max_dist[rid]));

        const vbool<Nx> vmask_first_hit = rmin <= rmax;

        return movemask(vmask_first_hit) & (((size_t)1 << N)-1);
      }

      template<int N, int Nx>
      __forceinline size_t intersectRayRobust(const typename BVHN<N>::AlignedNode* __restrict__ node, size_t rid, const NearFarPreCompute& nf) const
      {
        const vfloat<Nx> bminX = vfloat<Nx>(*(const vfloat<N>*)((const char*)&node->lower_x + nf.nearX));
        const vfloat<Nx> bminY = vfloat<Nx>(*(const vfloat<N>*)((const char*)&node->lower_x + nf.nearY));
        const vfloat<Nx> bminZ = vfloat<Nx>(*(const vfloat<N>*)((const char*)&node->lower_x + nf.nearZ));
        const vfloat<Nx> bmaxX = vfloat<Nx>(*(const vfloat<N>*)((const char*)&node->lower_x + nf.farX));
        const vfloat<Nx> bmaxY = vfloat<Nx>(*(const vfloat<N>*)((const char*)&node->lower_x + nf.farY));
        const vfloat<Nx> bmaxZ = vfloat<Nx>(*(const vfloat<N>*)((const char*)&node->lower_x + nf.farZ));
        
        const vfloat<Nx> rminX = (bminX - vfloat<Nx>(org_rdir.x[rid])) * vfloat<Nx>(rdir.x[rid]);
        const vfloat<Nx> rminY = (bminY - vfloat<Nx>(org_rdir.y[rid])) * vfloat<Nx>(rdir.y[rid]);
        const vfloat<Nx> rminZ = (bminZ - vfloat<Nx>(org_rdir.z[rid])) * vfloat<Nx>(rdir.z[rid]);
        const vfloat<Nx> rmaxX = (bmaxX - vfloat<Nx>(org_rdir.x[rid])) * vfloat<Nx>(rdir.x[rid]);
        const vfloat<Nx> rmaxY = (bmaxY - vfloat<Nx>(org_rdir.y[rid])) * vfloat<Nx>(rdir.y[rid]);
        const vfloat<Nx> rmaxZ = (bmaxZ - vfloat<Nx>(org_rdir.z[rid])) * vfloat<Nx>(rdir.z[rid]);
        const float round_down = 1.0f-2.0f*float(ulp); // FIXME: use per instruction rounding for AVX512
        const float round_up   = 1.0f+2.0f*float(ulp);
        const vfloat<Nx> rmin  = round_down*max(rminX, rminY, rminZ, vfloat<Nx>(min_dist[rid]));
        const vfloat<Nx> rmax  = round_up  *min(rmaxX, rmaxY, rmaxZ, vfloat<Nx>(max_dist[rid]));

        const vbool<Nx> vmask_first_hit = rmin <= rmax;

        return movemask(vmask_first_hit) & (((size_t)1 << N)-1);
      }

      template<int N>
        __forceinline size_t intersectFast(const typename BVHN<N>::AlignedNode* __restrict__ node, size_t bid, const NearFarPreCompute& nf) const
      {
        char *ptr = (char*)&node->lower_x + bid*sizeof(float);
        const vfloat<K> bminX = *(const float*)(ptr + nf.nearX);
        const vfloat<K> bminY = *(const float*)(ptr + nf.nearY);
        const vfloat<K> bminZ = *(const float*)(ptr + nf.nearZ);
        const vfloat<K> bmaxX = *(const float*)(ptr + nf.farX);
        const vfloat<K> bmaxY = *(const float*)(ptr + nf.farY);
        const vfloat<K> bmaxZ = *(const float*)(ptr + nf.farZ);

        const vfloat<K> rminX = msub(bminX, rdir.x, org_rdir.x);
        const vfloat<K> rminY = msub(bminY, rdir.y, org_rdir.y);
        const vfloat<K> rminZ = msub(bminZ, rdir.z, org_rdir.z);
        const vfloat<K> rmaxX = msub(bmaxX, rdir.x, org_rdir.x);
        const vfloat<K> rmaxY = msub(bmaxY, rdir.y, org_rdir.y);
        const vfloat<K> rmaxZ = msub(bmaxZ, rdir.z, org_rdir.z);
        
        const vfloat<K> rmin  = maxi(rminX, rminY, rminZ, min_dist);
        const vfloat<K> rmax  = mini(rmaxX, rmaxY, rmaxZ, max_dist);

        const vbool<K> vmask_first_hit = rmin <= rmax;

        return movemask(vmask_first_hit);
      }

      template<int N>
      __forceinline size_t intersectRobust(const typename BVHN<N>::AlignedNode* __restrict__ node, size_t bid, const NearFarPreCompute& nf) const
      {
        char *ptr = (char*)&node->lower_x + bid*sizeof(float);
        const vfloat<K> bminX = *(const float*)(ptr + nf.nearX);
        const vfloat<K> bminY = *(const float*)(ptr + nf.nearY);
        const vfloat<K> bminZ = *(const float*)(ptr + nf.nearZ);
        const vfloat<K> bmaxX = *(const float*)(ptr + nf.farX);
        const vfloat<K> bmaxY = *(const float*)(ptr + nf.farY);
        const vfloat<K> bmaxZ = *(const float*)(ptr + nf.farZ);
        
        const vfloat<K> rminX = (bminX - org_rdir.x) * rdir.x;
        const vfloat<K> rminY = (bminY - org_rdir.y) * rdir.y;
        const vfloat<K> rminZ = (bminZ - org_rdir.z) * rdir.z;
        const vfloat<K> rmaxX = (bmaxX - org_rdir.x) * rdir.x;
        const vfloat<K> rmaxY = (bmaxY - org_rdir.y) * rdir.y;
        const vfloat<K> rmaxZ = (bmaxZ - org_rdir.z) * rdir.z;
        
        const float round_down = 1.0f-2.0f*float(ulp); // FIXME: use per instruction rounding for AVX512
        const float round_up   = 1.0f+2.0f*float(ulp);
        const vfloat<K> rmin  = round_down*max(rminX, rminY, rminZ, vfloat<K>(min_dist));
        const vfloat<K> rmax  = round_up  *min(rmaxX, rmaxY, rmaxZ, vfloat<K>(max_dist));

        const vbool<K> vmask_first_hit = rmin <= rmax;

        return movemask(vmask_first_hit);
      }

    };
    
    /* Optimized frustum test. We calculate t=(p-org)/dir in ray/box
       * intersection. We assume the rays are split by octant, thus
       * dir intervals are either positive or negative in each
       * dimension.

         Case 1: dir.min >= 0 && dir.max >= 0:
           t_min = (p_min - org_max) / dir_max = (p_min - org_max)*rdir_min = p_min*rdir_min - org_max*rdir_min
           t_max = (p_max - org_min) / dir_min = (p_max - org_min)*rdir_max = p_max*rdir_max - org_min*rdir_max

         Case 2: dir.min < 0 && dir.max < 0:
           t_min = (p_max - org_min) / dir_min = (p_max - org_min)*rdir_max = p_max*rdir_max - org_min*rdir_max
           t_max = (p_min - org_max) / dir_max = (p_min - org_max)*rdir_min = p_min*rdir_min - org_max*rdir_min
      */

    template<int N, int Nx, int K, bool robust>
      struct Frustum
      {
        __forceinline Frustum() {}

        __forceinline Frustum(const vbool<K>& valid,
                              const Vec3vf<K>& org,
                              const Vec3vf<K>& rdir,
                              const vfloat<K>& ray_tnear,
                              const vfloat<K>& ray_tfar)
        {
          init(valid,org,rdir,ray_tnear,ray_tfar);
        }

        __forceinline void init(const vbool<K>& valid,
                                const Vec3vf<K>& org,
                                const Vec3vf<K>& rdir,
                                const vfloat<K>& ray_tnear,
                                const vfloat<K>& ray_tfar)
        {
          const Vec3fa reduced_min_org( reduce_min(select(valid,org.x,pos_inf)),
                                        reduce_min(select(valid,org.y,pos_inf)),
                                        reduce_min(select(valid,org.z,pos_inf)) );
          const Vec3fa reduced_max_org( reduce_max(select(valid,org.x,neg_inf)),
                                        reduce_max(select(valid,org.y,neg_inf)),
                                        reduce_max(select(valid,org.z,neg_inf)) );
          
          const Vec3fa reduced_min_rdir( reduce_min(select(valid,rdir.x,pos_inf)),
                                         reduce_min(select(valid,rdir.y,pos_inf)),
                                         reduce_min(select(valid,rdir.z,pos_inf)) );
          const Vec3fa reduced_max_rdir( reduce_max(select(valid,rdir.x,neg_inf)),
                                         reduce_max(select(valid,rdir.y,neg_inf)),
                                         reduce_max(select(valid,rdir.z,neg_inf)) );

          const float reduced_min_dist = reduce_min(select(valid,ray_tnear,vfloat<K>(pos_inf)));
          const float reduced_max_dist = reduce_max(select(valid,ray_tfar ,vfloat<K>(neg_inf)));

          init(reduced_min_org, reduced_max_org, reduced_min_rdir, reduced_max_rdir, reduced_min_dist, reduced_max_dist);
        }

        __forceinline void init(const Vec3fa& reduced_min_org,
                                const Vec3fa& reduced_max_org,
                                const Vec3fa& reduced_min_rdir,
                                const Vec3fa& reduced_max_rdir,
                                const float reduced_min_dist,
                                const float reduced_max_dist)
        {
          min_rdir = select(ge_mask(reduced_min_rdir, Vec3fa(zero)), reduced_min_rdir, reduced_max_rdir);
          max_rdir = select(ge_mask(reduced_min_rdir, Vec3fa(zero)), reduced_max_rdir, reduced_min_rdir);

          if (!robust)
          {
            min_org_rdir = min_rdir * select(ge_mask(reduced_min_rdir, Vec3fa(zero)), reduced_max_org, reduced_min_org);
            max_org_rdir = max_rdir * select(ge_mask(reduced_min_rdir, Vec3fa(zero)), reduced_min_org, reduced_max_org);
          }
          else
          {
            min_org_rdir = select(ge_mask(reduced_min_rdir, Vec3fa(zero)), reduced_max_org, reduced_min_org);
            max_org_rdir = select(ge_mask(reduced_min_rdir, Vec3fa(zero)), reduced_min_org, reduced_max_org);
          }

          min_dist = reduced_min_dist;
          max_dist = reduced_max_dist;
          
          nf = NearFarPreCompute(min_rdir,N);
        }

        __forceinline unsigned int intersectFast(const typename BVHN<N>::AlignedNode* __restrict__ node, vfloat<Nx>& dist) const
        {
          const vfloat<Nx> bminX = *(const vfloat<N>*)((const char*)&node->lower_x + nf.nearX);
          const vfloat<Nx> bminY = *(const vfloat<N>*)((const char*)&node->lower_x + nf.nearY);
          const vfloat<Nx> bminZ = *(const vfloat<N>*)((const char*)&node->lower_x + nf.nearZ);
          const vfloat<Nx> bmaxX = *(const vfloat<N>*)((const char*)&node->lower_x + nf.farX);
          const vfloat<Nx> bmaxY = *(const vfloat<N>*)((const char*)&node->lower_x + nf.farY);
          const vfloat<Nx> bmaxZ = *(const vfloat<N>*)((const char*)&node->lower_x + nf.farZ);
                    
          const vfloat<Nx> fminX = msub(bminX, vfloat<Nx>(min_rdir.x), vfloat<Nx>(min_org_rdir.x));
          const vfloat<Nx> fminY = msub(bminY, vfloat<Nx>(min_rdir.y), vfloat<Nx>(min_org_rdir.y));
          const vfloat<Nx> fminZ = msub(bminZ, vfloat<Nx>(min_rdir.z), vfloat<Nx>(min_org_rdir.z));
          const vfloat<Nx> fmaxX = msub(bmaxX, vfloat<Nx>(max_rdir.x), vfloat<Nx>(max_org_rdir.x));
          const vfloat<Nx> fmaxY = msub(bmaxY, vfloat<Nx>(max_rdir.y), vfloat<Nx>(max_org_rdir.y));
          const vfloat<Nx> fmaxZ = msub(bmaxZ, vfloat<Nx>(max_rdir.z), vfloat<Nx>(max_org_rdir.z));
          
          const vfloat<Nx> fmin  = maxi(fminX, fminY, fminZ, vfloat<Nx>(min_dist)); 
          dist = fmin;
          const vfloat<Nx> fmax  = mini(fmaxX, fmaxY, fmaxZ, vfloat<Nx>(max_dist));
          const vbool<Nx> vmask_node_hit = fmin <= fmax;
          size_t m_node = movemask(vmask_node_hit) & (((size_t)1 << N)-1);
          return m_node;          
        }

        __forceinline unsigned int intersectRobust(const typename BVHN<N>::AlignedNode* __restrict__ node, vfloat<Nx>& dist) const
        {
          const vfloat<Nx> bminX = *(const vfloat<N>*)((const char*)&node->lower_x + nf.nearX);
          const vfloat<Nx> bminY = *(const vfloat<N>*)((const char*)&node->lower_x + nf.nearY);
          const vfloat<Nx> bminZ = *(const vfloat<N>*)((const char*)&node->lower_x + nf.nearZ);
          const vfloat<Nx> bmaxX = *(const vfloat<N>*)((const char*)&node->lower_x + nf.farX);
          const vfloat<Nx> bmaxY = *(const vfloat<N>*)((const char*)&node->lower_x + nf.farY);
          const vfloat<Nx> bmaxZ = *(const vfloat<N>*)((const char*)&node->lower_x + nf.farZ);
                    
          const vfloat<Nx> fminX = (bminX - vfloat<Nx>(min_org_rdir.x)) * vfloat<Nx>(min_rdir.x);
          const vfloat<Nx> fminY = (bminY - vfloat<Nx>(min_org_rdir.y)) * vfloat<Nx>(min_rdir.y);
          const vfloat<Nx> fminZ = (bminZ - vfloat<Nx>(min_org_rdir.z)) * vfloat<Nx>(min_rdir.z);
          const vfloat<Nx> fmaxX = (bmaxX - vfloat<Nx>(max_org_rdir.x)) * vfloat<Nx>(max_rdir.x);
          const vfloat<Nx> fmaxY = (bmaxY - vfloat<Nx>(max_org_rdir.y)) * vfloat<Nx>(max_rdir.y);
          const vfloat<Nx> fmaxZ = (bmaxZ - vfloat<Nx>(max_org_rdir.z)) * vfloat<Nx>(max_rdir.z);
          
          const float round_down = 1.0f-2.0f*float(ulp); // FIXME: use per instruction rounding for AVX512
          const float round_up   = 1.0f+2.0f*float(ulp);
          const vfloat<Nx> fmin  = max(fminX, fminY, fminZ, vfloat<Nx>(min_dist)); 
          dist = fmin;
          const vfloat<Nx> fmax  = min(fmaxX, fmaxY, fmaxZ, vfloat<Nx>(max_dist));
          const vbool<Nx> vmask_node_hit = (round_down*fmin <= round_up*fmax);
          size_t m_node = movemask(vmask_node_hit) & (((size_t)1 << N)-1);
          return m_node;          
        }

        __forceinline unsigned int intersect(const typename BVHN<N>::AlignedNode* __restrict__ node, vfloat<Nx>& dist) const
        {
          if (robust) return intersectRobust(node,dist);
          else        return intersectFast(node,dist);
        }

        __forceinline void updateMaxDist(const vfloat<K>& ray_tfar) {
          max_dist = reduce_max(ray_tfar);
        }

        NearFarPreCompute nf;

        Vec3fa min_rdir; 
        Vec3fa max_rdir;
        
        Vec3fa min_org_rdir; 
        Vec3fa max_org_rdir; 

        float min_dist;
        float max_dist;
      };
  }
}