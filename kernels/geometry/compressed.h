#pragma once
#include <iostream>

#include "compressed_node.h"
#include "compressed_leaf.h"
#include "compressed_help.h"

#include "../bvh/node_intersector1.h"




namespace embree
{
    namespace isa
    {
        namespace compressed
        {

            // uncompressed nodes, only quantization
            // 3 bits for X, 3 bits for y, 2 bits for z
            // bounds per child stored and handled separately
            typedef Node<flavor::non,quantization::uni,quantization::uni,3,3,2> quantizedUniform;
            typedef Node<flavor::non,quantization::man,quantization::man2,3,3,2> quantizedNonUniform;

            // compressed nodes + quantization
            // 3 bits for X, 3 bits for y, 2 bits for z
            // 4 child bounding-boxes combined (outer and inner bounds)
            typedef Node<flavor::com,quantization::uni,quantization::uni,3,3,2> compressedUniform;
            typedef Node<flavor::com,quantization::man,quantization::man2,3,3,2> compressedNonUniform;

            // halfslab-compressed nodes + quantization
            // 3 bits for X, 3 bits for y, 2 bits for z
            // 4 child bounding-boxes combined (inner stored, outer reused from parent)
            typedef Node<flavor::mid,quantization::zero,quantization::uni,3,3,2> halfSlabUniform;
            typedef Node<flavor::mid,quantization::zero,quantization::man2,3,3,2> halfSlabNonUniform;

            // uncompressed not-quantized node
            // x, y, and z: each 32 bit float
            typedef Node<flavor::ref,quantization::zero,quantization::zero,32,32,32> fullPrecision;

            typedef quantTris<4> leaf_t;

            template<typename node_t, bool use_grid = false, bool use_leaf = false, bool use_box = true>
                class CompressedBVH
                {
                public:

                    CompressedBVH(const SubdivPatch1Base* patches,
                                  const unsigned x0, const unsigned x1, const unsigned y0, const unsigned y1, const unsigned subdiv_level,
                                  Scene* scene, node_t* nodes, leaf_t* leaves, Vec3f* grid, BBox3fa* bounds_o = nullptr)
                        : geomID(patches->geomID()), primID(patches->primID()), nodes(nodes), leaves(leaves), grid(grid)
                        {

                            const unsigned width = grid_width = x1-x0+1;
                            const unsigned height = y1-y0+1;

                            // prepare temporary storage
                            unsigned temp_size = width * height;
                            dynamic_large_stack_array(float,local_grid_u,temp_size,32*32*sizeof(float));
                            dynamic_large_stack_array(float,local_grid_v,temp_size,32*32*sizeof(float));
                            dynamic_large_stack_array(float,local_grid_x,temp_size,32*32*sizeof(float));
                            dynamic_large_stack_array(float,local_grid_y,temp_size,32*32*sizeof(float));
                            dynamic_large_stack_array(float,local_grid_z,temp_size,32*32*sizeof(float));


                            // compute vertex grid (+displacement)
                            evalGrid(patches[0],x0,x1,y0,y1,patches->grid_u_res,patches->grid_v_res,
                                     local_grid_x,local_grid_y,local_grid_z,
                                     local_grid_u,local_grid_v,scene->get<SubdivMesh>(patches->geomID()));

                            // copy to v array
                            std::vector<Vec3f> v(width*height);
                            for (size_t i = 0; i < width*height; ++i)
                                v[i] = Vec3f(local_grid_x[i], local_grid_y[i], local_grid_z[i]);


                            // indices of corner entries
                            unsigned i00 = 0;                   // u0 v0
                            unsigned i10 = width-1;             // u1 v0
                            unsigned i01 = width*(height-1);    // u0 v1
                            unsigned i11 = width*height-1;      // u1 v1

                            // store corner uvs and uv-stepping
                            uv[0] = Vec2f(local_grid_u[i00], local_grid_v[i00]);
                            uv[1] = Vec2f(local_grid_u[i11] - uv[0].x, local_grid_v[i11] - uv[0].y);

                            const float num_edges = 1 << subdiv_level;
                            rcp_edges = 1.f / num_edges;

                            // local coordinate frame corner vertices
                            Vec3fa frame_vertex00;
                            Vec3fa frame_vertex10;
                            Vec3fa frame_vertex01;
                            Vec3fa frame_vertex11;

                            // use displaced vertices for frame alignment
                            if (!use_leaf) {
                                frame_vertex00 = v[i00];
                                frame_vertex10 = v[i10];
                                frame_vertex01 = v[i01];
                                frame_vertex11 = v[i11];
                            }

                            // use un-displaced vertices for frame alignment
                            else {

                                const bool applyDisplacement = false;
                                evalGrid(patches[0],x0,x1,y0,y1,patches->grid_u_res,patches->grid_v_res,
                                         local_grid_x,local_grid_y,local_grid_z,
                                         local_grid_u,local_grid_v,scene->get<SubdivMesh>(patches->geomID()), applyDisplacement);

                                frame_vertex00 = Vec3fa(local_grid_x[i00], local_grid_y[i00], local_grid_z[i00]);
                                frame_vertex10 = Vec3fa(local_grid_x[i10], local_grid_y[i10], local_grid_z[i10]);
                                frame_vertex01 = Vec3fa(local_grid_x[i01], local_grid_y[i01], local_grid_z[i01]);
                                frame_vertex11 = Vec3fa(local_grid_x[i11], local_grid_y[i11], local_grid_z[i11]);
                            }


                            // local xy-sheared space with orthogonal z
                            const Vec3fa vx = normalize(frame_vertex10 - frame_vertex00 + frame_vertex11 - frame_vertex01);
                            const Vec3fa vy = normalize(frame_vertex01 - frame_vertex00 + frame_vertex11 - frame_vertex10);
                            const Vec3fa vz = normalize(cross(vx, vy));

                            const LinearSpace3f world = LinearSpace3f(vx, vy, vz);
                            space = world.inverse();


                            // get local projection
                            const Vec3f v00 = xfmPoint(space, v[i00]);
                            const Vec3f v10 = xfmPoint(space, v[i10]);
                            const Vec3f v01 = xfmPoint(space, v[i01]);
                            const Vec3f v11 = xfmPoint(space, v[i11]);

                            std::vector<Eigen::Vector2f> source(4);
                            source[0] = Eigen::Vector2f(v00.x, v00.y);
                            source[1] = Eigen::Vector2f(v10.x, v10.y);
                            source[2] = Eigen::Vector2f(v01.x, v01.y);
                            source[3] = Eigen::Vector2f(v11.x, v11.y);

                            std::vector<Eigen::Vector2f> target(4);
                            target[0] = Eigen::Vector2f(-1.f, -1.f);
                            target[1] = Eigen::Vector2f( 1.f, -1.f);
                            target[2] = Eigen::Vector2f(-1.f,  1.f);
                            target[3] = Eigen::Vector2f( 1.f,  1.f);

                            // check proper alignment of grid vertices
                            bool patchOK = !use_grid && true;
                            BBox3f lBox = empty;
                            for (size_t y = 0; y < height-1; ++y){
                                for (size_t x = 0; x < width-1; ++x) {
                                    const Vec3f v00 = xfmPoint(space, v[y*width+x]);
                                    const Vec3f v10 = xfmPoint(space, v[y*width+x+1]);
                                    const Vec3f v01 = xfmPoint(space, v[(y+1)*width+x]);
                                    const Vec3f v11 = xfmPoint(space, v[(y+1)*width+x+1]);

                                    lBox.extend(v00);
                                    lBox.extend(v10);
                                    lBox.extend(v01);
                                    lBox.extend(v11);

                                    if (v00.x > v10.x || v01.x > v11.x ||
                                        v00.y > v01.y || v10.y > v11.y)
                                        patchOK = false;
                                }
                            }

                            BBox3f pBox = empty;

                            // check if projections is applicable on patch
                            if (patchOK == true) {
                                proj = ComputeLinearEstimate(source, target);

                                for (auto it = v.begin(); it != v.end(); ++it) {
                                    const Vec3f p = project(xfmPoint(space, *it),proj);

                                    // check if projection breaks something
                                    if (!std::isfinite(p.x) || !std::isfinite(p.y) ||
                                        p.x < -1.5f || p.x > 1.5f ||
                                        p.y < -1.5f || p.y > 1.5f) {
                                        patchOK = false;
                                        break;
                                    }

                                    pBox.extend(p);
                                }
                            }

                            if (patchOK == true) {

                                // find scaling based on projected bounding box
                                // and apply to actual projection
                                source[0] = Eigen::Vector2f(pBox.lower.x, pBox.lower.y);
                                source[1] = Eigen::Vector2f(pBox.upper.x, pBox.lower.y);
                                source[2] = Eigen::Vector2f(pBox.lower.x, pBox.upper.y);
                                source[3] = Eigen::Vector2f(pBox.upper.x, pBox.upper.y);

                                proj = ComputeLinearEstimate(source, target) * proj;
                            }
                            else {

                                // find scaling based on not-projected bounding box
                                // omit projection
                                source[0] = Eigen::Vector2f(lBox.lower.x, lBox.lower.y);
                                source[1] = Eigen::Vector2f(lBox.upper.x, lBox.lower.y);
                                source[2] = Eigen::Vector2f(lBox.lower.x, lBox.upper.y);
                                source[3] = Eigen::Vector2f(lBox.upper.x, lBox.upper.y);

                                proj = ComputeLinearEstimate(source, target);
                            }

                            BBox3f projBox = empty;
                            *bounds_o = empty;


                            // build temp box hierarchy 
                            const Eigen::Matrix3f iproj = proj.inverse();

                            std::vector<std::vector<BBox3f>> hierarchy;
                            hierarchy.emplace_back(std::vector<BBox3f>());

                            build_box_hierarchy((height-1)*(width-1), width, v, hierarchy);
                            std::reverse(hierarchy.begin(), hierarchy.end());

                            // build node hierarchy
                            unsigned curr = 0;
                            for (size_t i = 0; i < hierarchy.size() -1; ++i) {
                                for (size_t k = 0; k < hierarchy[i].size(); ++k) {

                                    //build compressed node
                                    nodes[curr].setAABB(hierarchy[i][k],
                                                        hierarchy[i+1][k*4+0],
                                                        hierarchy[i+1][k*4+1],
                                                        hierarchy[i+1][k*4+2],
                                                        hierarchy[i+1][k*4+3]);

                                    // update hierarchy data with again uncompressed node
                                    for (size_t m = 0; m < 4; ++m) {
                                        BBox3f compBox;
                                        nodes[curr].getAABB(compBox, hierarchy[i][k], m);
                                        hierarchy[i+1][k*4+m] = compBox;

                                        // update root box with reconstructed leaf data
                                        if (i == hierarchy.size() - 2) {


                                            BBox3f tmpBox = BBox3f(project(compBox.lower, iproj));
                                            tmpBox.extend(project(Vec3f(compBox.upper.x, compBox.lower.y, compBox.lower.z), iproj));
                                            tmpBox.extend(project(Vec3f(compBox.lower.x, compBox.upper.y, compBox.lower.z), iproj));
                                            tmpBox.extend(project(Vec3f(compBox.upper.x, compBox.upper.y, compBox.lower.z), iproj));
                                            tmpBox.extend(project(Vec3f(compBox.lower.x, compBox.lower.y, compBox.upper.z), iproj));
                                            tmpBox.extend(project(Vec3f(compBox.upper.x, compBox.lower.y, compBox.upper.z), iproj));
                                            tmpBox.extend(project(Vec3f(compBox.lower.x, compBox.upper.y, compBox.upper.z), iproj));
                                            tmpBox.extend(project(compBox.upper, iproj));

                                            // extend unprojected box -> unprojected later
                                            projBox.extend(compBox);

                                            // extend return bounding box by every possible uncompressed vertex
                                            bounds_o->extend(xfmPoint(world, tmpBox.lower));
                                            bounds_o->extend(xfmPoint(world, Vec3f(tmpBox.upper.x, tmpBox.lower.y, tmpBox.lower.z)));
                                            bounds_o->extend(xfmPoint(world, Vec3f(tmpBox.lower.x, tmpBox.upper.y, tmpBox.lower.z)));
                                            bounds_o->extend(xfmPoint(world, Vec3f(tmpBox.upper.x, tmpBox.upper.y, tmpBox.lower.z)));
                                            bounds_o->extend(xfmPoint(world, Vec3f(tmpBox.upper.x, tmpBox.lower.y, tmpBox.upper.z)));
                                            bounds_o->extend(xfmPoint(world, Vec3f(tmpBox.lower.x, tmpBox.upper.y, tmpBox.upper.z)));
                                            bounds_o->extend(xfmPoint(world, Vec3f(tmpBox.upper.x, tmpBox.upper.y, tmpBox.upper.z)));
                                            bounds_o->extend(xfmPoint(world, tmpBox.upper));
                                        }
                                    }

                                    ++curr;
                                }
                            }



                            // store frustum box
                            const Vec3f p00 = project(Vec3f(projBox.lower.x, projBox.lower.y, projBox.lower.z), iproj);
                            const Vec3f p10 = project(Vec3f(projBox.upper.x, projBox.lower.y, projBox.lower.z), iproj);
                            const Vec3f p01 = project(Vec3f(projBox.lower.x, projBox.upper.y, projBox.upper.z), iproj);
                            const Vec3f p11 = project(Vec3f(projBox.upper.x, projBox.upper.y, projBox.upper.z), iproj);

                            box[0] = projBox.lower.z;
                            box[1] = projBox.upper.z;
                            box[2] = p00.x;
                            box[3] = p00.y;
                            box[4] = p10.x;
                            box[5] = p10.y;
                            box[6] = p01.x;
                            box[7] = p01.y;
                            box[8] = p11.x;
                            box[9] = p11.y;


                            // store pizza box approximation
                            if (use_leaf) {

                                // find max extent
                                extent = 0.f;
                                for (size_t i = 0; i < hierarchy.back().size(); ++i) {
                                    const unsigned x = DecodeMorton2X(i);
                                    const unsigned y = DecodeMorton2Y(i);

                                    extent = fmaxf(extent, leaves[i].estimateExtent(hierarchy.back()[i],
                                                                                    project(xfmPoint(space, v[y*width+x]), proj),
                                                                                    project(xfmPoint(space, v[y*width+x+1]), proj),
                                                                                    project(xfmPoint(space, v[(y+1)*width+x]), proj),
                                                                                    project(xfmPoint(space, v[(y+1)*width+x+1]), proj)));

                                }

#define MAX_EXTENT 1.0f
                                extent = fminf(extent, MAX_EXTENT);

                                // store z-values regarding found extent
                                for (size_t i = 0; i < hierarchy.back().size(); ++i) {
                                    const unsigned x = DecodeMorton2X(i);
                                    const unsigned y = DecodeMorton2Y(i);

                                    leaves[i].setZ(hierarchy.back()[i],
                                                   project(xfmPoint(space, v[y*width+x]), proj),
                                                   project(xfmPoint(space, v[y*width+x+1]), proj),
                                                   project(xfmPoint(space, v[(y+1)*width+x]), proj),
                                                   project(xfmPoint(space, v[(y+1)*width+x+1]), proj),
                                                   extent);
                                }
                            }

                            // store uncompressed vertex grid
                            if (use_grid) {
                                for (size_t i = 0; i < width*height; ++i) {
                                    grid[i] = v[i];

                                }
                            }

                        }


                    template<typename Allocator>
                        static CompressedBVH* create(const SubdivPatch1Base* patches,
                                                     const unsigned x0, const unsigned x1, const unsigned y0, const unsigned y1, const unsigned subdiv_level,
                                                     Scene* scene, Allocator& alloc, BBox3fa* bounds_o = nullptr)
                        {
                            // number of uv /intervals
                            const unsigned num_u_cells = x1-x0;
                            const unsigned num_v_cells = y1-y0;

                            // number of leaf-lvel voxels, vertices, hierarchy nodes
                            const unsigned num_cells = num_u_cells * num_v_cells;
                            const unsigned num_verts = num_cells + num_u_cells + num_v_cells + 1;
                            const unsigned num_nodes = (static_cast<size_t>(pow(4.f, subdiv_level) -1) / 3);

                            // basic struct size + hierarchy size
                            const size_t core_size = sizeof(CompressedBVH);
                            const size_t node_size = num_nodes * (sizeof(node_t));

                            // size of leaf-level approximation / vertex size
                            const size_t cell_size = use_leaf ? num_cells * (sizeof(leaf_t)) : 0;
                            const size_t vert_size = use_grid ? num_verts * (sizeof(Vec3f)) : 0;

                            // offset for data access
                            const size_t node_offset = core_size;
                            const size_t cell_offset = node_offset + node_size;
                            const size_t vert_offset = cell_offset + cell_size;

                            void* cbvh = alloc(core_size + node_size + cell_size + vert_size);

                            node_t* nodes =  reinterpret_cast<node_t*>(static_cast<char*>(cbvh) + node_offset);
                            leaf_t* leaves = reinterpret_cast<leaf_t*>(static_cast<char*>(cbvh) + cell_offset);
                            Vec3f* grid =    reinterpret_cast<Vec3f*> (static_cast<char*>(cbvh) + vert_offset);

                            CompressedBVH* ret = new (cbvh) CompressedBVH(patches, x0, x1, y0, y1, subdiv_level, scene, nodes, leaves, grid, bounds_o);
                            ret->elems = num_nodes;
                            return ret;
                        }


                private:

                    void build_box_hierarchy(const unsigned size, const unsigned width,
                                             const std::vector<Vec3f>& vertices, std::vector<std::vector<BBox3f>>& boxes) {

                        // fill leaf-level boxes
                        for (size_t i = 0; i < size; ++i) {

                            const unsigned x = DecodeMorton2X(i);
                            const unsigned y = DecodeMorton2Y(i);

                            boxes[0].emplace_back(BBox3f(project(xfmPoint(space, vertices[y*width+x]), proj)));
                            boxes[0].back().extend(project(xfmPoint(space, vertices[y*width+x+1]), proj));
                            boxes[0].back().extend(project(xfmPoint(space, vertices[(y+1)*width+x]), proj));
                            boxes[0].back().extend(project(xfmPoint(space, vertices[(y+1)*width+x+1]), proj));
                        }


                        // reconstruct up to root
                        while (boxes.back().size() > 1) {
                            std::vector<BBox3f>& curr = boxes.back();
                            std::vector<BBox3f> currLevel;
                            for (auto it = curr.begin(); it != curr.end(); it+=4)
                                currLevel.push_back(merge(*it, *(it+1), *(it+2), *(it+3)));
                            boxes.push_back(currLevel);
                        }
                    }


                public:

                    const unsigned geomID;
                    const unsigned primID;


                    Vec2f uv[2];
                    float rcp_edges;

                    node_t* nodes;
                    leaf_t* leaves;
                    Vec3f* grid;
                    float extent;

                    unsigned elems;
                    unsigned grid_width;
                    // local orientation
                    LinearSpace3f space;

                    // projection of oriented box corners`s x/y components to -1/1 range
                    // -> rectify patch's xy footprint
                    Eigen::Matrix3f proj;

                    // sheered box
                    float box[10];
                    //char alignment[7];

                };




            template<typename node_t, bool use_grid = false, bool use_leaf = false, bool use_box = true>
                class CompressedBVHIntersector1
                {
                public:
                    typedef CompressedBVH<node_t> Primitive;

                    class Precalculations
                    {
                    public:
                        __forceinline Precalculations (const Ray& ray, const void* ptr) {};
                    };


                    /*! Intersect a ray with the primitive. */
                    static __forceinline void intersect(Precalculations& pre, RayHit& ray, IntersectContext* context, const Primitive* prim, size_t ty, size_t& lazy_node) {
                        STAT3(normal.trav_prims, 1, 1, 1);
                        
                        // rotate ray into local frame
                        const Vec3f lOrg = xfmPoint(static_cast<LinearSpace3fa>(prim->space), ray.org);
                        const Vec3f lDir = xfmVector(static_cast<LinearSpace3fa>(prim->space), ray.dir);

                        // intersect frustum box, get entry and exit distance
                        float near = ray.tnear();
                        float far = ray.tfar;
                        if (not intersect_frustum(prim->box, lOrg, lDir, near, far))
                            return;

                        //// build projected ray inside frustum box
                        
                        // origin: projected entry point
                        Vec3fa org = project(lOrg+lDir*near, prim->proj);

                        // target: projected exit point
                        const Vec3fa tar = project(lOrg+lDir*far, prim->proj);

                        // local direction not yet normalized
                        Vec3fa dir = tar - org;
                        float zFactor;
                        bool special = false;
                        float tfar = 0.f;
                        const Vec3fa test = abs(dir);


                        /// handle special cases
                        const static float g_epsilon = 1.0E-4f;
                        // local frame very small
                        // use frustum box hit distance if hit found
                        if (unlikely(test.x < g_epsilon && test.y < g_epsilon && test.z < g_epsilon)) {
                            dir.z = copysign(1.f, lDir.z);
                            org.z -= dir.z;
                            zFactor = FLT_MAX;
                            tfar = FLT_MAX;
                        }
                        // local frame very flat
                        // special case handled separately
                        else if (unlikely(test.z < g_epsilon)) {
                            special = true;
                            tfar = length(dir);
                            dir = normalize(dir);
                        }
                        // default, go ahead
                        else {

                            // normalized dir: not usable for distance calculations anymore
                            // instead: use zFactor with z-component, which is not projected, nor scaled, nor translated.
                            dir = normalize(dir);
                            zFactor = lDir.z / dir.z;
                            tfar = (ray.tfar - near) * zFactor;
                        }

                        // stack size sufficient for compressed subd lvl 5
                        unsigned stack[16];
                        BBox3fa boxStack[16];
                        int sp = 0;

                        stack[0] = 0;

                        // local frame box xy: [-1, 1] and z from leaf data
                        boxStack[0] = BBox3fa(Vec3fa(-1.f, -1.f, prim->box[0]),
                                              Vec3fa( 1.f,  1.f, prim->box[1]));

                        // simd ray
                        typedef TravRay<4,4,true> tmpRay;
                        const tmpRay travRay(org, dir, 0.f, tfar);
                        vfloat4 tNear;
                        vfloat4 tFar;


                        // traversal loop
                        while (sp >= 0) {

                            // get current node and parent box
                            const unsigned& curr = stack[sp];
                            const BBox3fa& pBox = boxStack[sp--];

                            // check for inner node
                            if (curr >= prim->elems) {

                                // patch leaf approx
                                if (use_leaf)
                                {

                                    // get leaf-level index
                                    const unsigned idx = curr - prim->elems;
                                    const unsigned simdIdx = idx & 3;
                                    const float box_near = tNear[simdIdx];

                                    if (box_near >= tfar) continue;

                                    const float box_far = tFar[simdIdx];
                                    const float dimZ = pBox.upper.z - pBox.lower.z;
                                    float z1, z2, z3, z4;


                                    // reconstruct patch approximation
#define REFIT_LEAF_Z
#ifdef REFIT_LEAF_Z
                                    const float range = (1.f + 2.f * prim->extent) * dimZ;
                                    const float dz = prim->leaves[idx].z.getDelta() * range;
                                    prim->leaves[idx].getZ(z1, z2, z3, z4, range, pBox.lower.z - dimZ*prim->extent);
#else
                                    const float dz = prim->leaves[idx].z.getDelta() * dimZ;
                                    prim->leaves[idx].getZ(z1, z2, z3, z4, dimZ,pBox.lower.z);
#endif


                                    // intersect patch approx
                                    float u, v;
                                    float t = tfar;
                                    if (intersect_patch(idx, prim->rcp_edges, dz, box_near, box_far, z1, z2, z3, z4, pBox, org, dir, u, v, t)) {
                                        ray.u = prim->uv[0].x + u * prim->uv[1].x;
                                        ray.v = prim->uv[0].y + v * prim->uv[1].y;

                                        // dummy value -> we use smooth normals
                                        ray.Ng = Vec3f(1.f, 0.f, 0.f);

                                        ray.geomID = prim->geomID;
                                        ray.primID = prim->primID;
                                        ray.instID = context->instID;

                                        // update distance in local frame
                                        tfar = t;

                                        if (unlikely(special == true)) {
                                            // for flat frames, project hit point into global space got get distance
                                            const Vec3fa p = project(org + dir * t, prim->proj.inverse());
                                            ray.tfar = length(p - lOrg);
                                        }
                                        else
                                            // store fixed distance
                                            ray.tfar = t / zFactor + near;
                                    }

                                }


                                // full-precision vertex grid
                                if (use_grid)
                                {
                                    const unsigned idx = curr - prim->elems;
                                    if (intersect_triangles(idx, prim->rcp_edges, prim->grid_width, prim->grid, ray)) {
                                        ray.Ng = Vec3f(1.f, 0.f, 0.f);
                                        ray.u = prim->uv[0].x + ray.u * prim->uv[1].x;
                                        ray.v = prim->uv[0].y + ray.v * prim->uv[1].y;
                                        ray.geomID = prim->geomID;
                                        ray.primID = prim->primID;
                                        ray.instID = context->instID;

                                        tfar = (ray.tfar - near) * zFactor;
                                    }

                                }

                                // no leaf, box approx
                                if (use_box)
                                {
                                    const unsigned idx = curr - prim->elems;
                                    const float is = tNear[idx & 3];


                                    // use bbox hit as surface approximation
                                    if (is <= tfar) {

                                        // get uv coords vom xy hit points

                                        const float x = (((org.x + dir.x * is) - pBox.lower.x) / (pBox.upper.x - pBox.lower.x) +
                                                         (float)DecodeMorton2X(idx)) * prim->rcp_edges;

                                        const float y = (((org.y + dir.y * is) - pBox.lower.y) / (pBox.upper.y - pBox.lower.y) +
                                                         (float)DecodeMorton2Y(idx)) * prim->rcp_edges;

                                        ray.u = prim->uv[0].x + x * prim->uv[1].x;
                                        ray.v = prim->uv[0].y + y * prim->uv[1].y;


                                        // dummy value -> we use use smooth normals
                                        ray.Ng = Vec3f(1.0, 0.0, 0.0);

                                        ray.geomID = prim->geomID;
                                        ray.primID = prim->primID;
                                        ray.instID = context->instID;

                                        // update distance in local frame
                                        tfar = is;

                                        if (unlikely(special == true)) {
                                            // for flat frames, project hit point into global space got get distance
                                            const Vec3fa p = project(org + dir * is, prim->proj.inverse());
                                            ray.tfar = length(p - lOrg);
                                        }
                                        else
                                            // store fixed distance
                                            ray.tfar = is / zFactor + near;
                                    }
                                }

                                continue;
                            }


                            // index of first child
                            const unsigned next = curr * 4 + 1;

                            // get, uncompress, and intersect node 
                            const node_t& cNode = prim->nodes[curr];
                            const BVH4::AlignedNode eNode = cNode.getNode(pBox);
                            size_t mask = intersectNodeRobust(&eNode, travRay, tNear, tFar);

                            // nothing hit
                            if (unlikely(mask == 0))
                                continue;

                            // one child hit
                            size_t r = __bscf(mask);
                            if (likely(mask==0)) {
                                boxStack[++sp] = eNode.bounds(r);

                                stack[sp] = next + r;

                                continue;
                            }

                            // two children hit
                            const size_t c0 = r;
                            const unsigned int d0 = ((unsigned int*)&tNear)[r];
                            r = __bscf(mask);
                            const unsigned int d1 = ((unsigned int*)&tNear)[r];
                            if (likely(mask==0)) {

                                //put uncompressed reference boxes and indices on stack
                                if (d0 > d1) {
                                    boxStack[++sp] = eNode.bounds(c0); stack[sp] = next + c0;
                                    boxStack[++sp] = eNode.bounds( r); stack[sp] = next +  r;
                                }
                                else {
                                    boxStack[++sp] = eNode.bounds( r); stack[sp] = next +  r;
                                    boxStack[++sp] = eNode.bounds(c0); stack[sp] = next + c0;
                                }
                                continue;
                            }

                            // three children hit
                            const size_t c1 = r;
                            r = __bscf(mask);
                            const unsigned int d2 = ((unsigned int*)&tNear)[r];
                            if (likely(mask==0)) {

                                // calculate stack position
                                const unsigned shift0 = sp + 1 + (d0 <= d1) + (d0 <= d2);
                                const unsigned shift1 = sp + 1 + (d1 <  d0) + (d1 <= d2);
                                const unsigned shift2 = sp + 1 + (d2 <  d0) + (d2 <  d1);

                                // put uncompressed reference boxes on stack
                                boxStack[shift0] = eNode.bounds(c0);
                                boxStack[shift1] = eNode.bounds(c1);
                                boxStack[shift2] = eNode.bounds(r);

                                // put node indices on stack
                                stack[shift0] = next + c0;
                                stack[shift1] = next + c1;
                                stack[shift2] = next + r;

                                sp += 3;
                                continue;
                            }

                            // four children hit
                            const size_t c2 = r;
                            r = __bscf(mask);
                            const unsigned int d3 = ((unsigned int*)&tNear)[r];

                            // calculate stack position
                            const unsigned shift0 = sp + 1 + (d0 <= d1) + (d0 <= d2) + (d0 <= d3);
                            const unsigned shift1 = sp + 1 + (d1 <  d0) + (d1 <= d2) + (d1 <= d3);
                            const unsigned shift2 = sp + 1 + (d2 <  d0) + (d2 <  d1) + (d2 <= d3);
                            const unsigned shift3 = sp + 1 + (d3 <  d0) + (d3 <  d1) + (d3 <  d2);

                            // put uncompressed reference boxes on stack
                            boxStack[shift0] = eNode.bounds(c0);
                            boxStack[shift1] = eNode.bounds(c1);
                            boxStack[shift2] = eNode.bounds(c2);
                            boxStack[shift3] = eNode.bounds(r);

                            // put node indices on stack
                            stack[shift0] = next + c0;
                            stack[shift1] = next + c1;
                            stack[shift2] = next + c2;
                            stack[shift3] = next + r;

                            sp += 4;
                        }
                                
                    }

                    static __forceinline bool occluded( Precalculations& pre, Ray& ray, IntersectContext* context, const Primitive* prim, size_t ty, size_t& lazy_node) {
                        return true;
                    }

                };



// All quantization and compression types of our previous approach, not used currently.
// Our new approach is based on the first one, only.
            typedef CompressedBVHIntersector1<compressedNonUniform> SubdivPatch1Oriented_CompressedNonUniformIntersector1;
//#define COMPRESSED_USE_ALL
#ifdef COMPRESSED_USE_ALL
            typedef CompressedBVHIntersector1<quantizedUniform>     SubdivPatch1Oriented_QuantizedUniformIntersector1;
            typedef CompressedBVHIntersector1<quantizedNonUniform>  SubdivPatch1Oriented_QuantizedNonUniformIntersector1;
            typedef CompressedBVHIntersector1<compressedUniform>    SubdivPatch1Oriented_CompressedUniformIntersector1;
            typedef CompressedBVHIntersector1<halfSlabUniform>      SubdivPatch1Oriented_HalfSlabUniformIntersector1;
            typedef CompressedBVHIntersector1<halfSlabNonUniform>   SubdivPatch1Oriented_HalfSlabNonUniformIntersector1;
#endif

            typedef CompressedBVHIntersector1<fullPrecision>        SubdivPatch1Oriented_FullPrecisionIntersector1;

// Our new node types
            typedef CompressedBVHIntersector1<compressedNonUniform, true, false, false> SubdivPatch1cGridIntersector1; // compressed grid
            typedef CompressedBVHIntersector1<compressedNonUniform, false, true, false> SubdivPatch1cLeafIntersector1; // compressed with leaf approximation
            typedef CompressedBVHIntersector1<compressedNonUniform, false, false, true> SubdivPatch1cBoxIntersector1;  // compressed with voxels
        }
    }
}



// vim: set foldmethod=marker :
