// ======================================================================== //
// Copyright 2009-2018 Intel Corporation                                    //
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

#include "../common/math/random_sampler.h"
#include "../common/core/differential_geometry.h"
#include "../common/tutorial/tutorial_device.h"
#include "../common/tutorial/scene_device.h"

namespace embree {

extern "C" ISPCScene* g_ispc_scene;
extern "C" bool g_changed;
extern "C" int g_instancing_mode;

extern "C" bool g_adjustedCoherentBench;
extern "C" bool g_adjustedIncoherentBench;
Vec3fa* g_randomRayOrgs = nullptr;
Vec3fa* g_randomRayDirs = nullptr;
unsigned g_width = 0;


/* scene data */
RTCDevice g_device = nullptr;
RTCScene g_scene = nullptr;
bool g_subdiv_mode = false;

#define SPP 1

//#define FIXED_EDGE_TESSELLATION_VALUE 3
#define FORCE_FIXED_SUBD_LEVEL

#define MAX_EDGE_LEVEL 64.0f
#define MIN_EDGE_LEVEL  4.0f
#define LEVEL_FACTOR   64.0f

inline float updateEdgeLevel( ISPCSubdivMesh* mesh, const Vec3fa& cam_pos, const unsigned int e0, const unsigned int e1)
{
  const Vec3fa v0 = mesh->positions[0][mesh->position_indices[e0]];
  const Vec3fa v1 = mesh->positions[0][mesh->position_indices[e1]];
  const Vec3fa edge = v1-v0;
  const Vec3fa P = 0.5f*(v1+v0);
  const Vec3fa dist = cam_pos - P;
  return max(min(LEVEL_FACTOR*(0.5f*length(edge)/length(dist)),MAX_EDGE_LEVEL),MIN_EDGE_LEVEL);
}


void updateEdgeLevelBuffer( ISPCSubdivMesh* mesh, const Vec3fa& cam_pos, unsigned int startID, unsigned int endID )
{
  for (unsigned int f=startID; f<endID;f++) {
    unsigned int e = mesh->face_offsets[f];
    unsigned int N = mesh->verticesPerFace[f];
    if (N == 4) /* fast path for quads */
      for (unsigned int i=0; i<4; i++)
        mesh->subdivlevel[e+i] =  updateEdgeLevel(mesh,cam_pos,e+(i+0),e+(i+1)%4);
    else if (N == 3) /* fast path for triangles */
      for (unsigned int i=0; i<3; i++)
        mesh->subdivlevel[e+i] =  updateEdgeLevel(mesh,cam_pos,e+(i+0),e+(i+1)%3);
    else /* fast path for general polygons */
      for (unsigned int i=0; i<N; i++)
        mesh->subdivlevel[e+i] =  updateEdgeLevel(mesh,cam_pos,e+(i+0),e+(i+1)%N);
  }
}

#if defined(ISPC)
void updateSubMeshEdgeLevelBufferTask (int taskIndex, int threadIndex,  ISPCSubdivMesh* mesh, const Vec3fa& cam_pos )
{
  const unsigned int size = mesh->numFaces;
  const unsigned int startID = ((taskIndex+0)*size)/taskCount;
  const unsigned int endID   = ((taskIndex+1)*size)/taskCount;
  updateEdgeLevelBuffer(mesh,cam_pos,startID,endID);
}
void updateMeshEdgeLevelBufferTask (int taskIndex, int threadIndex,  ISPCScene* scene_in, const Vec3fa& cam_pos )
{
  ISPCGeometry* geometry = g_ispc_scene->geometries[taskIndex];
  if (geometry->type != SUBDIV_MESH) return;
  ISPCSubdivMesh* mesh = (ISPCSubdivMesh*) geometry;
  unsigned int geomID = mesh->geom.geomID;
  if (mesh->numFaces < 10000) {
    updateEdgeLevelBuffer(mesh,cam_pos,0,mesh->numFaces);
    rtcUpdateGeometryBuffer(geometry->geometry, RTC_BUFFER_TYPE_LEVEL, 0);
  }
  rtcCommitGeometry(geometry->geometry);
}
#endif

void updateEdgeLevels(ISPCScene* scene_in, const Vec3fa& cam_pos)
{
  /* first update small meshes */
#if defined(ISPC)
  parallel_for(size_t(0),size_t( scene_in->numGeometries ),[&](const range<size_t>& range) {
    const int threadIndex = (int)TaskScheduler::threadIndex();
    for (size_t i=range.begin(); i<range.end(); i++)
      updateMeshEdgeLevelBufferTask((int)i,threadIndex,scene_in,cam_pos);
  }); 
#endif

  /* now update large meshes */
  for (unsigned int g=0; g<scene_in->numGeometries; g++)
  {
    ISPCGeometry* geometry = g_ispc_scene->geometries[g];
    if (geometry->type != SUBDIV_MESH) continue;
    ISPCSubdivMesh* mesh = (ISPCSubdivMesh*) geometry;
#if defined(ISPC)
    if (mesh->numFaces < 10000) continue;
    parallel_for(size_t(0),size_t( (mesh->numFaces+4095)/4096 ),[&](const range<size_t>& range) {
    const int threadIndex = (int)TaskScheduler::threadIndex();
    for (size_t i=range.begin(); i<range.end(); i++)
      updateSubMeshEdgeLevelBufferTask((int)i,threadIndex,mesh,cam_pos);
  }); 
#else
    updateEdgeLevelBuffer(mesh,cam_pos,0,mesh->numFaces);
#endif
    rtcUpdateGeometryBuffer(geometry->geometry, RTC_BUFFER_TYPE_LEVEL, 0);
    rtcCommitGeometry(geometry->geometry);
  }
}

// compute normals after traversal
bool g_use_smooth_normals = true;
void device_key_pressed_handler(int key)
{
  if (key == 110 /*n*/) g_use_smooth_normals = !g_use_smooth_normals;
  else device_key_pressed_default(key);
}

RTCScene convertScene(ISPCScene* scene_in)
{
  for (unsigned int i=0; i<scene_in->numGeometries; i++)
  {
    ISPCGeometry* geometry = scene_in->geometries[i];
    if (geometry->type == SUBDIV_MESH) {
      g_subdiv_mode = true; break;
    }
  }

  RTCScene scene_out = ConvertScene(g_device, g_ispc_scene, RTC_BUILD_QUALITY_MEDIUM);

  /* commit individual objects in case of instancing */
  if (g_instancing_mode == ISPC_INSTANCING_SCENE_GEOMETRY || g_instancing_mode == ISPC_INSTANCING_SCENE_GROUP)
  {
    for (unsigned int i=0; i<scene_in->numGeometries; i++) {
      if (scene_in->geomID_to_scene[i]) rtcCommitScene(scene_in->geomID_to_scene[i]);
    }
  }

  /* commit changes to scene */
  return scene_out;
}


void postIntersectGeometry(const Ray& ray, DifferentialGeometry& dg, ISPCGeometry* geometry, int& materialID)
{
  if (geometry->type == TRIANGLE_MESH)
  {
    ISPCTriangleMesh* mesh = (ISPCTriangleMesh*) geometry;
    materialID = mesh->geom.materialID;
  }
  else if (geometry->type == QUAD_MESH)
  {
    ISPCQuadMesh* mesh = (ISPCQuadMesh*) geometry;
    materialID = mesh->geom.materialID;
  }
  else if (geometry->type == SUBDIV_MESH)
  {
    ISPCSubdivMesh* mesh = (ISPCSubdivMesh*) geometry;
    materialID = mesh->geom.materialID;
  }
  else if (geometry->type == CURVES)
  {
    ISPCHairSet* mesh = (ISPCHairSet*) geometry;
    materialID = mesh->geom.materialID;
  }
  else if (geometry->type == GROUP) {
    unsigned int geomID = ray.geomID; {
      postIntersectGeometry(ray,dg,((ISPCGroup*) geometry)->geometries[geomID],materialID);
    }
  }
  else
    assert(false);
}

AffineSpace3fa calculate_interpolated_space (ISPCInstance* instance, float gtime)
{
  if (instance->numTimeSteps == 1)
    return AffineSpace3fa(instance->spaces[0]);

   /* calculate time segment itime and fractional time ftime */
  const int time_segments = instance->numTimeSteps-1;
  const float time = gtime*(float)(time_segments);
  const int itime = clamp((int)(floor(time)),(int)0,time_segments-1);
  const float ftime = time - (float)(itime);
  return (1.0f-ftime)*AffineSpace3fa(instance->spaces[itime+0]) + ftime*AffineSpace3fa(instance->spaces[itime+1]);
}

inline int postIntersect(const Ray& ray, DifferentialGeometry& dg)
{
  int materialID = 0;
  unsigned int instID = ray.instID; {
    unsigned int geomID = ray.geomID; {
      ISPCGeometry* geometry = nullptr;
      if (g_instancing_mode == ISPC_INSTANCING_SCENE_GEOMETRY || g_instancing_mode == ISPC_INSTANCING_SCENE_GROUP) {
        ISPCInstance* instance = g_ispc_scene->geomID_to_inst[instID];
        geometry = g_ispc_scene->geometries[instance->geom.geomID];
      } else {
        geometry = g_ispc_scene->geometries[geomID];
      }
      postIntersectGeometry(ray,dg,geometry,materialID);
    }
  }

  if (g_instancing_mode != ISPC_INSTANCING_NONE)
  {
    unsigned int instID = ray.instID;
    {
      /* get instance and geometry pointers */
      ISPCInstance* instance = g_ispc_scene->geomID_to_inst[instID];

      /* convert normals */
      //AffineSpace3fa space = (1.0f-ray.time())*AffineSpace3fa(instance->space0) + ray.time()*AffineSpace3fa(instance->space1);
      AffineSpace3fa space = calculate_interpolated_space(instance,ray.time());
      dg.Ng = xfmVector(space,dg.Ng);
      dg.Ns = xfmVector(space,dg.Ns);
    }
  }

  return materialID;
}

inline Vec3fa face_forward(const Vec3fa& dir, const Vec3fa& _Ng) {
  const Vec3fa Ng = _Ng;
  return dot(dir,Ng) < 0.0f ? Ng : neg(Ng);
}



/* task that renders a single screen tile */
Vec3fa renderPixelStandard(float x, float y, const ISPCCamera& camera, RayStats& stats)
{
  /* initialize sampler */
  RandomSampler sampler;
  RandomSampler_init(sampler, (int)x, (int)y, 0);

  /* initialize ray */
  Ray ray(Vec3fa(camera.xfm.p), Vec3fa(normalize(x*camera.xfm.l.vx + y*camera.xfm.l.vy + camera.xfm.l.vz)), 0.0f, inf, RandomSampler_get1D(sampler));

  /* intersect ray with scene */
  RTCIntersectContext context;
  rtcInitIntersectContext(&context);
  context.flags = g_iflags_coherent;
  rtcIntersect1(g_scene,&context,RTCRayHit_(ray));
  RayStats_addRay(stats);

  /* shade background black */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) {
    return Vec3fa(0.0f);
  }

  /* shade all rays that hit something */
  Vec3fa color = Vec3fa(0.5f);

  /* compute differential geometry */
  DifferentialGeometry dg;
  dg.geomID = ray.geomID;
  dg.primID = ray.primID;
  dg.u = ray.u;
  dg.v = ray.v;
  dg.P  = ray.org+ray.tfar*ray.dir;
  dg.Ng = ray.Ng;
  dg.Ns = ray.Ng;

  if (g_use_smooth_normals)
    if (ray.geomID != RTC_INVALID_GEOMETRY_ID) // FIXME: workaround for ISPC bug, location reached with empty execution mask
    {
      Vec3fa dPdu,dPdv;
      unsigned int geomID = ray.geomID; {
        rtcInterpolate1(rtcGetGeometry(g_scene,geomID),ray.primID,ray.u,ray.v,RTC_BUFFER_TYPE_VERTEX,0,nullptr,&dPdu.x,&dPdv.x,3);
      }
      dg.Ns = cross(dPdv,dPdu);
    }

  int materialID = postIntersect(ray,dg);
  dg.Ng = face_forward(ray.dir,normalize(dg.Ng));
  dg.Ns = face_forward(ray.dir,normalize(dg.Ns));

  /* shade */
  if (g_ispc_scene->materials[materialID]->type == MATERIAL_OBJ) {
    ISPCOBJMaterial* material = (ISPCOBJMaterial*) g_ispc_scene->materials[materialID];
    color = Vec3fa(material->Kd);
  }

  return color*dot(neg(ray.dir),dg.Ns);
}


// current renderPixel function
typedef Vec3fa (* renderPixelFunc)(float x, float y, const ISPCCamera& camera, RayStats& stats);
renderPixelFunc renderPixel;

Vec3fa renderPixelCoherent(float x, float y, const ISPCCamera& camera, RayStats& stats) {
  /* initialize sampler */
  RandomSampler sampler;
  RandomSampler_init(sampler, (int)x, (int)y, 0);

  /* initialize ray */
  Ray ray(Vec3fa(camera.xfm.p), Vec3fa(normalize(x*camera.xfm.l.vx + y*camera.xfm.l.vy + camera.xfm.l.vz)), 0.0f, inf, RandomSampler_get1D(sampler));

  /* intersect ray with scene */
  RTCIntersectContext context;
  rtcInitIntersectContext(&context);
  context.flags = g_iflags_coherent;
  rtcIntersect1(g_scene,&context,RTCRayHit_(ray));
  RayStats_addRay(stats);

  /* shade background black */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) {
    return Vec3fa(0.0f);
  }
  /* shade all rays that hit something */
  return Vec3fa(0.5f);

  /* nothign else */

}


Vec3fa renderPixelIncoherent(float x, float y, const ISPCCamera& camera, RayStats& stats)
{

  /* initialize sampler */
  RandomSampler sampler;
  RandomSampler_init(sampler, (int)x, (int)y, 0);

  /* initialize ray */
  const int idx = static_cast<int>(y) * g_width + static_cast<int>(x);
  Ray ray(g_randomRayOrgs[idx], g_randomRayDirs[idx], 0.0f, inf, RandomSampler_get1D(sampler));

  /* intersect ray with scene */
  RTCIntersectContext context;
  rtcInitIntersectContext(&context);
  context.flags = g_iflags_incoherent;
  rtcIntersect1(g_scene,&context,RTCRayHit_(ray));
  RayStats_addRay(stats);

  /* shade background black */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) {
      return Vec3fa(0.0f);
  }
  /* shade all rays that hit something */
  return Vec3fa(0.5f);

  /* nothign else */
}

void makeRandomRay(Vec3fa& org, Vec3fa& dir, const BBox3fa box) {

  Vec3fa pos1, pos2, diam;
  diam = box.upper - box.lower;


  float X = drand48();
  float Y = drand48();
  float Z = drand48();

  pos1 = Vec3fa( X*diam.x + box.lower.x,
                 Y*diam.y + box.lower.y,
                 Z*diam.z + box.lower.z);

  X = drand48();
  Y = drand48();
  Z = drand48();

  pos2 = Vec3fa( X*diam.x + box.lower.x,
                 Y*diam.y + box.lower.y,
                 Z*diam.z + box.lower.z);

  
  org = pos1;
  dir = normalize(pos2 - pos1);
};

void prepareRandomRays(const unsigned int width, const unsigned int height) {
  g_randomRayOrgs = static_cast<Vec3fa*>(malloc(sizeof(Vec3fa)*width*height));
  g_randomRayDirs = static_cast<Vec3fa*>(malloc(sizeof(Vec3fa)*width*height));
  g_width = width;

  BBox3fa global_bounding_box = empty;

  ISPCGeometry** geometries = g_ispc_scene->geometries;

  for (size_t i = 0 ; i < g_ispc_scene->numGeometries; ++i) {
      if (geometries[i]->type == SUBDIV_MESH) {

          ISPCSubdivMesh* mesh = (ISPCSubdivMesh*)geometries[i];
          Vec3fa* vertices = *mesh->positions;

          for (size_t k = 0; k < mesh->numVertices; ++k) {
              global_bounding_box.extend(vertices[k]);

          }
      }
      if (geometries[i]->type == TRIANGLE_MESH) {
          ISPCTriangleMesh* mesh = (ISPCTriangleMesh*)geometries[i];
          Vec3fa* vertices = *mesh->positions;

          for (size_t k = 0; k < mesh->numVertices; ++k) {
              global_bounding_box.extend(vertices[k]);
          }
      }

  }

  for (size_t y = 0; y < height; ++y)
      for (size_t x = 0; x < width; ++x)
          makeRandomRay(g_randomRayOrgs[y*width+x], g_randomRayDirs[y*width+x], global_bounding_box);

}





void renderTileStandard(int taskIndex,
                        int threadIndex,
                        int* pixels,
                        const unsigned int width,
                        const unsigned int height,
                        const float time,
                        const ISPCCamera& camera,
                        const int numTilesX,
                        const int numTilesY)
{
  const int t = taskIndex;
  const unsigned int tileY = t / numTilesX;
  const unsigned int tileX = t - tileY * numTilesX;
  const unsigned int x0 = tileX * TILE_SIZE_X;
  const unsigned int x1 = min(x0+TILE_SIZE_X,width);
  const unsigned int y0 = tileY * TILE_SIZE_Y;
  const unsigned int y1 = min(y0+TILE_SIZE_Y,height);

  for (unsigned int y=y0; y<y1; y++) for (unsigned int x=x0; x<x1; x++)
  {
    Vec3fa color = renderPixel((float)x,(float)y,camera,g_stats[threadIndex]);

    /* write color to framebuffer */
    unsigned int r = (unsigned int) (255.0f * clamp(color.x,0.0f,1.0f));
    unsigned int g = (unsigned int) (255.0f * clamp(color.y,0.0f,1.0f));
    unsigned int b = (unsigned int) (255.0f * clamp(color.z,0.0f,1.0f));
    pixels[y*width+x] = (b << 16) + (g << 8) + r;
  }
}

/* task that renders a single screen tile */
void renderTileTask (int taskIndex, int threadIndex, int* pixels,
                         const unsigned int width,
                         const unsigned int height,
                         const float time,
                         const ISPCCamera& camera,
                         const int numTilesX,
                         const int numTilesY)
{
  renderTile(taskIndex,threadIndex,pixels,width,height,time,camera,numTilesX,numTilesY);
}

Vec3fa old_p;



/* called by the C++ code for initialization */
extern "C" void device_init (char* cfg)
{
  /* create new Embree device */
  g_device = rtcNewDevice(cfg);
  error_handler(nullptr,rtcGetDeviceError(g_device));

  /* set error handler */
  rtcSetDeviceErrorFunction(g_device,error_handler,nullptr);

  /* set start render mode */
  renderTile = renderTileStandard;
  key_pressed_handler = device_key_pressed_handler;
  old_p = Vec3fa(1E10);

  // set renderPixel function
  if (g_adjustedCoherentBench) 
    renderPixel = renderPixelCoherent;
  else if (g_adjustedIncoherentBench)
    renderPixel = renderPixelIncoherent;
  else
    renderPixel = renderPixelStandard;

}

/* called by the C++ code to render */
extern "C" void device_render (int* pixels,
                           const unsigned int width,
                           const unsigned int height,
                           const float time,
                           const ISPCCamera& camera)
{
  bool camera_changed = g_changed; g_changed = false;




  /* create scene */
  if (g_scene == nullptr) {
    g_scene = convertScene(g_ispc_scene);
    if (g_subdiv_mode) updateEdgeLevels(g_ispc_scene, camera.xfm.p);
    rtcCommitScene (g_scene);
    old_p = camera.xfm.p;
  }

  else
  {
    /* check if camera changed */
    if (ne(camera.xfm.p,old_p)) {
      camera_changed = true;
      old_p = camera.xfm.p;
    }

#ifndef FORCE_FIXED_SUBD_LEVEL
    /* update edge levels if camera changed */
    if (camera_changed && g_subdiv_mode) {
      updateEdgeLevels(g_ispc_scene,camera.xfm.p);
      rtcCommitScene (g_scene);
    }
#endif
  }

  /* create random rays */
  if (g_adjustedIncoherentBench && (g_randomRayOrgs == nullptr || g_randomRayDirs == nullptr))
      prepareRandomRays(width, height);

  /* render image */
  const int numTilesX = (width +TILE_SIZE_X-1)/TILE_SIZE_X;
  const int numTilesY = (height+TILE_SIZE_Y-1)/TILE_SIZE_Y;
  parallel_for(size_t(0),size_t(numTilesX*numTilesY),[&](const range<size_t>& range) {
    const int threadIndex = (int)TaskScheduler::threadIndex();
    for (size_t i=range.begin(); i<range.end(); i++)
      renderTileTask((int)i,threadIndex,pixels,width,height,time,camera,numTilesX,numTilesY);
  }); 
  //rtcDebug();
}

/* called by the C++ code for cleanup */
extern "C" void device_cleanup ()
{
  rtcReleaseScene (g_scene); g_scene = nullptr;
  rtcReleaseDevice(g_device); g_device = nullptr;

  free(g_randomRayOrgs);
  free(g_randomRayDirs);
}

} // namespace embree
