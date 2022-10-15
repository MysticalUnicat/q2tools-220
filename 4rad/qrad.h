/*
===========================================================================
Copyright (C) 1997-2006 Id Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
===========================================================================
*/

#include "cmdlib.h"
#include "mathlib.h"
#include "bspfile.h"
#include "polylib.h"
#include "threads.h"
#include "lbmlib.h"

#ifdef WIN32
#include <windows.h>
#endif

struct SH1 {
  float f[12];
};

static inline void SH1_Clear(struct SH1 *sh1) { memset(sh1, 0, sizeof(*sh1)); }

static inline struct SH1 SH1_FromDirectionalLight(const float direction[3], const float color[3]) {
  return (struct SH1){.f[0] = color[0],
                      .f[1] = direction[0] * color[0],
                      .f[2] = direction[1] * color[0],
                      .f[3] = direction[2] * color[0],
                      .f[4] = color[1],
                      .f[5] = direction[0] * color[1],
                      .f[6] = direction[1] * color[1],
                      .f[7] = direction[2] * color[1],
                      .f[8] = color[2],
                      .f[9] = direction[0] * color[2],
                      .f[10] = direction[1] * color[2],
                      .f[11] = direction[2] * color[2]};
}

static inline struct SH1 SH1_Add(const struct SH1 sh1_left, const struct SH1 sh1_right) {
  return (struct SH1){.f[0] = sh1_left.f[0] + sh1_right.f[0],
                      .f[1] = sh1_left.f[1] + sh1_right.f[1],
                      .f[2] = sh1_left.f[2] + sh1_right.f[2],
                      .f[3] = sh1_left.f[3] + sh1_right.f[3],
                      .f[4] = sh1_left.f[4] + sh1_right.f[4],
                      .f[5] = sh1_left.f[5] + sh1_right.f[5],
                      .f[6] = sh1_left.f[6] + sh1_right.f[6],
                      .f[7] = sh1_left.f[7] + sh1_right.f[7],
                      .f[8] = sh1_left.f[8] + sh1_right.f[8],
                      .f[9] = sh1_left.f[9] + sh1_right.f[9],
                      .f[10] = sh1_left.f[10] + sh1_right.f[10],
                      .f[11] = sh1_left.f[11] + sh1_right.f[11]};
}

static inline struct SH1 SH1_Scale(const struct SH1 sh1, float scale) {
  return (struct SH1){.f[0] = sh1.f[0] * scale,
                      .f[1] = sh1.f[1] * scale,
                      .f[2] = sh1.f[2] * scale,
                      .f[3] = sh1.f[3] * scale,
                      .f[4] = sh1.f[4] * scale,
                      .f[5] = sh1.f[5] * scale,
                      .f[6] = sh1.f[6] * scale,
                      .f[7] = sh1.f[7] * scale,
                      .f[8] = sh1.f[8] * scale,
                      .f[9] = sh1.f[9] * scale,
                      .f[10] = sh1.f[10] * scale,
                      .f[11] = sh1.f[11] * scale};
}

static inline struct SH1 SH1_ColorScale(const struct SH1 sh1, const float scale[3]) {
  return (struct SH1){.f[0] = sh1.f[0] * scale[0],
                      .f[1] = sh1.f[1] * scale[0],
                      .f[2] = sh1.f[2] * scale[0],
                      .f[3] = sh1.f[3] * scale[0],
                      .f[4] = sh1.f[4] * scale[1],
                      .f[5] = sh1.f[5] * scale[1],
                      .f[6] = sh1.f[6] * scale[1],
                      .f[7] = sh1.f[7] * scale[1],
                      .f[8] = sh1.f[8] * scale[2],
                      .f[9] = sh1.f[9] * scale[2],
                      .f[10] = sh1.f[10] * scale[2],
                      .f[11] = sh1.f[11] * scale[2]};
}

static inline void SH1_Sample(const struct SH1 sh1, const float direction[3], float output_color[3]) {
  // https://grahamhazel.com/blog/
  for(int component = 0; component < 3; component++) {
    float r0 = sh1.f[component * 4 + 0];
    float r1[3];
    r1[0] = sh1.f[component * 4 + 1];
    r1[1] = sh1.f[component * 4 + 2];
    r1[2] = sh1.f[component * 4 + 3];
    float r1_length_sq = r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2];
    if(r1_length_sq <= __FLT_EPSILON__ || r0 <= __FLT_EPSILON__) {
      output_color[component] = 0;
      continue;
    }
    float r1_length = sqrt(r1_length_sq);
    float one_over_r1_length = 1.0f / r1_length;
    float r1_normalized[3];
    r1_normalized[0] = r1[0] * one_over_r1_length;
    r1_normalized[1] = r1[1] * one_over_r1_length;
    r1_normalized[2] = r1[2] * one_over_r1_length;
    float q = 0.5f *
              (1 + r1_normalized[0] * direction[0] + r1_normalized[1] * direction[1] + r1_normalized[2] * direction[2]);
    float r1_length_over_r0 = r1_length / r0;
    float p = 1 + 2 * r1_length_over_r0;
    float a = (1 - r1_length_over_r0) / (1 + r1_length_over_r0);
    output_color[component] = r0 * (1 + (1 - a) * (p + 1) * powf(q, p)) * 0.25;
  }
}

typedef enum {
    emit_surface,
    emit_point,
    emit_spotlight,
    emit_sky
} emittype_t;

typedef struct directlight_s {
    struct directlight_s *next;
    emittype_t type;

    union {
        struct {
            vec3_t normal;
            winding_t *winding;
        } surface;
        struct {
            vec3_t normal;
            float dot;
        } spotlight;
        struct {
            vec3_t normal;
        } sky;
    };

    float intensity;
    int32_t style;
    float wait;
    float adjangle;
    int32_t falloff;
    vec3_t origin;
    vec3_t color;
    //vec3_t normal; // for surfaces and spotlights
    //float stopdot; // for spotlights
    dplane_t *plane;
    dleaf_t *leaf;
    dleaf_tx *leafX;
    int32_t nodenum;
} directlight_t;

// the sum of all tranfer->transfer values for a given patch
// should equal exactly 0x10000, showing that all radiance
// reaches other patches
typedef struct
{
    uint16_t patch;
    uint16_t transfer;
} transfer_t;

#define MAX_PATCHES             65535
#define MAX_PATCHES_QBSP        4000000 // qb: extended limit

#define LMSTEP                  16 // qb: lightmap step.  Default generates 1/16 of texture scale.
#define QBSP_LMSTEP             4  // qb: higher res lightmap

#define DEFAULT_SMOOTHING_VALUE 44.0
#define DEFAULT_NUDGE_VALUE     0.25

typedef struct patch_s {
    winding_t *winding;
    struct patch_s *next; // next in face
    int32_t numtransfers;
    transfer_t *transfers;
    byte *trace_hit;

    int32_t nodenum;

    int32_t cluster; // for pvs checking
    vec3_t origin;
    dplane_t *plane;

    qboolean sky;

    vec3_t totallight; // accumulated by radiosity
                       // does NOT include light
                       // accounted for by direct lighting
    float area;
    int32_t faceNumber;

    // illuminance * reflectivity = radiosity
    vec3_t reflectivity;
    vec3_t baselight; // emissivity only

    // each style 0 lightmap sample in the patch will be
    // added up to get the average illuminance of the entire patch
    struct SH1 samplelight;
    int32_t samples; // for averaging direct light
} patch_t;

extern patch_t *face_patches[MAX_MAP_FACES_QBSP];
extern entity_t *face_entity[MAX_MAP_FACES_QBSP];
extern vec3_t face_offset[MAX_MAP_FACES_QBSP]; // for rotating bmodels
extern patch_t * patches;//[MAX_PATCHES_QBSP];
extern unsigned num_patches;

extern int32_t leafparents[MAX_MAP_LEAFS_QBSP];
extern int32_t nodeparents[MAX_MAP_NODES_QBSP];

extern float lightscale;

extern qboolean use_qbsp;

void MakeShadowSplits(void);

//==============================================

void BuildVisMatrix(void);
qboolean CheckVisBit(unsigned p1, unsigned p2);

//==============================================

extern float ambient, maxlight;

void LinkPlaneFaces(void);

extern float grayscale;
extern float saturation;
extern qboolean extrasamples;
extern qboolean dicepatches;
extern int32_t numbounce;
extern qboolean noblock;
extern qboolean noedgefix;

extern directlight_t *directlights[MAX_MAP_LEAFS_QBSP];

extern byte nodehit[MAX_MAP_NODES_QBSP];

void BuildLightmaps(void);

void BuildFacelights(int32_t facenum);

void FinalLightFace(int32_t facenum);
qboolean PvsForOrigin(vec3_t org, byte *pvs);

int32_t PointInNodenum(vec3_t point);
int32_t TestLine(vec3_t start, vec3_t stop);
int32_t TestLine_color(int32_t node, vec3_t start, vec3_t stop, vec3_t occluded);
int32_t TestLine_r(int32_t node, vec3_t start, vec3_t stop);

void CreateDirectLights(void);

dleaf_t *RadPointInLeaf(vec3_t point);
dleaf_tx *RadPointInLeafX(vec3_t point);

extern dplane_t backplanes[MAX_MAP_PLANES_QBSP];
extern int32_t fakeplanes; // created planes for origin offset
extern int32_t maxdata, step;

extern float subdiv;

extern float direct_scale;
extern float entity_scale;

extern qboolean sun;
extern qboolean sun_alt_color;
extern vec3_t sun_pos;
extern float sun_main;
extern float sun_ambient;
extern vec3_t sun_color;

extern float smoothing_threshold;
extern float smoothing_value;
extern float sample_nudge;
extern int32_t num_smoothing;

extern int32_t refine_amt, refine_setting;
extern int32_t PointInLeafnum(vec3_t point);
extern void MakeTnodes(dmodel_t *bm);
extern void MakePatches(void);
extern void SubdividePatches(void);
extern void PairEdges(void);
extern void CalcTextureReflectivity(void);
extern byte *dlightdata_ptr;
extern byte dlightdata_raw[MAX_MAP_LIGHTING_QBSP];

extern float sunradscale;
