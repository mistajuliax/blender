/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version. 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2005 Blender Foundation.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/nodes/shader/pbr_prefiltering.c
 *  \ingroup nodes
 */

#include "node_shader_util.h"
#include "pbr_prefiltering.h"

#include "BLI_math.h"
#include "BLI_threads.h"

#include "BKE_context.h"
#include "BKE_report.h"
 
#include "../../editors/include/ED_screen.h"
#include "../../editors/include/ED_node.h"
 
#include "../../imbuf/intern/IMB_filter.h"

#include "GPU_draw.h"

#include "RE_render_ext.h"

#include "../../windowmanager/WM_api.h"
#include "../../windowmanager/WM_types.h"

/* ----------- Threaded Filtering -------------*/
#if 0
typedef struct FilteringQueue {
	SpinLock spin;
	ImBuf *ibuf;
	int miplevel;
	int x, y;
} FilteringQueue;
#endif

typedef struct FilteringThread {
	ImBuf *ibuf;
	int threadNbr;
	int threadTot;
} FilteringThread;

static float getRoughnessFromMip(int miplevel, int miptot)
{
	const float mipScale = 1.2f;
	const float mipOffset = 1.0f;
	return pow(2, (float)(miplevel - (miptot - 1) + mipOffset) / mipScale );
}

static unsigned int getSampleNbrFromMip(int miplevel, int UNUSED(miptot))
{
	/* Variable Sample to reduce compute time :
	 * - Start at 64 samples for low roughness
	 * - 2^13 is sufficient for low res mipmaps */
	unsigned int sampleNbr = (unsigned int)1 << (6 + miplevel*2); /* */
	CLAMP_MAX(sampleNbr, (1 << 13) ); 
	return sampleNbr;
}

static void *do_filtering_thread(void *data_v)
{
	FilteringThread *data = (FilteringThread *) data_v;
	float texelVect[3] = {0.0f}, roughness;
	int miplvl;
	unsigned int sampleNbr;
	ImBuf *ibuf = data->ibuf;

	/* Start with filtering the 1st mipmap not the 0 base level */
	for (miplvl = 1; miplvl < ibuf->miptot; ++miplvl)
	{
		ImBuf *mipbuf = IMB_getmipmap(ibuf, miplvl);
		int x = 0, y = 0;

		roughness = getRoughnessFromMip(miplvl, ibuf->miptot);
		sampleNbr = getSampleNbrFromMip(miplvl, ibuf->miptot);

		/* Each thread compute the rows that are multiple of the thread number */
		for (y = data->threadNbr; y < mipbuf->y; y += data->threadTot)
		{
			for (x = 0; x < mipbuf->x; ++x)
			{
				int channels = 4;
				float *pixel;

				if (ibuf->channels != 0)
					channels = ibuf->channels;

				getVectorFromCubemap(x, y, mipbuf->x, &texelVect, true);

				pixel = mipbuf->rect_float + channels * (y*mipbuf->x + x);

				prefilterCubemapTexel(&texelVect, ibuf, roughness, sampleNbr, pixel);
			}
		}
	}

	return NULL;
}


static void do_threaded_filtering(ImBuf *ibuf)
{
	FilteringThread *handles;
	ListBase threads;
	int i, tot_thread = BLI_system_thread_count();

#if 0
	FilteringQueue queue;
	/* Initializing queue */
	BLI_spin_init(&queue.spin);

	queue.ibuf = ibuf;
	queue.miplevel = 1; 
	queue.x = 0;
	queue.y = 0;
#endif

	handles = MEM_callocN(sizeof(FilteringThread) * tot_thread, "Filtering threaded handles");

	if (tot_thread > 1)
		BLI_init_threads(&threads, do_filtering_thread, tot_thread);

	for (i = 0; i < tot_thread; i++) {
		FilteringThread *handle = &handles[i];

		handle->threadNbr = i;
		handle->threadTot = tot_thread;
		handle->ibuf = ibuf;

		if (tot_thread > 1)
			BLI_insert_thread(&threads, handle);
	}

	if (tot_thread > 1)
		BLI_end_threads(&threads);
	else
		do_filtering_thread(handles);

	MEM_freeN(handles);
}


/* ------------- PBR filtering functions ---------------- */

static void getHammersleyPoint( unsigned int i, const unsigned int nbrSample, float *hamPoint )
{
	float phi;

	hamPoint[0] = (float)i / (float)nbrSample ;

	/* From http://holger.dammertz.org/stuff/notes_HammersleyOnHemisphere.html
	 * Radical Inverse : Van der Corput */
	i = ((i >> 16) | (i << 16));
    i = (((i & 0xaaaaaaaa) >> 1) | ((i & 0x55555555) << 1));
    i = (((i & 0xcccccccc) >> 2) | ((i & 0x33333333) << 2));
    i = (((i & 0xf0f0f0f0) >> 4) | ((i & 0x0f0f0f0f) << 4));
    i = (((i & 0xff00ff00) >> 8) | ((i & 0x00ff00ff) << 8));

	hamPoint[1] =  2.3283064365386963e-10 * (float)(i); /* 0x100000000 */

	phi = 2.0f * M_PI * hamPoint[1];
	hamPoint[2] = cos(phi);
	hamPoint[3] = sin(phi);
}

void getVectorFromLatLong(float u, float v, float *vec)
{
    vec[0] = (float)( -sin(v*M_PI) * cos(u*2.0f*M_PI) ); /* blender X */
 	vec[2] = (float)( sin(v*M_PI) * sin(u*2.0f*M_PI) ); /* blender Y */
 	vec[1] = (float)( -cos(v*M_PI) ); /* blender Z */
}

void getLatLongFromVector(float *vect, float *u, float *v)
{
	*u = (atan2(-vect[2], vect[0]) + M_PI) / (2.0f*M_PI);
	*v = atan2(vect[1], hypot(vect[0], -vect[2])) / M_PI + 0.5f;
}

void getVectorFromCubemap(int x, int y, int mapwidth, float *vec, bool seamFixup)
{
	float cubevect[3];
	int cubesize = mapwidth/3;
	int cubeface;
	float u = (float)(x%cubesize);
	float v = (float)(y%cubesize);

	if (x >= 2*cubesize)
		cubeface = 2;
	else if (x >= cubesize)
		cubeface = 1;
	else
		cubeface = 0;

	if (y < cubesize)
		cubeface += 3;

	if (seamFixup) {
		/* Code from Nvtt : http://code.google.com/p/nvidia-texture-tools/source/browse/trunk/src/nvtt/CubeSurface.cpp		
		 * transform from [0..res - 1] to [-1 .. 1], match up edges exactly. */
		u = (2.0f * u / ((float)cubesize - 1.0f) ) - 1.0f;
		v = (2.0f * v / ((float)cubesize - 1.0f) ) - 1.0f;
	}
	else {
		/* transform from [0..res - 1] to [- (1 - 1 / res) .. (1 - 1 / res)]
		 * (+ 0.5f is for texel center addressing) */
		u = (2.0f * (u + 0.5f) / (float)(cubesize) ) - 1.0f;
		v = (2.0f * (v + 0.5f) / (float)(cubesize) ) - 1.0f;
	}


	if (cubeface == 0) {
		float face_vect[3] = {u, v, -1.0f };
		copy_v3_v3(cubevect, face_vect);
	}
	else if (cubeface == 1) {
		float face_vect[3] = {-1.0f, v, -u};
		copy_v3_v3(cubevect, face_vect);
	}
	else if (cubeface == 2) {
		float face_vect[3] = {-u, v, 1.0f};
		copy_v3_v3(cubevect, face_vect);
	}
	else if (cubeface == 3) {
		float face_vect[3] = {v, -1.0f, u};
		copy_v3_v3(cubevect, face_vect);
	}
	else if (cubeface == 4) {
		float face_vect[3] = {-v, 1.0f, u};
		copy_v3_v3(cubevect, face_vect);
	}
	else {
		float face_vect[3] = {1.0f, v, u};
		copy_v3_v3(cubevect, face_vect);
	}

	normalize_v3(cubevect);
	copy_v3_v3(vec, cubevect);
}

void getCubemapFromVector(float *vec, int cubesize, float *u, float *v, bool seamFixup)
{
	int cubeface;
	float uv[2];
	float cubefacebiasU;
	float cubefacebiasV;
	float absX = fabs(vec[0]);
	float absY = fabs(vec[2]);
	float absZ = fabs(vec[1]);
	float maxCoord;


	if (absX >= absY && absX >= absZ) {
		maxCoord = absX;
		if (vec[0] > 0)
			cubeface = 5;
		else
			cubeface = 1;
	}
	else if (absY >= absX && absY >= absZ) {
		maxCoord = absY;
		if (vec[2] > 0)
			cubeface = 2;
		else
			cubeface = 0;
	}
	else {
		maxCoord = absZ;
		if (vec[1] > 0)
			cubeface = 4;
		else
			cubeface = 3;
	}

	/* divide through by max coord so face vector lies on cube face */
	mul_v3_fl(vec,1.0f/maxCoord);

	/* TODO : use an array with all face coord, would be much readable */
	if (cubeface == 0) {
		uv[0] = vec[0];
		uv[1] = vec[1];
		cubefacebiasU = 0.0f;
		cubefacebiasV = (float)cubesize;
	}
	else if (cubeface == 1) {
		uv[0] = -vec[2];
		uv[1] = vec[1];
		cubefacebiasU = (float)cubesize;
		cubefacebiasV = (float)cubesize;
	}
	else if (cubeface == 2) {
		uv[0] = -vec[0];
		uv[1] = vec[1];
		cubefacebiasU = (float)cubesize*2;
		cubefacebiasV = (float)cubesize;
	}
	else if (cubeface == 3) {
		uv[0] = vec[2];
		uv[1] = vec[0];
		cubefacebiasU = 0.0f;
		cubefacebiasV = 0.0f;
	}
	else if (cubeface == 4) {
		uv[0] = vec[2];
		uv[1] = -vec[0];
		cubefacebiasU = (float)cubesize;
		cubefacebiasV = 0.0f;
	}
	else {
		uv[0] = vec[2];
		uv[1] = vec[1];
		cubefacebiasU = (float)cubesize*2;
		cubefacebiasV = 0.0f;
	}

	if (seamFixup) {
		*u = ceil((uv[0] + 1.0f) * ((float)cubesize - 1.0f) * 0.5f );
		*v = ceil((uv[1] + 1.0f) * ((float)cubesize - 1.0f) * 0.5f );
	}
	else {
		*u = ceil((uv[0] + 1.0f) * ((float)cubesize) * 0.5f - 0.5f);
		*v = ceil((uv[1] + 1.0f) * ((float)cubesize) * 0.5f - 0.5f);
	}

	*u += cubefacebiasU;
	*v += cubefacebiasV;
	*u /= (float)cubesize*3.0f - 1.0f;
	*v /= (float)cubesize*2.0f - 1.0f;
}

/* Adapted from http://blog.selfshadow.com/publications/s2013-shading-course/karis/s2013_pbs_epic_notes_v2.pdf */
MINLINE void importance_sample_GGX(float Xi[2], float wT[3], float wB[3], float wN[3], float Hn[3], float a2)
{
	float tmp[3], H[3];

	float CosTheta = sqrt( (1 - Xi[0]) / ( 1 + (a2 - 1) * Xi[0] ) );
	float SinTheta = sqrt( 1 - CosTheta * CosTheta );

	H[0] = SinTheta* Xi[2];
	H[1] = SinTheta* Xi[3];
	H[2] = CosTheta;

	// Tangent to world space
	mul_v3_v3fl(Hn, wT, H[0]);
	mul_v3_v3fl(tmp, wB, H[1]);
	add_v3_v3(Hn, tmp);
	mul_v3_v3fl(tmp, wN, H[2]);
	add_v3_v3(Hn, tmp);
}

/* Adapted from http://blog.selfshadow.com/publications/s2013-shading-course/karis/s2013_pbs_epic_notes_v2.pdf */
void prefilterCubemapTexel(float *texelVect, struct ImBuf *ibuf, float roughness, int nbrSample, float *outColor)
{

	float FilteredColor[3] = {0.0f};
	float UpVector[3] = {0.0f};
	float sample[4];
	float Weight = 0.0f;
	float wT[3], wB[3], wN[3], Xi[2], Hn[3], Ln[3];
	float u, v, a, a2, nl;
	float dot;
	int i;

	a = roughness*roughness;
	CLAMP_MIN(a, 1e-4f);
	a2 = a*a;

	copy_v3_v3(wN, texelVect);

	/* Computing Tangent space */
	if (fabs(wN[1]) < 0.99999)
		UpVector[1] = 1.0f;
	else
		UpVector[0] = 1.0f;
	cross_v3_v3v3(wT,UpVector,wN);
	normalize_v3(wT);
	cross_v3_v3v3(wB,wN,wT);

	/* Importance sampling */
	for (i = 0; i < nbrSample; i++) {

		getHammersleyPoint(i, nbrSample, &Xi);
		importance_sample_GGX(Xi, wT, wB, wN, Hn, a2);

		/* Reflect */
		dot = dot_v3v3(wN, Hn);
		mul_v3_v3fl(Ln, Hn, 2*dot);
		sub_v3_v3(Ln, wN);

		nl = dot_v3v3(wN, Ln);
		CLAMP(nl, 0, 1);

		if (nl > 0) {
			getCubemapFromVector(&Ln, ibuf->y/2, &u, &v, true);
			ibuf_sample(ibuf, u, v, (1.0f / (float)ibuf->x), (1.0f / (float)ibuf->y), sample);

			mul_v3_fl(&sample, nl);
			add_v3_v3(FilteredColor, &sample);

			Weight += nl;
		}
	}

	CLAMP_MIN(Weight, 0.001f);
	mul_v3_v3fl(outColor, FilteredColor, 1.0f / Weight);
}

void prefilterIMB(struct ImBuf *ibuf, struct Image *ima)
{
	int miplevel = 0;
	bool do_filtering = false;
	int i;

	if (ibuf->mipmap[0] == NULL) {
		IMB_makemipmap(ibuf, 0);
		do_filtering = true;
	}
	else {
		IMB_remakemipmap(ibuf, 0);
		ibuf->userflags &= ~IB_MIPMAP_INVALID;
		do_filtering = true;
	}

	if (do_filtering && ibuf->miptot > 2) {

		/* Force gpu texture update */
		GPU_free_image(ima);

		/* Delete the last 2 mip 2x1 is too small for a cubemap 
		 * and 3x2 does not containt enough data */
		for (i = 0; i < 2; ++i)
		{
			miplevel = ibuf->miptot - 1;
			if (ibuf->mipmap[miplevel])
				IMB_freeImBuf(ibuf->mipmap[miplevel]);
			ibuf->mipmap[miplevel] = NULL;
			ibuf->miptot--;
		}

		do_threaded_filtering(ibuf);
	}
}

/* ----------------- Conversion functions --------------------- */

bool convertEquirectToCubemap(struct Image *ima, struct ImBuf *ibuf_ima)
{
	ImBuf *ibuf;
	unsigned int x, y;

	/* Force re-filtering and updating the gpu texture */
	GPU_free_image(ima);
	ibuf_ima->userflags |= IB_MIPMAP_INVALID;

	/* Cubemap is an ibuf in the original ibuf */
	if (!ibuf_ima->cubemap_ibuf) {
		/* Size of the cubemap is the nearest power of 2 from half the height of the Equirectangular map */
		unsigned int cubesize = (unsigned int)1 << (int)round(log(ibuf_ima->y/2)/log(2));

		ibuf_ima->cubemap_ibuf = IMB_allocImBuf(cubesize*3, cubesize*2, 32, IB_rect | IB_rectfloat);

		if (!ibuf_ima->cubemap_ibuf)
			return false;
		
		IMB_refImBuf(ibuf_ima->cubemap_ibuf);
	}

	ibuf = ibuf_ima->cubemap_ibuf;

	/* Create Cubemap */
	for (y = 0; y < ibuf->y; ++y)
	{
		for (x = 0; x < ibuf->x; ++x)
		{
			float *pixel, u, v;
			float vec[3] = {0.0f};

			pixel = ibuf->rect_float + ibuf->channels * (y*ibuf->x + x);

			getVectorFromCubemap(x, y, ibuf->x, &vec, false);
			getLatLongFromVector(&vec, &u, &v);
			ibuf_sample(ibuf_ima, u, v, (1.0f / (float)ibuf_ima->x), (1.0f / (float)ibuf_ima->y), pixel);
		}
	}

	IMB_rect_from_float(ibuf);

	return true;
}

/* ********************** Shader Node Button ********************** */

static int node_shader_prefilter_envmap_poll(bContext *C)
{
	SpaceNode *snode = CTX_wm_space_node(C);
	bNode *node;

	/* see if we have a shader script node in context */
	node = CTX_data_pointer_get_type(C, "node", &RNA_ShaderNodePbrShader).data;

	if (!node && snode && snode->edittree)
		node = nodeGetActive(snode->edittree);

	if (node && node->type == SH_NODE_PBR_SHADER) {
		NodePbrShader *pbrNode = node->storage;

		if (node->id) {
			return ED_operator_node_editable(C);
		}
	}

	return 0;
}

static int node_shader_prefilter_envmap_exec(bContext *C, wmOperator *op)
{
	Main *bmain = CTX_data_main(C);
	Scene *scene = CTX_data_scene(C);
	SpaceNode *snode = CTX_wm_space_node(C);
	PointerRNA nodeptr = CTX_data_pointer_get_type(C, "node", &RNA_ShaderNodePbrShader);
	bNodeTree *ntree = NULL;
	bNode *node = NULL;
	int result = 0;

	/* get node */
	if (nodeptr.data) {
		ntree = nodeptr.id.data;
		node = nodeptr.data;
	}
	else if (snode && snode->edittree) {
		ntree = snode->edittree;
		node = nodeGetActive(snode->edittree);
	}

	if (node && node->type == SH_NODE_PBR_SHADER) {
		NodePbrShader *pbrNode = node->storage;

		if (node->id) {
			Image *ima = (Image *)node->id;
			ImageUser *iuser = NULL;
			ImBuf *ibuf = BKE_image_acquire_ibuf(ima, iuser, NULL), *ibuf_ima = NULL;

			if (!ibuf || !ima || !ima->envmap_ok) {
				if (ibuf && ima) 
					BKE_image_release_ibuf(ima, ibuf, NULL);
				BKE_report(op->reports, RPT_ERROR, "Envmap is not ready");
				return OPERATOR_CANCELLED;
			}

			if (ima->pbr_projection != SHD_PROJ_CUBEMAP) {
				if (ibuf->cubemap_ibuf) {
					ibuf_ima = ibuf;
					ibuf = ibuf->cubemap_ibuf;
				}
				else {
					BKE_image_release_ibuf(ima, ibuf, NULL);
					BKE_report(op->reports, RPT_ERROR, "Envmap was not converted");
					return OPERATOR_CANCELLED;
				}
			}


			if (ibuf->y*1.5f == ibuf->x) {
				prefilterIMB(ibuf, ima);
				BKE_report(op->reports, RPT_INFO, "Filtering done!");
				result = true;
			}
			else {
				BKE_report(op->reports, RPT_ERROR, "Envmap is not a cubemap");
			}

			if (ibuf_ima)
				BKE_image_release_ibuf(ima, ibuf_ima, NULL);
			else
				BKE_image_release_ibuf(ima, ibuf, NULL);
		}

		ED_node_tag_update_nodetree(bmain, ntree);
	}


	return (result) ? OPERATOR_FINISHED : OPERATOR_CANCELLED;
}

void NODE_OT_prefilter_envmap(wmOperatorType *ot)
{
	/* identifiers */
	ot->name = "Prefilter Envmap";
	ot->description = "Precompute roughness level into mipmaps.";
	ot->idname = "NODE_OT_prefilter_envmap";

	/* api callbacks */
	ot->exec = node_shader_prefilter_envmap_exec;
	// ot->invoke = objects_bake_render_invoke;
	// ot->modal = objects_bake_render_modal;
	ot->poll = node_shader_prefilter_envmap_poll;
}
