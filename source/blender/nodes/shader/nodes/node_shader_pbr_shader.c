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

#include "../node_shader_util.h"
#include "../pbr_prefiltering.h"
#include "../pbr_spherical_harmonics.h"

#include "IMB_colormanagement.h"
#include "BLI_math.h"
#include "GPU_draw.h"
#include "../../imbuf/intern/IMB_filter.h"

#include <math.h>

/* **************** OUTPUT ******************** */

static bNodeSocketTemplate sh_node_pbr_shader_in[] = {
	{ 	SOCK_RGBA, 1, N_("Albedo"),			  0.8f, 0.8f, 0.8f, 1.0f, 0.0f, 1.0f},
	{ 	SOCK_RGBA, 1, N_("Specular"),	      0.04f, 0.04f, 0.04f, 1.0f, 0.0f, 1.0f, PROP_NONE},
	{ 	SOCK_FLOAT, 1, N_("Roughness"),		  0.3f, 0.3f, 0.3f, 0.3f, 0.0f, 1.0f, PROP_NONE},
	{	SOCK_VECTOR, 1, N_("Normal"),	      0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, PROP_NONE, SOCK_HIDE_VALUE},
	{ 	SOCK_RGBA, 1, N_("Occlusion"),		  1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f, PROP_NONE},
	{ 	SOCK_VECTOR, 1, N_("EnvMap Mapping"), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 360.0f, PROP_NONE, SOCK_HIDE_VALUE},
	{ 	SOCK_FLOAT, 1, N_("EnvMap Horizon Fade"), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 4.0f, PROP_NONE},
	{ 	SOCK_RGBA, 1, N_("EnvMap Fresnel Intensity"), 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f, PROP_NONE},
	{ 	SOCK_RGBA, 1, N_("EnvMap Intensity"), 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 10.0f, PROP_NONE},
	{	-1, 0, ""	}
};

#define SH_PBR_INPUT_ALBEDO 	0
#define SH_PBR_INPUT_SPECULAR 	1
#define SH_PBR_INPUT_ROUGHNESS 	2
#define SH_PBR_INPUT_NORMAL 	3
#define SH_PBR_INPUT_OCCLUSION 	4
#define SH_PBR_INPUT_IBL_ROTAT 	5
#define SH_PBR_INPUT_IBL_FADE 	6
#define SH_PBR_INPUT_IBL_FRESN  7
#define SH_PBR_INPUT_IBL_POWER 	8

static bNodeSocketTemplate sh_node_pbr_shader_out[] = {
	{	SOCK_RGBA, 0, N_("Combined")},
	{	SOCK_RGBA, 0, N_("SH")},
	{	SOCK_RGBA, 0, N_("IBL Spec")},
	{	SOCK_RGBA, 0, N_("Lighting")},
	{	-1, 0, ""	}
};

#define SH_PBR_OUTPUT_COMBINED 		0
#define SH_PBR_OUTPUT_SH 			1
#define SH_PBR_OUTPUT_IBL 			2
#define SH_PBR_OUTPUT_LIGHTING 		3

static void node_shader_init_pbr_shader(bNodeTree *UNUSED(ntree), bNode *node)
{
	NodePbrShader *pbr = MEM_callocN(sizeof(NodePbrShader), "NodePbrShader");
	default_tex_mapping(&pbr->base.tex_mapping, TEXMAP_TYPE_POINT);
	default_color_mapping(&pbr->base.color_mapping);
	pbr->iuser.frames = 1;
	pbr->iuser.sfra = 1;
	pbr->iuser.fie_ima = 2;
	pbr->iuser.ok = 1;
	pbr->idMat = NULL;
	pbr->energy_conservation = 1;
	node->storage = pbr;
	node->custom1 = SH_NODE_MAT_DIFF | SH_NODE_MAT_SPEC;
}

static GPUNodeLink *gpu_get_input_link(GPUNodeStack *in)
{
	if (in->link)
		return in->link;
	else
		return GPU_uniform(in->vec);
}

static int node_shader_gpu_pbr_shader(GPUMaterial *gpumat, bNode *node, bNodeExecData *UNUSED(execdata), GPUNodeStack *in, GPUNodeStack *out)
{
	NodePbrShader *PbrNode = node->storage;
	Image *ima = (Image *)node->id;
	bNodeSocket *sock;
	ImageUser *iuser = NULL;
	GPUNodeLink *outSh, *outIbl, *outLight, *outCombined;
	int isdata = false;
	int ret = 0, i = 0;
	float scale = 511.0f, one = 1.0f, maxLod = 0.0f, cubeSize = 0.0f, cubeSizeSquared = 0.0f;
	float shCoef[9][3] = { { 0.0f } };
	float rotMat[4][4];
	bool ibl_ok = true;
	bool use_importance_sample = (ima != NULL && ima->envmap_sampling != SHD_PBR_SAMPLE_PRECALC);

	/* Envmap Rotation : Copying matrix from mapping node if one is connected */
	unit_m4(rotMat);
	for (sock = node->inputs.first; sock; sock = sock->next){
		if(i == SH_PBR_INPUT_IBL_ROTAT && sock->link){
			bNodeLink *mappingLink = sock->link;
			if(mappingLink->fromnode){
				bNode *mappingNode = mappingLink->fromnode;
				if(mappingNode->typeinfo->type == SH_NODE_MAPPING){
					TexMapping *texmap = mappingNode->storage;
					memcpy(rotMat, texmap->mat, 4 * 4 * sizeof(float));
				}
				break;
			}
		}else{
			i++;
		}
	}

	/* Energy Conservation */
	if (PbrNode->energy_conservation) {
		GPU_link(gpumat, "mix_sub", GPU_uniform(&one), gpu_get_input_link(&in[SH_PBR_INPUT_ALBEDO]), gpu_get_input_link(&in[SH_PBR_INPUT_SPECULAR]), &in[SH_PBR_INPUT_ALBEDO].link);
		GPU_link(gpumat, "rgba_clamp", in[SH_PBR_INPUT_ALBEDO].link, &in[SH_PBR_INPUT_ALBEDO].link);
	}

	/* Use world normal by default */
	if (!in[SH_PBR_INPUT_NORMAL].link){
		GPU_link(gpumat, "world_normal", GPU_builtin(GPU_VIEW_NORMAL), GPU_builtin(GPU_INVERSE_VIEW_MATRIX), &in[SH_PBR_INPUT_NORMAL].link);
	}

	/* Conversion to cubemap, Spherical harmonics */
	if (ima) {
		ImBuf *ibuf = BKE_image_acquire_ibuf(ima, iuser, NULL);
		ImBuf *ibuf_ima = NULL;

		if (ibuf) {
			bool projection_has_changed = (ima->pbr_last_projection && (ima->pbr_last_projection != ima->pbr_projection + 1) );
			bool sampling_has_changed = (ima->envmap_ok != ima->envmap_sampling + 1);
			
			/* Converting to cubemap */
			if (ima->pbr_projection != SHD_PROJ_CUBEMAP) {

				if (!ibuf->cubemap_ibuf || projection_has_changed) {
					ibl_ok = convertEquirectToCubemap(ima, ibuf);
				}

				if (ibl_ok) {
					ibuf_ima = ibuf;
					ibuf = ibuf->cubemap_ibuf;
				}
			}

			/* Test cubemap size */
			if (ibl_ok && (ibuf->y*1.5f == ibuf->x)) { 

				/* If no data or data is outdated */
				if (!ima->SH_Coefs[0][0] || projection_has_changed) {
					/* Compute SH & Store SH */
					SHFilter(ibuf, &shCoef);
					memcpy(ima->SH_Coefs, shCoef, 9 * 3 * sizeof(float));
				}
				else{
					/* Get from Image Data */
					memcpy(shCoef, ima->SH_Coefs, 9 * 3 * sizeof(float));
				}

				/* Determine MaxLOD based on image dimentions */
				cubeSize = (float)ibuf->x / 3.0f;
				maxLod = log(cubeSize)/log(2);
				cubeSizeSquared = cubeSize*cubeSize;

				if (!use_importance_sample)
					maxLod = ibuf->miptot;
				
				/* Force gpu texture update */
				if (projection_has_changed || sampling_has_changed)
					GPU_free_image(ima);
				
				/* Tag to not redo calculations for this type of projection */
				ima->pbr_last_projection = ima->pbr_projection + 1;
				ima->envmap_ok = ima->envmap_sampling + 1;
			}
			else {
				printf("Incorrect envmap size\n");
				ima->pbr_last_projection = 0;
				ima->envmap_ok = 0;
				ibl_ok = false;
			}

			if (ibuf_ima)
				BKE_image_release_ibuf(ima, ibuf_ima, NULL);
			else
				BKE_image_release_ibuf(ima, ibuf, NULL);
		}
	}
	
	/* Image Based Lighting : Spherical Harmonic (Diffuse lighting) */
	ret |= !GPU_link(gpumat, "irradianceFromSH", gpu_get_input_link(&in[SH_PBR_INPUT_NORMAL]), GPU_uniform(shCoef[0]), GPU_uniform(shCoef[1]), GPU_uniform(shCoef[2]), GPU_uniform(shCoef[3]), GPU_uniform(shCoef[4]), GPU_uniform(shCoef[5]), GPU_uniform(shCoef[6]), GPU_uniform(shCoef[7]), GPU_uniform(shCoef[8]), GPU_uniform(&rotMat), &outSh);
	ret |= !GPU_link(gpumat, "mix_mult", GPU_uniform(&one), gpu_get_input_link(&in[SH_PBR_INPUT_ALBEDO]), outSh, &outSh);
	ret |= !GPU_link(gpumat, "mix_mult", GPU_uniform(&one), gpu_get_input_link(&in[SH_PBR_INPUT_OCCLUSION]), outSh, &outSh);
	/* End of Spherical Harmonics */


	/* Image Based Lighting : Specular Envmap */
	if (!ima || !ibl_ok){
		float black[4] = {0.0f};
		ret |= !GPU_link(gpumat, "set_rgba", GPU_uniform(&black), &outIbl);
	}else{
		if (!ima->pbr_envmap) {
			/* TODO : remove this leak */
			ima->pbr_envmap = BKE_add_envmap();
		}
		ima->pbr_envmap->ima = ima;
		ima->pbr_envmap->stype = ENV_LOAD;
		ima->pbr_envmap->gpu_use_mips = !use_importance_sample;

		if (!use_importance_sample) 
			ret |= !GPU_link(gpumat, "node_pbr_ibl_approx_cubemap",	GPU_envmap(ima->pbr_envmap, iuser, isdata), gpu_get_input_link(&in[SH_PBR_INPUT_ROUGHNESS]),  gpu_get_input_link(&in[SH_PBR_INPUT_NORMAL]),  gpu_get_input_link(&in[SH_PBR_INPUT_IBL_FRESN]),  gpu_get_input_link(&in[SH_PBR_INPUT_SPECULAR]), gpu_get_input_link(&in[SH_PBR_INPUT_IBL_FADE]),  GPU_builtin(GPU_VIEW_NORMAL), GPU_builtin(GPU_VIEW_POSITION), GPU_builtin(GPU_INVERSE_VIEW_MATRIX), GPU_uniform(&maxLod), GPU_uniform(&cubeSize), GPU_uniform(&rotMat), &outIbl);
		else if (ima->envmap_sampling == SHD_PBR_SAMPLE_32)
			ret |= !GPU_link(gpumat, "node_pbr_ibl_32_cubemap",	GPU_envmap(ima->pbr_envmap, iuser, isdata), gpu_get_input_link(&in[SH_PBR_INPUT_ROUGHNESS]), gpu_get_input_link(&in[SH_PBR_INPUT_NORMAL]), gpu_get_input_link(&in[SH_PBR_INPUT_IBL_FRESN]),  gpu_get_input_link(&in[SH_PBR_INPUT_SPECULAR]), gpu_get_input_link(&in[SH_PBR_INPUT_IBL_FADE]),  GPU_builtin(GPU_VIEW_NORMAL), GPU_builtin(GPU_VIEW_POSITION), GPU_builtin(GPU_INVERSE_VIEW_MATRIX), GPU_uniform(&cubeSizeSquared), GPU_uniform(&rotMat), &outIbl);
		else
			ret |= !GPU_link(gpumat, "node_pbr_ibl_64_cubemap",	GPU_envmap(ima->pbr_envmap, iuser, isdata), gpu_get_input_link(&in[SH_PBR_INPUT_ROUGHNESS]), gpu_get_input_link(&in[SH_PBR_INPUT_NORMAL]), gpu_get_input_link(&in[SH_PBR_INPUT_IBL_FRESN]),  gpu_get_input_link(&in[SH_PBR_INPUT_SPECULAR]), gpu_get_input_link(&in[SH_PBR_INPUT_IBL_FADE]),  GPU_builtin(GPU_VIEW_NORMAL), GPU_builtin(GPU_VIEW_POSITION), GPU_builtin(GPU_INVERSE_VIEW_MATRIX), GPU_uniform(&cubeSizeSquared), GPU_uniform(&rotMat), &outIbl);
	}
	ret |= !GPU_link(gpumat, "mix_mult", GPU_uniform(&one), gpu_get_input_link(&in[SH_PBR_INPUT_OCCLUSION]), outIbl, &outIbl);
	/* End of Specular Envmap */


	/* Analytic Lighting (lamps)*/
	if (PbrNode->idMat) {
		Material *ma = (Material *)PbrNode->idMat;
		GPUShadeInput shi;
		GPUShadeResult shr;
	
		ma->spec_shader = MA_SPEC_GGX;
		ma->diff_shader = MA_DIFF_LAMBERT;
		
		GPU_shadeinput_set(gpumat, ma, &shi);

		/* Diffuse Color */
		ret |= !GPU_link(gpumat, "set_rgb", gpu_get_input_link(&in[0]), &shi.rgb);
		/* Specular Color/Intensity */
		ret |= !GPU_link(gpumat, "set_rgb", gpu_get_input_link(&in[1]), &shi.specrgb);
		/* Diffuse Intensity */
		ret |= !GPU_link(gpumat, "set_value", GPU_uniform(&one), &shi.refl);
		/* Roughness : Range converted from 0-1 to 1-511*/
		ret |= !GPU_link(gpumat, "set_value", gpu_get_input_link(&in[2]), &shi.har);
		ret |= !GPU_link(gpumat, "math_subtract", GPU_uniform(&one), shi.har, &shi.har);
		ret |= !GPU_link(gpumat, "math_multiply", GPU_uniform(&scale), shi.har, &shi.har);

		if (in[SH_PBR_INPUT_NORMAL].link) {
			GPUNodeLink *tmp;
			ret |= !GPU_link(gpumat, "world_to_viewspace", gpu_get_input_link(&in[SH_PBR_INPUT_NORMAL]), GPU_builtin(GPU_VIEW_MATRIX), &shi.vn);
			ret |= !GPU_link(gpumat, "vec_math_normalize", shi.vn, &shi.vn, &tmp);
		}
		
		/* Compute shading */
		GPU_shaderesult_set(&shi, &shr);

		outLight = shr.combined;
	}
	else {
		float black[4] = {0.0f};
		ret |= !GPU_link(gpumat, "set_rgba", GPU_uniform(&black), &outLight);
	}
	/* End of Analytic Lighting */

	/* Combine lights passes */
	ret |= !GPU_link(gpumat, "mix_add", GPU_uniform(&one), outIbl, outSh, &outCombined);
	ret |= !GPU_link(gpumat, "mix_mult", GPU_uniform(&one), outCombined, gpu_get_input_link(&in[SH_PBR_INPUT_IBL_POWER]), &outCombined);
	ret |= !GPU_link(gpumat, "mix_add", GPU_uniform(&one), outCombined, outLight, &outCombined);

	/* Ouput */
	out[SH_PBR_OUTPUT_COMBINED].link = outCombined;
	out[SH_PBR_OUTPUT_SH].link = outSh;
	out[SH_PBR_OUTPUT_IBL].link = outIbl;
	out[SH_PBR_OUTPUT_LIGHTING].link = outLight;


	/* Should return false if any GPU_link fail */
	return !ret;
}

/* node type definition */
void register_node_type_sh_pbr_shader(void)
{
	static bNodeType ntype;

	sh_node_type_base(&ntype, SH_NODE_PBR_SHADER, "Pbr Shader", NODE_CLASS_SHADER, 0);
	node_type_compatibility(&ntype, NODE_OLD_SHADING);
	node_type_socket_templates(&ntype, sh_node_pbr_shader_in, sh_node_pbr_shader_out);
	node_type_init(&ntype, node_shader_init_pbr_shader);
	node_type_storage(&ntype, "NodePbrShader", node_free_standard_storage, node_copy_standard_storage);
	node_type_gpu(&ntype, node_shader_gpu_pbr_shader);

	nodeRegisterType(&ntype);
}