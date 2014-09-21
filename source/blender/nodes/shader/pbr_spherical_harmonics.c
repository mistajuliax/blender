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

/** \file blender/nodes/shader/pbr_spherical_harmonics.c
 *  \ingroup nodes
 */

#include "node_shader_util.h"
#include "pbr_prefiltering.h"
#include "pbr_spherical_harmonics.h"

#include "BLI_math.h"

/* Adapted from sebastien lagarde's ATI Cubemapgen */ 

static float AreaElement( float x, float y )
{
	return atan2(x * y, sqrt(x * x + y * y + 1));
}

static void getTexelVectorSolidAngle(int x, int y, int width, int height, float *texelVect)
{
	float u, v, x0, y0, x1, y1;
	float halfTexelSize[2];
	int face_width = height/2;

	u = (float)(x%face_width);
	v = (float)(y%face_width);

	getVectorFromCubemap(x, y, width, texelVect, false);

	/* transform from [0..res - 1] to [- (1 - 1 / res) .. (1 - 1 / res)]
	 (+ 0.5f is for texel center addressing) */
	u = (2.0f * (u + 0.5f) / ((float)(face_width)) ) - 1.0f;
	v = (2.0f * (v + 0.5f) / ((float)(face_width)) ) - 1.0f;

 	/* Solid angle weight approximation :
	 * U and V are the -1..1 texture coordinate on the current face.
	 * Get projected area for this texel */
	halfTexelSize[0] = 1.0f / (float)face_width;
	x0 = u - halfTexelSize[0];
	y0 = v - halfTexelSize[0];
	x1 = u + halfTexelSize[0];
	y1 = v + halfTexelSize[0];

	texelVect[3] = AreaElement(x0, y0) - AreaElement(x0, y1) - AreaElement(x1, y0) + AreaElement(x1, y1);
}

void SHFilter(ImBuf *ibuf, float *SHCoef)
{
	const int NUM_SH_COEFFICIENT = 9;
	int channels = 4, ofs;
	char *rect;
	float texelVect[4] = {0.0f};
	float SHdir[9] = {0.0f};
	float weightAccum = 0.0f;
	float weight = 0.0f;
	int x,y,i;
	float r,g,b, xV,yV,zV;

	if (ibuf->channels != 0)
		channels = ibuf->channels;

	/* Tets if valid cubemap */
	if (3*ibuf->y/2 != ibuf->x) {
		printf("SH: Incorrect envmap size\n");
		return;
	}

	for (y = 0; y < ibuf->y; y++)
	{
		for (x = 0; x < ibuf->x; x++)
		{
			ofs = y*ibuf->x + x;

			getTexelVectorSolidAngle(x, y, ibuf->x, ibuf->y, &texelVect);

			weight = texelVect[3];   

			r = xV = texelVect[0];
			g = yV = texelVect[1];
			b = zV = texelVect[2];

		    SHdir[0] = (float)(0.282095f);

		    SHdir[1] = (float)(-0.488603f * yV * 2.0f/3.0f);
		    SHdir[2] = (float)(0.488603f * zV * 2.0f/3.0f);
		    SHdir[3] = (float)(-0.488603f * xV * 2.0f/3.0f);

		    SHdir[4] = (float)(1.092548f * xV * yV * 1.0f/4.0f);
		    SHdir[5] = (float)(-1.092548f * yV * zV * 1.0f/4.0f);
		    SHdir[6] = (float)(0.315392f * (3.0f * zV * zV - 1.0f) * 1.0f/4.0f);
		    SHdir[7] = (float)(-1.092548f * xV * zV * 1.0f/4.0f);
		    SHdir[8] = (float)(0.546274f * (xV * xV - yV * yV) * 1.0f/4.0f);
		    
			if (ibuf->rect_float) {
				if (ibuf->channels > 2) {
					r = ibuf->rect_float[channels*ofs+0];
					g = ibuf->rect_float[channels*ofs+1];
					b = ibuf->rect_float[channels*ofs+2];
				}
				else {
					r = g = b = ibuf->rect_float[ofs];
				}
			}else{
				rect = (char *)(ibuf->rect + ofs);
				r = ((float)rect[0])*(1.0f/255.0f);
				g = ((float)rect[1])*(1.0f/255.0f);
				b = ((float)rect[2])*(1.0f/255.0f);
			}

			for (i = 0; i < NUM_SH_COEFFICIENT; i++)
			{
				SHCoef[i*3+0] += r * SHdir[i] * weight;
				SHCoef[i*3+1] += g * SHdir[i] * weight;
				SHCoef[i*3+2] += b * SHdir[i] * weight;
			}

			weightAccum += weight;
		}
	}
	
	/* Normalization - The sum of solid angle should be equal to the solid angle of the sphere (4 PI), so
	 * normalize in order our weightAccum exactly match 4 PI. */

	for (i = 0; i < NUM_SH_COEFFICIENT; ++i)
	{
	 	SHCoef[i*3+0] *= 4.0f * M_PI / weightAccum;
	 	SHCoef[i*3+1] *= 4.0f * M_PI / weightAccum;
	 	SHCoef[i*3+2] *= 4.0f * M_PI / weightAccum;
	}
}
