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

/** \file blender/nodes/shader/pbr_prefiltering.h
 *  \ingroup nodes
 */

#include "IMB_imbuf.h"
#include "BLI_math.h"
#include <math.h>

void getVectorFromLatLong(float u, float v, float *vec);
void getLatLongFromVector(float *vect, float *u, float *v);
void getVectorFromCubemap(int x, int y, int mapwidth, float *vec, bool seamFixup);
void getCubemapFromVector(float *vec, int cubesize, float *u, float *v, bool seamFixup);
void prefilterCubemapTexel(float *texelVect, struct ImBuf *ibuf, float roughness, int nbrSample, float *outColor);
void prefilterIMB(struct ImBuf *ibuf, struct Image *ima);
bool convertEquirectToCubemap(struct Image *ima, struct ImBuf *ibuf_ima);

