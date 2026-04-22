/*!
 * \copyright   This file is part of the public version of the AREPO code.
 * \copyright   Copyright (C) 2009-2019, Max-Planck Institute for Astrophysics
 * \copyright   Developed by Volker Springel (vspringel@MPA-Garching.MPG.DE) and
 *              contributing authors.
 * \copyright   Arepo is free software: you can redistribute it and/or modify
 *              it under the terms of the GNU General Public License as published by
 *              the Free Software Foundation, either version 3 of the License, or
 *              (at your option) any later version.
 *
 *              Arepo is distributed in the hope that it will be useful,
 *              but WITHOUT ANY WARRANTY; without even the implied warranty of
 *              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *              GNU General Public License for more details.
 *
 *              A copy of the GNU General Public License is available under
 *              LICENSE as part of this program.  See also
 *              <https://www.gnu.org/licenses/>.
 *
 * \file        src/fof/fof_seeding_registry.c
 * \date        20/10/2025
 * \details     contains functions:
 *                      void fof_seeding_init(int RestartFlag)
 *                      void fof_seeding_io(int modus)   
 *                      void halo_seed_registry_init(HaloSeedRegistry *r, int restart_flag)    
 *                      void halo_seed_registry_free(HaloSeedRegistry *r)
 *                      void halo_seed_registry_grow(HaloSeedRegistry *r)  
 *                      int  halo_is_seeded(HaloSeedRegistry *r, MyIDType id)
 *                      void halo_mark_seeded(HaloSeedRegistry *r, MyIDType id)
 *          
 */

#include <gsl/gsl_math.h>
#include <inttypes.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "fof_seeding.h"

HaloSeedRegistry HaloSeeds;

typedef enum {
    RESTART_READ  = 0,
    RESTART_WRITE = 1
} RestartMode;

void fof_seeding_init(int RestartFlag)
{
    halo_seed_registry_init(&HaloSeeds, RestartFlag);
}

void halo_seed_registry_init(HaloSeedRegistry *r, int RestartFlag)
{
    r->n   = 0;
    r->max = 0;

    if(RestartFlag == 0)
    {
        r->max = 10000;
        r->ids = mymalloc_movable(&r->ids,
                                  "HaloSeedIDs",
                                  r->max * sizeof(MyIDType));
    }
    else
        r->ids = NULL;
}

void halo_seed_registry_free(HaloSeedRegistry *r)
{
    if(r->ids)
        myfree_movable(r->ids);

    r->ids = NULL;
    r->n = r->max = 0;
}

void halo_seed_registry_grow(HaloSeedRegistry *r)
{
    if(r->ids == NULL || r->max==0)
    {
        r->max = 128;

        r->ids = mymalloc_movable(&r->ids,
                                 "HaloSeedIDs",
                                 r->max * sizeof(MyIDType));
        return;
    }

    r->max = 2 * r->max + 64;

    myrealloc_movable(&r->ids,
                      r->max * sizeof(MyIDType));
}

int halo_is_seeded(HaloSeedRegistry *r, MyIDType id)
{
    int lo = 0, hi = r->n - 1;

    while(lo <= hi)
    {
        int mid = (lo + hi) / 2;

        if(r->ids[mid] == id)
            return 1;
        else if(r->ids[mid] < id)
            lo = mid + 1;
        else
            hi = mid - 1;
    }

    return 0;
}

void halo_mark_seeded(HaloSeedRegistry *r, MyIDType id)
{
    if(r->n == r->max)
        halo_seed_registry_grow(r);

    int pos = r->n;

    while(pos > 0 && r->ids[pos-1] > id)
    {
        r->ids[pos] = r->ids[pos-1];
        pos--;
    }

    r->ids[pos] = id;
    r->n++;
}

void halo_seed_registry_pack(HaloSeedRegistry *r, MyIDType **buf, int *n)
{
    *buf = r->ids;
    *n   = r->n;
}

void halo_seed_registry_unpack(HaloSeedRegistry *r, MyIDType *buf, int n)
{
     if(r->ids)
        myfree_movable(r->ids);

    r->n = n;
    r->max = n;

    r->ids = mymalloc_movable(&r->ids,
                              "HaloSeedIDs",
                              n * sizeof(MyIDType));

    memcpy(r->ids, buf, n * sizeof(MyIDType));
}