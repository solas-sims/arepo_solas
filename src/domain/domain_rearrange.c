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
 * \file        src/domain_rearrange.c
 * \date        05/2018
 * \brief       Rearranges particle and cell arrays and gets rid of inactive
 *              particles.
 * \details     contains functions:
 *                void domain_rearrange_particle_sequence(void)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 05.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../mesh/voronoi/voronoi.h"
#include "domain.h"

/*! \brief Gets rid of inactive/eliminated cells and particles.
 *
 *  Cells that were de-refined or turned into star particles are kept in the
 *  SphP array, but flagged as inactive until this point. This routine cleans
 *  up these arrays in order to make sure only active particles/cells are
 *  exported.
 *
 *  \return void
 */
void domain_rearrange_particle_sequence(void)
{
#ifdef USE_SFR
  if(Stars_converted)
    {
      struct particle_data psave;
      peanokey key;

      for(int i = 0; i < NumGas; i++)
        if(P[i].Type != 0) /*If not a gas particle, swap to the end of the list */
          {
            psave = P[i];
            key   = Key[i];

            P[i]    = P[NumGas - 1];
            SphP[i] = SphP[NumGas - 1];
            Key[i]  = Key[NumGas - 1];

            P[NumGas - 1]   = psave;
            Key[NumGas - 1] = key;

#ifdef STARS
            SPP(NumGas - 1).PID = NumGas - 1;
#endif
            NumGas--;
            i--;
          }
      /* Now we have rearranged the particles,
       * we don't need to do it again unless there are more stars
       */
      Stars_converted = 0;
    }
#endif /* #ifdef USE_SFR */

#ifdef REFINEMENT_MERGE_CELLS
  int i, count_elim, count_gaselim;

  count_elim    = 0;
  count_gaselim = 0;

  for(i = 0; i < NumPart; i++)
    if((P[i].Mass == 0 && P[i].ID == 0) || (P[i].Type == 4 && P[i].Mass == 0))
      {
        if(P[i].Type == 0)
          {
            P[i]    = P[NumGas - 1];
            SphP[i] = SphP[NumGas - 1];
            Key[i]  = Key[NumGas - 1];

            P[NumGas - 1]   = P[NumPart - 1];
            Key[NumGas - 1] = Key[NumPart - 1];

#ifdef STARS
            if(P[NumGas - 1].Type == 4)
              SPP(NumGas - 1).PID = NumGas - 1;
#endif

#ifdef BLACKHOLES
            if(P[NumGas - 1].Type == 5)
              BPP(NumGas - 1).PID = NumGas - 1;
#endif

            NumGas--;
            count_gaselim++;
          }
        else
          {
            /* If the ELIMINATED particle (P[i], before being overwritten below) is itself a
             * black hole or star, its own BhP[]/SP[] slot must be reclaimed here -- via
             * swap-with-last on BhP[]/SP[] itself, decrementing NumBhs/NumStars -- BEFORE
             * P[i] is overwritten, since P[i].BhID/SID is the only way to find that slot.
             * Without this, NumBhs/NumStars never shrinks even though the particle is gone
             * (e.g. after bh_merger() zeroes out a merged-away BH's Mass/ID expecting this
             * function to clean it up), leaving an orphaned BhP[]/SP[] entry that's still
             * counted -- exactly the kind of mismatch that made fill_write_buffer()'s A_BH
             * scan walk off the end of the array in the original bug this mirrors. */
#ifdef BLACKHOLES
            if(P[i].Type == 5)
              {
                int bhid = P[i].BhID;
                printf(
                    "DEBUG_BHSWAP: Task=%d eliminating P[%d] ID=%llu bhid=%d NumBhs=%d -- source BhP[NumBhs-1=%d]: "
                    "PID=%d TimeBinBh=%d ID-of-owner=%llu\n",
                    ThisTask, i, (unsigned long long)P[i].ID, bhid, NumBhs, NumBhs - 1, BhP[NumBhs - 1].PID,
                    BhP[NumBhs - 1].TimeBinBh, (unsigned long long)P[BhP[NumBhs - 1].PID].ID);
                myflush(stdout);
                BhP[bhid] = BhP[NumBhs - 1];
                P[BhP[bhid].PID].BhID = bhid;
                NumBhs--;
                printf("DEBUG_BHSWAP: Task=%d -> BhP[%d] now TimeBinBh=%d PID=%d ID-of-owner=%llu, new NumBhs=%d\n", ThisTask,
                       bhid, BhP[bhid].TimeBinBh, BhP[bhid].PID, (unsigned long long)P[BhP[bhid].PID].ID, NumBhs);
                myflush(stdout);
              }
#endif
#ifdef STARS
            if(P[i].Type == 4)
              {
                int sid = P[i].SID;
                SP[sid] = SP[NumStars - 1];
                P[SP[sid].PID].SID = sid;
                NumStars--;
              }
#endif

            /* non-gas particle being eliminated (a Type 4 star with zero mass, or any
             * non-gas type matching Mass==0 && ID==0): swap the true last particle into
             * this slot before NumPart shrinks, mirroring the gas branch above. Without
             * this, NumPart simply shrinks without removing this slot's stale data,
             * silently dropping whatever particle genuinely occupied the last slot
             * instead (and leaving the eliminated particle's own stale, still-flagged-
             * for-elimination data sitting within the new, smaller NumPart bound) */
            P[i]   = P[NumPart - 1];
            Key[i] = Key[NumPart - 1];

#ifdef STARS
            if(P[i].Type == 4)
              SPP(i).PID = i;
#endif

#ifdef BLACKHOLES
            if(P[i].Type == 5)
              BPP(i).PID = i;
#endif
          }

        NumPart--;
        i--;
        count_elim++;
      }

  int count[2] = {count_elim, count_gaselim};
  int tot[2] = {0, 0}, nelem = 2;

  MPI_Allreduce(count, tot, nelem, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("DOMAIN: Eliminated %d derefined/swallowed gas cells.\n", tot[1]);
      myflush(stdout);
    }

  All.TotNumPart -= tot[0];
  All.TotNumGas -= tot[1];

#endif /* #ifdef REFINEMENT_MERGE_CELLS */
}