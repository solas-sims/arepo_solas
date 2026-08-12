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
 * \file        src/fof/fof_distribute.c
 * \date        05/2018
 * \brief       Communication and reordering routines for FoF.
 * \details     contains functions:
 *                void fof_subfind_exchange(MPI_Comm Communicator)
 *                void fof_reorder_PS(int *Id, int Nstart, int N)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 24.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"
#include "../subfind/subfind.h"
#include "fof.h"

#ifdef FOF

/*! \brief Redistributes the particles according to what is stored in
 *         PS[].TargetTask, and PS[].TargetIndex.
 *
 *  \param[in] Communicator MPI communicator.
 *
 *  \return void
 */
void fof_subfind_exchange(MPI_Comm Communicator)
{
  int nimport, nexport;
  int i, j, n, p, type, ngrp, target;
  int max_load, max_loadsph, load;
  struct particle_data *partBuf;
  struct subfind_data *subBuf;
  struct sph_particle_data *sphBuf;
#ifdef BLACKHOLES
  struct Bh_Particle_Data *bhBuf;
  int max_loadbh;
#endif /* #ifdef BLACKHOLES */
#ifdef STARS
  Star_Particle_Data *starBuf;
  int max_loadstar;
#endif /* #ifdef STARS */

  int CommThisTask, CommNTask;

  MPI_Comm_size(Communicator, &CommNTask);
  MPI_Comm_rank(Communicator, &CommThisTask);

  int old_AllMaxPart    = All.MaxPart;
  int old_AllMaxPartSph = All.MaxPartSph;

  for(type = 0; type < NTYPES; type++)
    {
      size_t ExportSpace = 0.5 * (FreeBytes); /* we will try to grab at most half of the still available memory  */
      size_t PartSpace   = sizeof(struct particle_data) + sizeof(struct subfind_data) + sizeof(struct sph_particle_data);
#ifdef BLACKHOLES
      PartSpace += sizeof(struct Bh_Particle_Data);
#endif /* #ifdef BLACKHOLES */
#ifdef STARS
      PartSpace += sizeof(Star_Particle_Data);
#endif /* #ifdef STARS */
      if(PartSpace > ExportSpace)
        terminate("seems like we have insufficient storage, PartSpace=%lld ExportSpace=%lld", (long long)PartSpace,
                  (long long)ExportSpace);

      int glob_flag = 0;

      do
        {
          /* Every earlier type's pass through this same do-while loop can have bulk-moved any
           * particle locally via the two block memmove() calls below (closing the gap left by
           * that type's own exports, then opening a gap for that type's own imports) -- neither
           * of those knows about BhP[]/SP[] and neither updates the PID/SID back-reference, so
           * by the time this cycle's own export scan runs, BhP[]/SP[].PID for a BH/star that
           * merely got carried along by an earlier type's memmove can be stale. The forward
           * reference, P[].BhID/SID, is never affected by a mere position change (it travels
           * inside the P[] struct itself), so it remains the authoritative source to resync
           * from here, before this cycle's own export logic below reads any PID/SID back
           * references via a swap-with-last. */
#ifdef BLACKHOLES
          for(n = 0; n < NumPart; n++)
            if(P[n].Type == 5)
              BhP[P[n].BhID].PID = n;
#endif /* #ifdef BLACKHOLES */
#ifdef STARS
          for(n = 0; n < NumPart; n++)
            if(P[n].Type == 4)
              SP[P[n].SID].PID = n;
#endif /* #ifdef STARS */

          for(n = 0; n < CommNTask; n++)
            {
              Send_count[n] = 0;
            }

          ptrdiff_t AvailableSpace = ExportSpace; /* this must be a type that can become negative */

          for(n = 0; n < NumPart; n++)
            {
              if(AvailableSpace < 0)
                break;

              if(P[n].Type == type && PS[n].TargetTask != CommThisTask)
                {
                  target = PS[n].TargetTask;

                  if(target < 0 || target >= CommNTask)
                    terminate("n=%d targettask=%d", n, target);

                  AvailableSpace -= PartSpace;

                  Send_count[target]++;
                }
            }

          MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);

          for(j = 0, nimport = 0, nexport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < CommNTask; j++)
            {
              nexport += Send_count[j];
              nimport += Recv_count[j];

              if(j > 0)
                {
                  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
                }
            }

          /* for resize */
          load = (NumPart + nimport - nexport);
          MPI_Allreduce(&load, &max_load, 1, MPI_INT, MPI_MAX, Communicator);

          if(type == 0)
            {
              load = (NumGas + nimport - nexport);
              MPI_Allreduce(&load, &max_loadsph, 1, MPI_INT, MPI_MAX, Communicator);
            }

#ifdef BLACKHOLES
          if(type == 5)
            {
              load = (NumBhs + nimport - nexport);
              MPI_Allreduce(&load, &max_loadbh, 1, MPI_INT, MPI_MAX, Communicator);
            }
#endif /* #ifdef BLACKHOLES */

#ifdef STARS
          if(type == 4)
            {
              load = (NumStars + nimport - nexport);
              MPI_Allreduce(&load, &max_loadstar, 1, MPI_INT, MPI_MAX, Communicator);
            }
#endif /* #ifdef STARS */

          partBuf = (struct particle_data *)mymalloc_movable(&partBuf, "partBuf", nexport * sizeof(struct particle_data));
          subBuf  = (struct subfind_data *)mymalloc_movable(&subBuf, "subBuf", nexport * sizeof(struct subfind_data));
          if(type == 0)
            sphBuf = (struct sph_particle_data *)mymalloc_movable(&sphBuf, "sphBuf", nexport * sizeof(struct sph_particle_data));
#ifdef BLACKHOLES
          if(type == 5)
            bhBuf = (struct Bh_Particle_Data *)mymalloc_movable(&bhBuf, "bhBuf", nexport * sizeof(struct Bh_Particle_Data));
#endif /* #ifdef BLACKHOLES */
#ifdef STARS
          if(type == 4)
            starBuf = (Star_Particle_Data *)mymalloc_movable(&starBuf, "starBuf", nexport * sizeof(Star_Particle_Data));
#endif /* #ifdef STARS */

          for(i = 0; i < CommNTask; i++)
            {
              Send_count[i] = 0;
            }

          AvailableSpace = ExportSpace; /* this must be allowed to become negative */

          int nstay         = 0;
          int delta_numpart = 0;
          int delta_numgas  = 0;

          for(n = 0; n < NumPart; n++)
            {
              if(AvailableSpace < 0)
                break;

              if(P[n].Type == type && PS[n].TargetTask != CommThisTask)
                {
                  target = PS[n].TargetTask;

                  AvailableSpace -= PartSpace;

                  partBuf[Send_offset[target] + Send_count[target]] = P[n];
                  subBuf[Send_offset[target] + Send_count[target]]  = PS[n];

                  if(P[n].Type == 0)
                    {
                      sphBuf[Send_offset[target] + Send_count[target]] = SphP[n];
                      delta_numgas++;
                    }

#ifdef BLACKHOLES
                  if(P[n].Type == 5)
                    {
                      bhBuf[Send_offset[target] + Send_count[target]] = BhP[P[n].BhID];

                      /* reclaim this BH's local BhP[] slot via swap-with-last, mirroring
                       * domain_exchange()'s own handling of a BH leaving this task; P[n].BhID
                       * is read fresh here (not cached earlier), so if an earlier export this
                       * same scan already relocated this particle's BhP[] entry via the swap
                       * below, we pick up the correct, current slot */
                      int bhid = P[n].BhID;
                      BhP[bhid] = BhP[NumBhs - 1];
                      P[BhP[bhid].PID].BhID = bhid;
                      NumBhs--;
                    }
#endif /* #ifdef BLACKHOLES */

#ifdef STARS
                  if(P[n].Type == 4)
                    {
                      starBuf[Send_offset[target] + Send_count[target]] = SP[P[n].SID];

                      /* reclaim this star's local SP[] slot via swap-with-last, same idea as BhP above */
                      int sid = P[n].SID;
                      SP[sid] = SP[NumStars - 1];
                      P[SP[sid].PID].SID = sid;
                      NumStars--;
                    }
#endif /* #ifdef STARS */

                  Send_count[target]++;
                  delta_numpart++;
                }
              else
                {
                  if(nstay != n)
                    {
                      /* now move P[n] to P[nstay] */

                      P[nstay]  = P[n];
                      PS[nstay] = PS[n];

                      if(P[nstay].Type == 0)
                        SphP[nstay] = SphP[n];

                      /* keep BhP[]/SP[].PID in sync with this particle's new position
                       * immediately, not just at the end of the function: the swap-with-last
                       * compaction above (for a *different* particle being exported later in
                       * this very scan) looks up "who currently owns this BhP[]/SP[] slot" via
                       * exactly this PID back-reference, so it must stay correct throughout the
                       * scan, not only once the scan has finished */
#ifdef BLACKHOLES
                      if(P[nstay].Type == 5)
                        BhP[P[nstay].BhID].PID = nstay;
#endif /* #ifdef BLACKHOLES */
#ifdef STARS
                      if(P[nstay].Type == 4)
                        SP[P[nstay].SID].PID = nstay;
#endif /* #ifdef STARS */
                    }

                  nstay++;
                }
            }

          if(delta_numgas > 0)
            if(delta_numpart != delta_numgas)
              terminate("delta_numpart=%d != delta_numgas=%d", delta_numpart, delta_numgas);

          /* now close gap (if present) */
          memmove(P + nstay, P + nstay + delta_numpart, (NumPart - (nstay + delta_numpart)) * sizeof(struct particle_data));
          memmove(PS + nstay, PS + nstay + delta_numpart, (NumPart - (nstay + delta_numpart)) * sizeof(struct subfind_data));

          if(delta_numgas > 0)
            if(NumGas - (nstay + delta_numgas) > 0)
              memmove(SphP + nstay, SphP + nstay + delta_numpart,
                      (NumGas - (nstay + delta_numgas)) * sizeof(struct sph_particle_data));

          NumPart -= delta_numpart;
          NumGas -= delta_numgas;

          /* do resize, but only increase arrays!! (otherwise data in ActiveParticleList etc. gets lost */
          if(max_load > (1.0 - ALLOC_TOLERANCE) * All.MaxPart)
            {
              All.MaxPart = max_load / (1.0 - 2 * ALLOC_TOLERANCE);
              reallocate_memory_maxpart();
              PS = (struct subfind_data *)myrealloc_movable(PS, All.MaxPart * sizeof(struct subfind_data));
            }

          if(type == 0)
            {
              if(max_loadsph > (1.0 - ALLOC_TOLERANCE) * All.MaxPartSph)
                {
                  All.MaxPartSph = max_loadsph / (1.0 - 2 * ALLOC_TOLERANCE);
                  reallocate_memory_maxpartsph();
                }
            }

#ifdef BLACKHOLES
          if(type == 5)
            {
              if(max_loadbh > (1.0 - ALLOC_TOLERANCE) * All.MaxPartBhs)
                {
                  All.MaxPartBhs = max_loadbh / (1.0 - 2 * ALLOC_TOLERANCE);
                  reallocate_memory_maxpartbhs();
                }
            }
#endif /* #ifdef BLACKHOLES */

#ifdef STARS
          if(type == 4)
            {
              if(max_loadstar > (1.0 - ALLOC_TOLERANCE) * All.MaxPartStars)
                {
                  All.MaxPartStars = max_loadstar / (1.0 - 2 * ALLOC_TOLERANCE);
                  reallocate_memory_maxpartstars();
                }
            }
#endif /* #ifdef STARS */

          /* create a gap behind the existing gas particles where we will insert the incoming particles */
          memmove(P + NumGas + nimport, P + NumGas, (NumPart - NumGas) * sizeof(struct particle_data));
          memmove(PS + NumGas + nimport, PS + NumGas, (NumPart - NumGas) * sizeof(struct subfind_data));

          for(i = 0; i < CommNTask; i++)
            Recv_offset[i] += NumGas;

          for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
            {
              target = CommThisTask ^ ngrp;

              if(target < CommNTask)
                {
                  if(Send_count[target] > 0 || Recv_count[target] > 0)
                    {
                      MPI_Sendrecv(partBuf + Send_offset[target], Send_count[target] * sizeof(struct particle_data), MPI_BYTE, target,
                                   TAG_PDATA, P + Recv_offset[target], Recv_count[target] * sizeof(struct particle_data), MPI_BYTE,
                                   target, TAG_PDATA, Communicator, MPI_STATUS_IGNORE);

                      MPI_Sendrecv(subBuf + Send_offset[target], Send_count[target] * sizeof(struct subfind_data), MPI_BYTE, target,
                                   TAG_KEY, PS + Recv_offset[target], Recv_count[target] * sizeof(struct subfind_data), MPI_BYTE,
                                   target, TAG_KEY, Communicator, MPI_STATUS_IGNORE);

                      if(type == 0)
                        MPI_Sendrecv(sphBuf + Send_offset[target], Send_count[target] * sizeof(struct sph_particle_data), MPI_BYTE,
                                     target, TAG_SPHDATA, SphP + Recv_offset[target],
                                     Recv_count[target] * sizeof(struct sph_particle_data), MPI_BYTE, target, TAG_SPHDATA,
                                     Communicator, MPI_STATUS_IGNORE);

#ifdef BLACKHOLES
                      /* Recv_offset[target] was shifted by "+= NumGas" just above (like every other
                       * type's receive offset); subtracting NumGas here recovers the plain prefix-sum
                       * position within this round's BH-only import batch, so incoming BhP data is
                       * appended right after the current NumBhs entries -- mirrors domain_exchange()'s
                       * offset_recv_bhs[0] = NumBhs convention. */
                      if(type == 5)
                        MPI_Sendrecv(bhBuf + Send_offset[target], Send_count[target] * sizeof(struct Bh_Particle_Data), MPI_BYTE, target,
                                     TAG_BHDATA, BhP + NumBhs + (Recv_offset[target] - NumGas),
                                     Recv_count[target] * sizeof(struct Bh_Particle_Data), MPI_BYTE, target, TAG_BHDATA, Communicator,
                                     MPI_STATUS_IGNORE);
#endif /* #ifdef BLACKHOLES */

#ifdef STARS
                      if(type == 4)
                        MPI_Sendrecv(starBuf + Send_offset[target], Send_count[target] * sizeof(Star_Particle_Data), MPI_BYTE, target,
                                     TAG_STARDATA, SP + NumStars + (Recv_offset[target] - NumGas),
                                     Recv_count[target] * sizeof(Star_Particle_Data), MPI_BYTE, target, TAG_STARDATA, Communicator,
                                     MPI_STATUS_IGNORE);
#endif /* #ifdef STARS */
                    }
                }
            }

#ifdef BLACKHOLES
          if(type == 5)
            {
              /* received BH particles always land at P[NumGas .. NumGas+nimport), the same "gap"
               * every non-gas type is inserted into above; BhP-side, they were just received into
               * BhP[NumBhs .. NumBhs+nimport) (old NumBhs). The offset between a received particle's
               * P[]-index and its BhP[]-index is therefore the single constant (NumBhs - NumGas). */
              int bh_offset = NumBhs - NumGas;

              for(p = NumGas; p < NumGas + nimport; p++)
                {
                  P[p].BhID              = bh_offset + p;
                  BhP[bh_offset + p].PID = p;
                }

              NumBhs += nimport;
            }
#endif /* #ifdef BLACKHOLES */

#ifdef STARS
          if(type == 4)
            {
              int star_offset = NumStars - NumGas;

              for(p = NumGas; p < NumGas + nimport; p++)
                {
                  P[p].SID              = star_offset + p;
                  SP[star_offset + p].PID = p;
                }

              NumStars += nimport;
            }
#endif /* #ifdef STARS */

          if(type == 0)
            NumGas += nimport;

          NumPart += nimport;

#ifdef BLACKHOLES
          if(type == 5)
            myfree_movable(bhBuf);
#endif /* #ifdef BLACKHOLES */

#ifdef STARS
          if(type == 4)
            myfree_movable(starBuf);
#endif /* #ifdef STARS */

          if(type == 0)
            myfree_movable(sphBuf);

          myfree_movable(subBuf);
          myfree_movable(partBuf);

          int loc_flag = 0;
          if(AvailableSpace < 0)
            loc_flag = 1;

          MPI_Allreduce(&loc_flag, &glob_flag, 1, MPI_INT, MPI_SUM, Communicator);
          if(glob_flag > 0 && CommThisTask == 0)
            {
              printf(
                  "FOF-DISTRIBUTE: Need to cycle in particle exchange due to memory shortage. type=%d glob_flag=%d ThisTask=%d "
                  "CommThisTask=%d   PartSpace=%lld  ExportSpace=%lld\n",
                  type, glob_flag, ThisTask, CommThisTask, (long long)PartSpace, (long long)ExportSpace);
              fflush(stdout);
            }
        }
      while(glob_flag);
    }

  /* if there was a temporary memory shortage during the exchange, we may had to increase the maximum allocations. Go back to smaller
   * values again if possible */

  load = NumPart;
  MPI_Allreduce(&load, &max_load, 1, MPI_INT, MPI_MAX, Communicator);
  max_load = max_load / (1.0 - 2 * ALLOC_TOLERANCE);
  if(max_load < old_AllMaxPart)
    max_load = old_AllMaxPart;
  if(max_load != All.MaxPart)
    {
      All.MaxPart = max_load;
      reallocate_memory_maxpart();
      PS = (struct subfind_data *)myrealloc_movable(PS, All.MaxPart * sizeof(struct subfind_data));
    }

  load = NumGas;
  MPI_Allreduce(&load, &max_loadsph, 1, MPI_INT, MPI_MAX, Communicator);
  max_loadsph = max_loadsph / (1.0 - 2 * ALLOC_TOLERANCE);
  if(max_loadsph < old_AllMaxPartSph)
    max_loadsph = old_AllMaxPartSph;
  if(max_loadsph != All.MaxPartSph)
    {
      All.MaxPartSph = max_loadsph;
      reallocate_memory_maxpartsph();
    }

  /* finally, let's also address the desired local order according to PS[].TargetIndex */

  struct fof_local_sort_data *mp;
  int *Id;

  if(NumGas)
    {
      mp = (struct fof_local_sort_data *)mymalloc("mp", sizeof(struct fof_local_sort_data) * NumGas);
      Id = (int *)mymalloc("Id", sizeof(int) * NumGas);

      for(i = 0; i < NumGas; i++)
        {
          mp[i].index       = i;
          mp[i].targetindex = PS[i].TargetIndex;
        }

      qsort(mp, NumGas, sizeof(struct fof_local_sort_data), fof_compare_local_sort_data_targetindex);

      for(i = 0; i < NumGas; i++)
        Id[mp[i].index] = i;

      reorder_gas(Id);

      for(i = 0; i < NumGas; i++)
        Id[mp[i].index] = i;

      fof_reorder_PS(Id, 0, NumGas);

      myfree(Id);
      myfree(mp);
    }

  if(NumPart - NumGas > 0)
    {
      mp = (struct fof_local_sort_data *)mymalloc("mp", sizeof(struct fof_local_sort_data) * (NumPart - NumGas));
      mp -= NumGas;

      Id = (int *)mymalloc("Id", sizeof(int) * (NumPart - NumGas));
      Id -= NumGas;

      for(i = NumGas; i < NumPart; i++)
        {
          mp[i].index       = i;
          mp[i].targetindex = PS[i].TargetIndex;
        }

      qsort(mp + NumGas, NumPart - NumGas, sizeof(struct fof_local_sort_data), fof_compare_local_sort_data_targetindex);

      for(i = NumGas; i < NumPart; i++)
        Id[mp[i].index] = i;

      reorder_particles(Id);

      for(i = NumGas; i < NumPart; i++)
        Id[mp[i].index] = i;

      fof_reorder_PS(Id, NumGas, NumPart);

      Id += NumGas;
      myfree(Id);
      mp += NumGas;
      myfree(mp);
    }

  /* P[].BhID/SID (the forward reference to this particle's BhP[]/SP[] entry) stays valid
   * throughout this whole function no matter how many times a particle's own P[]-array position
   * gets shifted -- it travels along with the rest of P[n] on every copy/memmove, and is only
   * ever reassigned to a genuinely different slot by the explicit swap-with-last/import logic
   * above, which keeps it correct by construction. The REVERSE reference, BhP[]/SP[].PID, is a
   * different story: it must equal the CURRENT P[]-index of whichever particle owns that slot,
   * and this function repositions particles locally in several bulk ways that don't know about
   * BhP[]/SP[] at all (the per-type "stay" compaction copies, and the two block memmove() calls
   * that close/open gaps for export/import). reorder_particles(), called just above, only
   * revisits particles that still need to move to reach the FINAL local order (it skips anything
   * already at its target index), so it does not reliably fix up PID for a particle whose
   * position only changed during this function's own internal reshuffling. Rebuilding every PID
   * from the authoritative
   * forward reference in one full local pass, unconditionally, is simpler and safer than trying
   * to patch up every individual repositioning site above. */
#ifdef BLACKHOLES
  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 5)
      BhP[P[i].BhID].PID = i;
#endif /* #ifdef BLACKHOLES */
#ifdef STARS
  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 4)
      SP[P[i].SID].PID = i;
#endif /* #ifdef STARS */
}

/*! \brief Reorders the elements in the PS array according to the indices given
 *         in the ID array.
 *
 *  \param[in, out] ID Array that specifies new index of element in PS array;
 *                  i.e. PS[i] -> PS[ ID[i] ].
 *  \param[in] Nstart Starting index in ID and PS arrays.
 *  \param[in] N Final element +1 in ID and PS arrays.
 *
 *  \return void
 */
void fof_reorder_PS(int *Id, int Nstart, int N)
{
  int i;
  struct subfind_data PSsave, PSsource;
  int idsource, idsave, dest;

  for(i = Nstart; i < N; i++)
    {
      if(Id[i] != i)
        {
          PSsource = PS[i];

          idsource = Id[i];
          dest     = Id[i];

          do
            {
              PSsave = PS[dest];
              idsave = Id[dest];

              PS[dest] = PSsource;
              Id[dest] = idsource;

              if(dest == i)
                break;

              PSsource = PSsave;
              idsource = idsave;

              dest = idsource;
            }
          while(1);
        }
    }
}

#endif /* #ifdef FOF */
