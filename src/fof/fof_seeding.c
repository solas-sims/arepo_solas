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
 * \file        src/fof/fof_seeding.c
 * \date        20/10/2025
 * \brief       Parallel friend of friends (FoF) group finder.
 * \details     contains functions:
 *                int fof_seeding_list(MyIDType *halo_ids, int max_ids)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 24.05.2018 Prepared file for public release -- Rainer Weinberger
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

#include "../domain/domain.h"
#include "../main/allvars.h"
#include "../main/proto.h"
#include "../blackholes/bh_proto.h"
// #include "../subfind/subfind.h"

#ifdef HALO_SEEDING
#ifndef FOF // Ensure that FOF is enabled if HALO_SEEDING is enabled.
#error "HALO_SEEDING requires FOF to be defined"
#endif /* #ifndef FOF */

#ifdef BLACKHOLE_SEEDING
#if !defined(BH_SEED_ON_MASS) && !defined(BH_SEED_ON_ZERO_METALLICITY)
#error "BLACKHOLE_SEEDING requires at least one seeding channel: BH_SEED_ON_MASS and/or BH_SEED_ON_ZERO_METALLICITY"
#endif
#ifdef BH_SEED_ON_ZERO_METALLICITY
#ifndef METALS
#error "BH_SEED_ON_ZERO_METALLICITY requires METALS to be defined (per-cell metal mass is needed to judge a halo pristine)"
#endif /* #ifndef METALS */
#endif /* #ifdef BH_SEED_ON_ZERO_METALLICITY */
#endif /* #ifdef BLACKHOLE_SEEDING */

#include "fof.h"
#include "fof_seeding.h"

static MyIDType *MinID=NULL;
static int *Head=NULL, *Len=NULL, *Next=NULL, *Tail=NULL, *MinIDTask=NULL;

/*! \brief Query used to fetch the total mass of a group from its owning task. */
struct halo_mass_query
{
  MyIDType MinID; /*!< MinID identifying the group */
  MyFloat Mass;   /*!< total FOF mass of the group (filled by the owning task) */
  int Task;       /*!< task owning the group */
};

static int fof_seeding_compare_query_task(const void *a, const void *b)
{
  if(((struct halo_mass_query *)a)->Task < ((struct halo_mass_query *)b)->Task)
    return -1;
  if(((struct halo_mass_query *)a)->Task > ((struct halo_mass_query *)b)->Task)
    return +1;
  return 0;
}

static int fof_seeding_compare_query_minid(const void *a, const void *b)
{
  if(((struct halo_mass_query *)a)->MinID < ((struct halo_mass_query *)b)->MinID)
    return -1;
  if(((struct halo_mass_query *)a)->MinID > ((struct halo_mass_query *)b)->MinID)
    return +1;
  return 0;
}

/*! \brief Tag all gas cells with the total FOF mass of their host halo.
 *
 *  Must be called while FOF_PList / FOF_GList / Group are still allocated
 *  (i.e. from within fof_seeding_list, after fof_finish_group_properties).
 *
 *  Every task queries the owning task of each group it has particles in
 *  (one query per FOF_GList entry), then walks FOF_PList (sorted by MinID)
 *  in tandem with the answered queries to fill SphP[].HostHaloMass. Gas
 *  cells not belonging to any catalogued group are set to 0.
 *
 *  \return void
 */
static void fof_seeding_tag_host_halo_mass(void)
{
  int i, j, ngrp, recvTask, nimport;

  /* reset all gas cells; those in halos are overwritten below */
  for(i = 0; i < NumGas; i++)
    SphP[i].HostHaloMass = 0;

  struct halo_mass_query *query =
      (struct halo_mass_query *)mymalloc("halo_mass_query", (NgroupsExt > 0 ? NgroupsExt : 1) * sizeof(struct halo_mass_query));

  for(i = 0; i < NgroupsExt; i++)
    {
      query[i].MinID = FOF_GList[i].MinID;
      query[i].Task  = FOF_GList[i].MinIDTask;
      query[i].Mass  = 0;
    }

  /* route the queries to the owning tasks */
  mysort(query, NgroupsExt, sizeof(struct halo_mass_query), fof_seeding_compare_query_task);

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;
  for(i = 0; i < NgroupsExt; i++)
    Send_count[query[i].Task]++;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct halo_mass_query *import =
      (struct halo_mass_query *)mymalloc("halo_mass_import", (nimport > 0 ? nimport : 1) * sizeof(struct halo_mass_query));

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++) /* note: ngrp == 0 handles the local (self) queries */
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&query[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct halo_mass_query), MPI_BYTE, recvTask,
                       TAG_DENS_A, &import[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct halo_mass_query), MPI_BYTE,
                       recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  /* answer the queries from the owned group catalogue (Group is sorted by MinID) */
  for(i = 0; i < nimport; i++)
    {
      int lo = 0, hi = Ngroups - 1, found = -1;

      while(lo <= hi)
        {
          int mid = (lo + hi) / 2;

          if(Group[mid].MinID == import[i].MinID)
            {
              found = mid;
              break;
            }
          else if(Group[mid].MinID < import[i].MinID)
            lo = mid + 1;
          else
            hi = mid - 1;
        }

      if(found < 0)
        terminate("FOF_SEEDING: halo mass query for MinID=%llu not found on task %d", (unsigned long long)import[i].MinID, ThisTask);

      import[i].Mass = Group[found].Mass;
    }

  /* return the answers; blocks come back in the order they were sent */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&import[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct halo_mass_query), MPI_BYTE, recvTask,
                       TAG_DENS_B, &query[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct halo_mass_query), MPI_BYTE,
                       recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  myfree(import);

  /* tandem walk: FOF_PList and the answered queries, both sorted by MinID */
  mysort(query, NgroupsExt, sizeof(struct halo_mass_query), fof_seeding_compare_query_minid);

  int start = 0;

  for(i = 0; i < NgroupsExt; i++)
    {
      while(start < NumPart && FOF_PList[start].MinID < query[i].MinID)
        start++;

      while(start < NumPart && FOF_PList[start].MinID == query[i].MinID)
        {
          int p = FOF_PList[start].Pindex;

          if(P[p].Type == 0 && !(P[p].Mass == 0 && P[p].ID == 0))
            SphP[p].HostHaloMass = query[i].Mass;

          start++;
        }
    }

  myfree(query);
}

/*! \brief Main routine to execute the friend of friends group finder.
 *
 *  Does a FOF search to find halos and seed them provided they satisfy the 
 *  seeding criteria.
 *
 *  \return number of seed events (global; events array identical on all tasks)
 */
int fof_seeding_list(HaloSeedEvent *events, int max_events)
{
  int i, start, lenloc, largestgroup;
  double t0, t1, cputime;

  TIMER_START(CPU_FOF);

  mpi_printf("FOF_SEEDING: Begin to compute FoF group catalogue...  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  /* check */
  for(i = 0; i < NumPart; i++)
    if((P[i].Mass == 0 && P[i].ID == 0))
    {
        mpi_printf("TASK %d: After domain_Decomposition ID=%d Mass=%g\n",
                   ThisTask, P[i].ID, P[i].Mass);
        continue;
                  //  terminate("This should not happen")
    }

  /* this structure will hold auxiliary information for each particle, needed only during group finding */
  PS = (struct subfind_data *)mymalloc_movable(&PS, "PS", All.MaxPart * sizeof(struct subfind_data));

  memset(PS, 0, NumPart * sizeof(struct subfind_data));

  /* First, we save the original location of the particles, in order to be able to revert to this layout later on */
  for(i = 0; i < NumPart; i++)
    {
      PS[i].OriginTask  = ThisTask;
      PS[i].OriginIndex = i;
    }

  fof_OldMaxPart    = All.MaxPart;
  fof_OldMaxPartSph = All.MaxPartSph;

  LinkL = fof_get_comoving_linking_length(); 

  mpi_printf("FOF_SEEDING: Comoving linking length: %g    (presently allocated=%g MB)\n", LinkL, AllocatedBytes / (1024.0 * 1024.0));

  MinID     = (MyIDType *)mymalloc("MinID", NumPart * sizeof(MyIDType));
  MinIDTask = (int *)mymalloc("MinIDTask", NumPart * sizeof(int));

  Head = (int *)mymalloc("Head", NumPart * sizeof(int));
  Len  = (int *)mymalloc("Len", NumPart * sizeof(int));
  Next = (int *)mymalloc("Next", NumPart * sizeof(int));
  Tail = (int *)mymalloc("Tail", NumPart * sizeof(int));

  timebin_make_list_of_active_particles_up_to_timebin(&TimeBinsGravity, All.HighestOccupiedTimeBin);

  construct_forcetree(0, 0, 1, All.HighestOccupiedTimeBin); /* build tree for all particles */

  /* initialize link-lists */
  for(i = 0; i < NumPart; i++)
    {
      Head[i] = Tail[i] = i;
      Len[i]            = 1;
      Next[i]           = -1;
      MinID[i]          = P[i].ID;
      MinIDTask[i]      = ThisTask;
    }

  /* call routine to find primary groups */
  cputime = fof_find_groups(MinID, Head, Len, Next, Tail, MinIDTask);
  mpi_printf("FOF_SEEDING: group finding took = %g sec\n", cputime);

#ifdef FOF_SECONDARY_LINK_TARGET_TYPES
  myfree(Father);
  myfree(Nextnode);
  myfree(Tree_Points);

  /* now rebuild the tree with all the types selected as secondary link targets */
  construct_forcetree(0, 0, 2, All.HighestOccupiedTimeBin);
#endif /* #ifdef FOF_SECONDARY_LINK_TARGET_TYPES */

#ifdef HIERARCHICAL_GRAVITY
  timebin_make_list_of_active_particles_up_to_timebin(&TimeBinsGravity, All.HighestActiveTimeBin);
#endif /* #ifdef HIERARCHICAL_GRAVITY */

  /* call routine to attach secondary particles/cells (gas, stars, BHs) to nearest primary (dm) groups;
   * without this, gas cells are never linked into FOF_PList/FOF_GList and HostHaloMass stays 0 everywhere */
  cputime = fof_find_nearest_dmparticle(MinID, Head, Len, Next, Tail, MinIDTask);
  mpi_printf("FOF_SEEDING: attaching gas and star particles to nearest dm particles took = %g sec\n", cputime);

  myfree(Father);
  myfree(Nextnode);
  myfree(Tree_Points);

  force_treefree();

  myfree(Tail);
  myfree(Next);
  myfree(Len);

  t0 = second();

  FOF_PList = (struct fof_particle_list *)mymalloc_movable(&FOF_PList, "FOF_PList", NumPart * sizeof(struct fof_particle_list));

  for(i = 0; i < NumPart; i++)
    {
      FOF_PList[i].MinID     = MinID[Head[i]];
      FOF_PList[i].MinIDTask = MinIDTask[Head[i]];
      FOF_PList[i].Pindex    = i;
    }

  myfree_movable(Head);
  myfree_movable(MinIDTask);
  myfree_movable(MinID);
    
  FOF_GList = (struct fof_group_list *)mymalloc_movable(&FOF_GList, "FOF_GList", sizeof(struct fof_group_list) * NumPart);

  fof_compile_catalogue(); 

  t1 = second();
  mpi_printf("FOF_SEEDING: compiling local group data and catalogue took = %g sec\n", timediff(t0, t1));

  MPI_Allreduce(&Ngroups, &TotNgroups, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  sumup_large_ints(1, &Nids, &TotNids);

  if(TotNgroups > 0)
    {
      int largestloc = 0;

      for(i = 0; i < NgroupsExt; i++)
        if(FOF_GList[i].LocCount + FOF_GList[i].ExtCount > largestloc)
          largestloc = FOF_GList[i].LocCount + FOF_GList[i].ExtCount;
      MPI_Allreduce(&largestloc, &largestgroup, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
  else
    largestgroup = 0;

  mpi_printf("FOF_SEEDING: Total number of FOF groups with at least %d particles: %d\n", FOF_GROUP_MIN_LEN, TotNgroups);
  mpi_printf("FOF_SEEDING: Largest FOF group has %d particles.\n", largestgroup);
  mpi_printf("FOF_SEEDING: Total number of particles in FOF groups: %lld\n", TotNids);

  t0 = second();

  MaxNgroups = 2 * imax(NgroupsExt, TotNgroups / NTask + 1);

  Group = (struct group_properties *) mymalloc_movable(&Group, "Group", sizeof(struct group_properties) * MaxNgroups);

  mpi_printf("FOF_SEEDING: group properties are now allocated.. (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  for(i = 0, start = 0; i < NgroupsExt; i++)
    {
      while(FOF_PList[start].MinID < FOF_GList[i].MinID)
        {
          start++;
          if(start > NumPart)
            terminate("start > NumPart");
        }

      if(FOF_PList[start].MinID != FOF_GList[i].MinID)
        terminate("ID mismatch");

      for(lenloc = 0; start + lenloc < NumPart;)
        if(FOF_PList[start + lenloc].MinID == FOF_GList[i].MinID)
          lenloc++;
        else
          break;

      Group[i].MinID     = FOF_GList[i].MinID;
      Group[i].MinIDTask = FOF_GList[i].MinIDTask;

      fof_compute_group_properties(i, start, lenloc);

      start += lenloc;
    }

  fof_exchange_group_data();

  fof_finish_group_properties();

  t1 = second();
  mpi_printf("FOF_SEEDING: computation of group properties took = %g sec\n", timediff(t0, t1));

  mpi_printf("FOF_SEEDING: Finished computing FoF groups.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  /* tag gas cells with the mass of their host halo (e.g. for halo-mass dependent star formation) */
  t0 = second();
  fof_seeding_tag_host_halo_mass();
  t1 = second();
  mpi_printf("FOF_SEEDING: tagging gas cells with host halo mass took = %g sec\n", timediff(t0, t1));

  /* Collect seed candidates among the groups owned by this task */
  int n_local = 0;

  HaloSeedEvent *local_events =
      (HaloSeedEvent *)mymalloc("local_seed_events", (Ngroups > 0 ? Ngroups : 1) * sizeof(HaloSeedEvent));

  for(int n = 0; n < Ngroups; n++)
    {
      if(halo_is_seeded(&HaloSeeds, Group[n].MinID))
        continue;

      /* never seed a halo that already hosts a black hole, regardless of which
       * channel below would otherwise trigger; Type 5 is a FOF_SECONDARY_LINK_TYPES
       * member so existing BHs are already counted in LenType[5] */
      if(Group[n].LenType[5] > 0)
        continue;

      int seed_this = 0;

#ifdef BH_SEED_ON_MASS
      if(Group[n].Mass >= All.MinHaloMassForFOFSeeding)
        seed_this = 1;
#endif /* #ifdef BH_SEED_ON_MASS */

#ifdef BH_SEED_ON_ZERO_METALLICITY
      /* whole-halo pristine check: even the most enriched gas cell must be
       * below threshold (MaxGasMetallicity == -1 means no gas at all, handled
       * by the MaxGasDens<0 deferral below) */
      if(!seed_this && Group[n].MaxGasMetallicity >= 0 &&
         Group[n].MaxGasMetallicity <= All.ZeroMetallicityThresholdForFOFSeeding)
        seed_this = 1;
#endif /* #ifdef BH_SEED_ON_ZERO_METALLICITY */

      if(!seed_this)
        continue;

      if(Group[n].MaxGasDens < 0)
        {
          /* halo satisfies a seeding channel but contains no gas: do not mark
             it seeded, so it gets another chance at the next FOF pass */
          printf("FOF_SEEDING: Task %d: group MinID=%llu (M=%g) has no gas cell, deferring seeding.\n", ThisTask,
                 (unsigned long long)Group[n].MinID, Group[n].Mass);
          continue;
        }

      local_events[n_local].HaloMinID  = Group[n].MinID;
      local_events[n_local].HaloMass   = Group[n].Mass;
      local_events[n_local].DonorID    = Group[n].MaxGasDensID;
      local_events[n_local].DonorTask  = Group[n].MaxGasDensTask;
      local_events[n_local].DonorIndex = Group[n].MaxGasDensIndex;
      n_local++;
    }

  /* Gather all seed events on all tasks: the donor cell may live on a
   * different task than the group, and applying the same global event list
   * everywhere keeps the seed registries identical across tasks. */
  int *counts  = (int *)mymalloc("seed_counts", NTask * sizeof(int));
  int *bcounts = (int *)mymalloc("seed_bcounts", NTask * sizeof(int));
  int *bdispls = (int *)mymalloc("seed_bdispls", NTask * sizeof(int));

  MPI_Allgather(&n_local, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);

  int n_global = 0;
  for(i = 0; i < NTask; i++)
    {
      bcounts[i] = counts[i] * sizeof(HaloSeedEvent);
      bdispls[i] = n_global * sizeof(HaloSeedEvent);
      n_global += counts[i];
    }

  if(n_global > max_events)
    terminate("FOF_SEEDING: too many halos to seed (%d > max_events=%d)", n_global, max_events);

  MPI_Allgatherv(local_events, n_local * sizeof(HaloSeedEvent), MPI_BYTE, events, bcounts, bdispls, MPI_BYTE, MPI_COMM_WORLD);

  /* every task records every seeded halo -> registries stay in lockstep */
  for(i = 0; i < n_global; i++)
    halo_mark_seeded(&HaloSeeds, events[i].HaloMinID);

  myfree(bdispls);
  myfree(bcounts);
  myfree(counts);
  myfree(local_events);

  myfree_movable(FOF_GList);
  myfree_movable(FOF_PList);

  myfree_movable(Group);   
  myfree_movable(PS);              

  TIMER_STOP(CPU_FOF);
  mpi_printf("FOF_SEEDING: All FOF related work finished, %d seed events identified.\n", n_global);

  return n_global; /* number of seed events, identical on all tasks */
}

#endif // #ifdef(HALO_SEEDING)
