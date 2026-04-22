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

#include "fof.h"
#include "fof_seeding.h"

static MyIDType *MinID=NULL;
static int *Head=NULL, *Len=NULL, *Next=NULL, *Tail=NULL, *MinIDTask=NULL;

/*! \brief Main routine to execute the friend of friends group finder.
 *
 *  Does a FOF search to find halos and seed them provided they satisfy the 
 *  seeding criteria.
 *
 *  \return list of halos to be seeded
 */
int fof_seeding_list(MyIDType *halo_ids, int max_ids)
{
  int i, start, lenloc, largestgroup;
  double t0, t1, cputime;

  TIMER_START(CPU_FOF);

  mpi_printf("FOF_SEEDING: Begin to compute FoF group catalogue...  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  /* check */
  for(i = 0; i < NumPart; i++)
    if((P[i].Mass == 0 && P[i].ID == 0))
      terminate("this should not happen");

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

  // timebin_make_list_of_active_particles_up_to_timebin(&TimeBinsGravity, All.HighestActiveTimeBin);
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

  /* Collect groups needing seeding on this task */
  int n_to_seed = 0;
  
  for(int n = 0; n < Ngroups; n++)
    {
      if(Group[n].Mass < All.MinHaloMassForFOFSeeding) continue;
      if(halo_is_seeded(&HaloSeeds, Group[n].MinID)) continue;
      halo_mark_seeded(&HaloSeeds, Group[n].MinID);

      if (n_to_seed >= max_ids) 
        terminate("Too many halos to seed on this task! Increase max_ids or reduce the number of halos being seeded.");
      
      halo_ids[n_to_seed++] = Group[n].MinID;
    }

  myfree_movable(FOF_GList);
  myfree_movable(FOF_PList);

  myfree_movable(Group);   
  myfree_movable(PS);              

  TIMER_STOP(CPU_FOF);
  mpi_printf("FOF_SEEDING: All FOF related work finished...\n");

  return n_to_seed; // Return the number of halos to seed on this task

}

#endif // #ifdef(HALO_SEEDING)
