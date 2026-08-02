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
 * \file        src/domain_DC_update.c
 * \date        05/2018
 * \brief       Algorithms for voronoi dynamic update
 * \details     contains functions:
 *                void domain_mark_in_trans_table(int i, int task)
 *                void domain_exchange_and_update_DC(void)
 *                int domain_compare_connection_ID(const void *a,
 *                  const void *b)
 *                int domain_compare_local_trans_data_ID(const void *a,
 *                  const void *b)
 *                int domain_compare_recv_trans_data_ID(const void *a,
 *                  const void *b)
 *                int domain_compare_recv_trans_data_oldtask(const void *a,
 *                  const void *b)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 17.05.2018 Prepared file for public release -- Rainer Weinberger
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

struct trans_data *trans_table;
int N_trans;

/*! \brief Data structure for local auxiliary translation table.
 */
static struct local_aux_trans_data
{
  MyIDType ID;
  int new_index;
} * local_trans_data;

/*! \brief Data structure for communicating the translation table.
 */
static struct aux_trans_data
{
  MyIDType ID;
  int old_task;
  int old_index;
  int new_index;
} * send_trans_data, *recv_trans_data;

/*! \brief Data structure for transcribing data.
 */
static struct aux_transscribe_data
{
  int old_index;
  int new_task;
  int new_index;
  int image_flags;
} * send_transscribe_data, *recv_transscribe_data;

/*! \brief Fill translation table.
 *
 *  Mark where cells are moved to and mark in DC accordingly to make sure
 *  they get communicated to the same task.
 *
 *  \param[in] i Index in P and SphP arrays.
 *  \param[in] task Task to which particle i is exported.
 *
 *  \return void
 */
void domain_mark_in_trans_table(int i, int task)
{
  if(Largest_Nvc > 0)
    {
      if(i < NumGas)
        {
          trans_table[i].ID       = P[i].ID;
          trans_table[i].new_task = task;

          int q = SphP[i].first_connection;

          while(q >= 0)
            {
              int qq = DC[q].next;
              if(q == qq)
                {
                  /* Confirmed via production diagnostics: this is a genuinely inactive cell
                   * (TimeBinSynchronized[P[i].TimeBinHydro] == 0) whose claimed connection
                   * range [first_connection, last_connection] should be stable and untouched
                   * -- AREPO's incremental connectivity scheme (voronoi_update_connectivity())
                   * only rebuilds/frees connections for cells in TimeBinsHydro.ActiveParticleList
                   * each step. Despite that, one of this cell's own connection slots partway
                   * through its range has been silently reallocated to a completely different
                   * particle's connection list (confirmed: DC[q].ID belongs to a different
                   * particle than P[i].ID) -- a real gap in that stability guarantee, not yet
                   * root-caused. Walking the (now foreign) slot's .next field routes into
                   * unrelated territory and can loop back on itself here.
                   *
                   * Stop walking this cell's chain at the point of corruption rather than
                   * crashing the whole run. Connections already walked before this point keep
                   * whatever .next/task assignment was already written to them above; slots
                   * from here to the recorded last_connection are left untouched, since they
                   * likely belong to a different particle now and overwriting them would
                   * corrupt that particle's own connection instead. This cell's connections
                   * will be correctly rebuilt from scratch the next time it becomes active
                   * (voronoi_update_connectivity() only rebuilds active cells' connections, so
                   * this defers the fix rather than resolving it -- if this fires often, the
                   * underlying slot-reclaim gap still needs a proper fix). */
                  printf(
                      "WARNING: DOMAIN_DC: Task=%d cell i=%d (ID=%llu, inactive, TimeBinHydro=%d) has a corrupted connection "
                      "chain -- slot q=%d (claimed range [%d,%d]) unexpectedly self-loops and now belongs to a different "
                      "particle (DC[q].ID=%llu) -- stopping this cell's chain walk early rather than crashing\n",
                      ThisTask, i, (unsigned long long)P[i].ID, P[i].TimeBinHydro, q, SphP[i].first_connection,
                      SphP[i].last_connection, (unsigned long long)DC[q].ID);
                  myflush(stdout);
                  break;
                }

              if((P[i].Mass == 0 && P[i].ID == 0) || P[i].Type != 0) /* this cell has been deleted or turned into a star */
                DC[q].next = -1;
              else
                DC[q].next = task; /* we will temporarily use the next variable to store the new task */

              if(q == SphP[i].last_connection)
                break;

              q = qq;
            }
        }
      else if(i < N_trans)
        trans_table[i].new_task = -1; /* this one has been removed by rerrange_particle_sequence() */
    }
}

/*! \brief Communicates connections.
 *
 *  This algorithms communicates Delauny connections and updates them on the
 *  new task.
 *
 *  \return void
 */
void domain_exchange_and_update_DC(void)
{
  double t0 = second();

#if !defined(GRAVITY_NOT_PERIODIC) && !defined(DO_NOT_RANDOMIZE_DOMAINCENTER) && defined(SELFGRAVITY)
  /* remove all image flags, after our box movement stunt they are all incorrect anyway */
  for(int i = 0; i < MaxNvc; i++)
    {
      DC[i].image_flags = 1;
    }
#endif /* #if !defined(GRAVITY_NOT_PERIODIC) && !defined(DO_NOT_RANDOMIZE_DOMAINCENTER) && defined(SELFGRAVITY) */

  /* first, we need to complete the translation table */
  for(int j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(int i = 0; i < N_trans; i++)
    if(trans_table[i].new_task >= 0)
      Send_count[trans_table[i].new_task]++;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  int nimport = 0, nexport = 0;
  Recv_offset[0] = Send_offset[0] = 0;

  for(int j = 0; j < NTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  send_trans_data = mymalloc("send_trans_data", nexport * sizeof(struct aux_trans_data));
  recv_trans_data = mymalloc("recv_trans_data", nimport * sizeof(struct aux_trans_data));

  for(int j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(int i = 0; i < N_trans; i++)
    {
      int task = trans_table[i].new_task;
      if(task >= 0)
        {
          send_trans_data[Send_offset[task] + Send_count[task]].ID        = trans_table[i].ID;
          send_trans_data[Send_offset[task] + Send_count[task]].old_index = i;
          send_trans_data[Send_offset[task] + Send_count[task]].old_task  = ThisTask;
          Send_count[task]++;
        }
    }

  /* exchange the data */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&send_trans_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct aux_trans_data), MPI_BYTE,
                       recvTask, TAG_DENS_B, &recv_trans_data[Recv_offset[recvTask]],
                       Recv_count[recvTask] * sizeof(struct aux_trans_data), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE);
    }

  /* let's now sort the incoming list according to ID */
  mysort(recv_trans_data, nimport, sizeof(struct aux_trans_data), domain_compare_recv_trans_data_ID);

  /* make an auxiliary list for the local particles that we will also sort according to ID */
  local_trans_data = mymalloc("local_trans_data", NumGas * sizeof(struct local_aux_trans_data));
  for(int i = 0; i < NumGas; i++)
    {
      local_trans_data[i].ID        = P[i].ID;
      local_trans_data[i].new_index = i;
    }
  mysort(local_trans_data, NumGas, sizeof(struct local_aux_trans_data), domain_compare_local_trans_data_ID);

  int i, j;
  /* now we go through and put in the new index for matching IDs */
  for(i = 0, j = 0; i < nimport && j < NumGas;)
    {
      if(recv_trans_data[i].ID < local_trans_data[j].ID)
        {
          recv_trans_data[i].new_index = -1; /* this particle has been eliminated */
          i++;
        }
      else if(recv_trans_data[i].ID > local_trans_data[j].ID)
        j++;
      else
        {
          recv_trans_data[i].new_index = local_trans_data[j].new_index;
          i++;
          j++;
        }
    }

  for(; i < nimport; i++)
    recv_trans_data[i].new_index = -1; /* this particle has been eliminated */

  myfree(local_trans_data);

  /* now order the received data by sending task, so that we can return it */
  mysort(recv_trans_data, nimport, sizeof(struct aux_trans_data), domain_compare_recv_trans_data_oldtask);

  /* return the data */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&recv_trans_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct aux_trans_data), MPI_BYTE,
                       recvTask, TAG_DENS_B, &send_trans_data[Send_offset[recvTask]],
                       Send_count[recvTask] * sizeof(struct aux_trans_data), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE);
    }

  /* now let's fill in the new_index entry into the translation table */
  for(int i = 0; i < nexport; i++)
    trans_table[send_trans_data[i].old_index].new_index = send_trans_data[i].new_index;

  myfree(recv_trans_data);
  myfree(send_trans_data);

  /* it's now time to transcribe the task and index fields in the DC list.
   *
   * DC[i].index < 0 marks a connection whose target cell has already been removed (see
   * the "cell has been removed" convention in voronoi_get_connected_particles()) -- its
   * slot can still be marked live (task >= 0) even though there's no real particle left
   * to look up in trans_table[] on the owning task. The second exchange round below
   * (search "count where they should go") already excludes these via `DC[i].index >= 0`
   * and lets them simply vanish from the rebuilt DC[] after the exchange; this first
   * round must exclude them the same way, and consistently across all three loops that
   * touch it (count / populate / copy-back), since they all rely on matching iteration
   * order over MaxNvc to stay in sync with send_transscribe_data[]/recv_transscribe_data[].
   * Skipping this guard here sends a removed cell's connection to its owning task, whose
   * trans_table[] lookup receives an out-of-range index and terminates. */
  for(int j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(int i = 0; i < MaxNvc; i++)
    {
      int task = DC[i].task;
      if(task >= 0 && DC[i].index >= 0)
        {
          if(task >= NTask)
            terminate("i=%d Nvc=%d MaxNvc=%d task=%d\n", i, Nvc, MaxNvc, task);

          Send_count[task]++;
        }
    }

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  nimport = nexport = 0;
  Recv_offset[0] = Send_offset[0] = 0;

  for(int j = 0; j < NTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  send_transscribe_data = mymalloc("send_transscribe_data", nexport * sizeof(struct aux_transscribe_data));
  recv_transscribe_data = mymalloc("recv_transscribe_data", nimport * sizeof(struct aux_transscribe_data));

  for(int j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(int i = 0; i < MaxNvc; i++)
    {
      int task = DC[i].task;
      if(task >= 0 && DC[i].index >= 0)
        {
          send_transscribe_data[Send_offset[task] + Send_count[task]].old_index   = DC[i].index;
          send_transscribe_data[Send_offset[task] + Send_count[task]].image_flags = DC[i].image_flags;
          Send_count[task]++;
        }
    }

  /* exchange the data */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&send_transscribe_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct aux_transscribe_data),
                       MPI_BYTE, recvTask, TAG_DENS_B, &recv_transscribe_data[Recv_offset[recvTask]],
                       Recv_count[recvTask] * sizeof(struct aux_transscribe_data), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE);
    }

  for(int i = 0; i < nimport; i++)
    {
      if(recv_transscribe_data[i].old_index >= N_trans)
        terminate("recv_transscribe_data[i].old_index >= N_trans");

      if(recv_transscribe_data[i].old_index < 0)
        terminate("recv_transscribe_data[i].old_index < 0");

      int old_index = recv_transscribe_data[i].old_index;

      recv_transscribe_data[i].new_task  = trans_table[old_index].new_task;
      recv_transscribe_data[i].new_index = trans_table[old_index].new_index;

#if !defined(GRAVITY_NOT_PERIODIC) && !defined(DO_NOT_RANDOMIZE_DOMAINCENTER) && defined(SELFGRAVITY)
      // Nothing to do here
#else  /* #if !defined(GRAVITY_NOT_PERIODIC) && !defined(DO_NOT_RANDOMIZE_DOMAINCENTER) && defined(SELFGRAVITY) */
      if(recv_transscribe_data[i].new_task >= 0)
        {
          if(trans_table[old_index].wrapped)
            {
              int bitflags = ffs(recv_transscribe_data[i].image_flags) - 1;
              int zbits    = (bitflags / 9);
              int ybits    = (bitflags - zbits * 9) / 3;
              int xbits    = bitflags - zbits * 9 - ybits * 3;

              if(trans_table[old_index].wrapped & 1)
                {
                  if(xbits == 1)
                    xbits = 0;
                  else if(xbits == 0)
                    xbits = 2;
                  else /* xbits == 2 */
                    terminate("b");
                }
              else if(trans_table[old_index].wrapped & 2)
                {
                  if(xbits == 1)
                    {
                      terminate("a");
                    }
                  else if(xbits == 0)
                    xbits = 1;
                  else /* xbits == 2 */
                    xbits = 0;
                }

              if(trans_table[old_index].wrapped & 4)
                {
                  if(ybits == 1)
                    ybits = 0;
                  else if(ybits == 0)
                    ybits = 2;
                  else
                    {
                      terminate("b");
                    }
                }
              else if(trans_table[old_index].wrapped & 8)
                {
                  if(ybits == 1)
                    {
                      terminate("a");
                    }
                  else if(ybits == 0)
                    ybits = 1;
                  else
                    ybits = 0;
                }

              if(trans_table[old_index].wrapped & 16)
                {
                  if(zbits == 1)
                    zbits = 0;
                  else if(zbits == 0)
                    zbits = 2;
                  else
                    {
                      terminate("b");
                    }
                }
              else if(trans_table[old_index].wrapped & 32)
                {
                  if(zbits == 1)
                    {
                      terminate("a");
                    }
                  else if(zbits == 0)
                    zbits = 1;
                  else
                    zbits = 0;
                }

              recv_transscribe_data[i].image_flags = (1 << (zbits * 9 + ybits * 3 + xbits));
            }
        }
#endif /* #if !defined(GRAVITY_NOT_PERIODIC) && !defined(DO_NOT_RANDOMIZE_DOMAINCENTER) && defined(SELFGRAVITY) #else */
    }

  /* now return the data */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&recv_transscribe_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct aux_transscribe_data),
                       MPI_BYTE, recvTask, TAG_DENS_B, &send_transscribe_data[Send_offset[recvTask]],
                       Send_count[recvTask] * sizeof(struct aux_transscribe_data), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE);
    }

  for(int j = 0; j < NTask; j++)
    Send_count[j] = 0;

  /* copy the results over to the DC structure. Must use the exact same condition as the
   * count/populate loops above (task >= 0 && index >= 0) -- this consumes
   * send_transscribe_data[] in the same per-task order it was produced in, so skipping a
   * different set of entries here would desync the two and misattribute results to the
   * wrong connections. DC[i].index itself hasn't been overwritten yet at this point, so
   * it's still safe to test here for the same "cell already removed" condition. */
  for(int i = 0; i < MaxNvc; i++)
    {
      int task = DC[i].task;
      if(task >= 0 && DC[i].index >= 0)
        {
          DC[i].task        = send_transscribe_data[Send_offset[task] + Send_count[task]].new_task;
          DC[i].index       = send_transscribe_data[Send_offset[task] + Send_count[task]].new_index;
          DC[i].image_flags = send_transscribe_data[Send_offset[task] + Send_count[task]].image_flags;
          Send_count[task]++;
        }
    }

  myfree(recv_transscribe_data);
  myfree(send_transscribe_data);

  /* now we can exchange the DC data. The task where each item should go is stored in 'next' at this point */
  for(int j = 0; j < NTask; j++)
    Send_count[j] = 0;

  /* count where they should go */
  for(int i = 0; i < MaxNvc; i++)
    {
      if(DC[i].task >= 0)
        {
          int task = DC[i].next;
          if(task >= 0)
            {
              if(task >= NTask)
                {
                  /* .next was never assigned a valid destination task this round -- this is a
                   * slot past the point where domain_mark_in_trans_table() stopped walking a
                   * corrupted connection chain (see the "corrupted connection chain" WARNING
                   * there), so .next still holds a stale leftover value (e.g. an old
                   * free-list/connection index) instead of a task. Treat it the same as an
                   * already-removed connection (DC[i].index < 0): drop it from the rebuilt
                   * DC[] below instead of sending it to a bogus task. The later loop that
                   * actually populates the send buffer also reads DC[i].next as a task without
                   * re-validating it, so this guard must run before that loop, not just here. */
                  printf(
                      "WARNING: DOMAIN_DC: Task=%d DC[%d] has out-of-range .next=%d (not a valid task in [0,%d)) -- "
                      "dropping this stale connection instead of sending it\n",
                      ThisTask, i, task, NTask);
                  myflush(stdout);
                  DC[i].index = -1;
                }
              else if(DC[i].index >= 0)
                Send_count[task]++;
            }
        }
    }

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  nimport = nexport = 0;
  Recv_offset[0] = Send_offset[0] = 0;

  for(int j = 0; j < NTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  /* make sure that we have enough room to store the new DC list */
  while(nimport > MaxNvc)
    {
      int old_MaxNvc = MaxNvc;
      Mesh.Indi.AllocFacNvc *= ALLOC_INCREASE_FACTOR;
      MaxNvc = Mesh.Indi.AllocFacNvc;
#ifdef VERBOSE
      printf("Task=%d: increase memory allocation, MaxNvc=%d Indi.AllocFacNvc=%g\n", ThisTask, MaxNvc, Mesh.Indi.AllocFacNvc);
#endif /* #ifdef VERBOSE */
      DC = myrealloc_movable(DC, MaxNvc * sizeof(connection));
      for(int n = old_MaxNvc; n < MaxNvc; n++)
        DC[n].task = -1;
    }

  connection *tmpDC = mymalloc("tmpDC", nexport * sizeof(connection));

  for(int j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(int i = 0; i < MaxNvc; i++)
    {
      if(DC[i].task >= 0)
        {
          int task = DC[i].next;

          if(task >= 0 && DC[i].index >= 0)
            tmpDC[Send_offset[task] + Send_count[task]++] = DC[i];
        }
    }

  /* exchange the connection information */

  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&tmpDC[Send_offset[recvTask]], Send_count[recvTask] * sizeof(connection), MPI_BYTE, recvTask, TAG_DENS_B,
                       &DC[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(connection), MPI_BYTE, recvTask, TAG_DENS_B,
                       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  myfree(tmpDC);

  Nvc = nimport;

  /* mark the remaining ones as available */
  for(int i = Nvc; i < MaxNvc - 1; i++)
    {
      DC[i].next = i + 1;
      DC[i].task = -1;
    }
  DC[MaxNvc - 1].next = -1;
  DC[MaxNvc - 1].task = -1;

  if(Nvc < MaxNvc)
    FirstUnusedConnection = Nvc;
  else
    FirstUnusedConnection = -1;

  /* now we need to connect the information to the particles, this we do via the IDs */

  local_trans_data = mymalloc("local_trans_data", NumGas * sizeof(struct local_aux_trans_data));
  for(int i = 0; i < NumGas; i++)
    {
      local_trans_data[i].ID        = P[i].ID;
      local_trans_data[i].new_index = i; /* is here used as rank of the particle */
    }
  mysort(local_trans_data, NumGas, sizeof(struct local_aux_trans_data), domain_compare_local_trans_data_ID);

  mysort(DC, Nvc, sizeof(connection), domain_compare_connection_ID);

  int last = -1;
  for(i = 0, j = 0; i < NumGas && j < Nvc; i++)
    {
      int k = local_trans_data[i].new_index;

      if(P[k].ID < DC[j].ID)
        {
          /* this particle has no connection information (new cell) */
          SphP[k].first_connection = -1;
          SphP[k].last_connection  = -1;
        }
      else if(P[k].ID == DC[j].ID)
        {
          SphP[k].first_connection = j;

          while(j < Nvc)
            {
              SphP[k].last_connection = j;

              if(last >= 0)
                DC[last].next = j;

              last = j;
              j++;
              if(j >= Nvc)
                break;
              if(P[k].ID != DC[j].ID)
                break;
            }
        }
      else
        {
          terminate("strange");
        }
    }

  for(; i < NumGas; i++)
    {
      int k                    = local_trans_data[i].new_index;
      SphP[k].first_connection = -1;
      SphP[k].last_connection  = -1;
    }

  if(last >= 0)
    DC[last].next = -1;

  myfree(local_trans_data);

  double t1 = second();
  mpi_printf("DOMAIN: done with rearranging connection information (took %g sec)\n", timediff(t0, t1));
}

/*! \brief Compare which ID is larger.
 *
 *  For connection data.
 *
 *  \param[in] a Pointer to first object.
 *  \param[in] b Pointer to second object.
 *
 *  \return (-1,0,1) -1 if a->ID is smaller.
 */
int domain_compare_connection_ID(const void *a, const void *b)
{
  if(((connection *)a)->ID < (((connection *)b)->ID))
    return -1;

  if(((connection *)a)->ID > (((connection *)b)->ID))
    return +1;

  return 0;
}

/*! \brief Compare which ID is larger.
 *
 *  For local_aux_trans_data.
 *
 *  \param[in] a Pointer to first object.
 *  \param[in] b Pointer to second object.
 *
 *  \return (-1,0,1) -1 if a->ID is smaller.
 */
int domain_compare_local_trans_data_ID(const void *a, const void *b)
{
  if(((struct local_aux_trans_data *)a)->ID < (((struct local_aux_trans_data *)b)->ID))
    return -1;

  if(((struct local_aux_trans_data *)a)->ID > (((struct local_aux_trans_data *)b)->ID))
    return +1;

  return 0;
}

/*! \brief Compare which ID is larger.
 *
 *  For aux_trans_data.
 *
 *  \param[in] a Pointer to first object.
 *  \param[in] b Pointer to second object.
 *
 *  \return (-1,0,1) -1 if a->ID is smaller.
 */
int domain_compare_recv_trans_data_ID(const void *a, const void *b)
{
  if(((struct aux_trans_data *)a)->ID < (((struct aux_trans_data *)b)->ID))
    return -1;

  if(((struct aux_trans_data *)a)->ID > (((struct aux_trans_data *)b)->ID))
    return +1;

  return 0;
}

/*! \brief Compare which old_task is larger.
 *
 *  For aux_trans_data.
 *
 *  \param[in] a Pointer to first object.
 *  \param[in] b Pointer to second object.
 *
 *  \return (-1,0,1) -1 if a->old_task is smaller.
 */
int domain_compare_recv_trans_data_oldtask(const void *a, const void *b)
{
  if(((struct aux_trans_data *)a)->old_task < (((struct aux_trans_data *)b)->old_task))
    return -1;

  if(((struct aux_trans_data *)a)->old_task > (((struct aux_trans_data *)b)->old_task))
    return +1;

  return 0;
}
