/*!
 * \file        src/sidm/sidm_tree.c
 * \brief       Allocation/free for the dedicated SIDM (DM-only) neighbour
 *              tree. Build and walk logic follow in later files -- this is
 *              step 1: get the node/index-array allocation compiling and
 *              correctly sized before tackling tree construction, which is
 *              the more bug-prone piece (see ngb_treebuild_construct's
 *              domain-topleaf/pseudo-particle handling for the scale of
 *              what that involves).
 *
 *              Sizing mirrors ngb_treeallocate() in src/ngbtree/ngbtree.c,
 *              but uses All.MaxPart rather than All.MaxPartSph -- there is
 *              no All.MaxPartDM in this codebase, so DM particles share
 *              the general per-task particle allocation, same convention
 *              the gravity tree itself uses.
 */

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"
#include "sidm_tree.h"

#ifdef SIDM

struct SidmNODE *SidmTree_Nodes;

int SidmTree_MaxPart;
int SidmTree_MaxNodes;
int SidmTree_NumNodes;
int SidmTree_FirstNonTopLevelNode;
int SidmTree_NextFreeNode;

int *SidmTree_Nextnode;
int *SidmTree_Father;
int *SidmTree_Marker;
int *SidmTree_DomainNodeIndex;

static int *SidmTree_Node_Tmp_Sibling;

static int sidm_treebuild_construct(void);
static int sidm_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z);
static void sidm_record_topnode_siblings(int no, int sib);
static void sidm_update_node_recursive(int no, int sib, int father, int *last, int mode);
static void sidm_exchange_topleafdata(void);
static inline unsigned long long sidm_double_to_int(double d);

/*! \brief Allocates memory for the SIDM DM-only tree.
 *
 *  No All.TotNumGas==0-style early-out here: unlike ngbtree, SIDM's
 *  tree existing is gated purely on the SIDM config flag, not on
 *  whether any particular particle count is currently nonzero on this
 *  task (a task with zero local DM particles still needs a valid,
 *  if trivially small, allocation to safely participate in the
 *  collective domain-topleaf exchange during tree construction).
 */
void sidm_treeallocate(void)
{
  if(SidmTree_MaxPart == 0)
    {
      SidmTree_MaxPart  = All.MaxPart;
      SidmTree_MaxNodes = (int)(All.NgbTreeAllocFactor * (All.MaxPart + BASENUMBER)) + NTopnodes;
    }

  if(SidmTree_Nodes)
    terminate("SIDM tree already allocated");

  SidmTree_DomainNodeIndex =
      (int *)mymalloc_movable(&SidmTree_DomainNodeIndex, "SidmTree_DomainNodeIndex", NTopleaves * sizeof(int));

  SidmTree_Nodes = (struct SidmNODE *)mymalloc_movable(&SidmTree_Nodes, "SidmTree_Nodes", (SidmTree_MaxNodes + 1) * sizeof(struct SidmNODE));
  SidmTree_Nodes -= SidmTree_MaxPart;

  SidmTree_Nextnode = (int *)mymalloc_movable(&SidmTree_Nextnode, "SidmTree_Nextnode", (SidmTree_MaxPart + NTopleaves) * sizeof(int));
  SidmTree_Father    = (int *)mymalloc_movable(&SidmTree_Father, "SidmTree_Father", SidmTree_MaxPart * sizeof(int));

  SidmTree_Marker = (int *)mymalloc_movable(&SidmTree_Marker, "SidmTree_Marker", (SidmTree_MaxNodes + SidmTree_MaxPart) * sizeof(int));

  mpi_printf("SIDM_TREE: allocated (MaxPart=%d MaxNodes=%d).\n", SidmTree_MaxPart, SidmTree_MaxNodes);
}

/*! \brief Frees the memory allocated for the SIDM DM-only tree.
 *
 *  Silently does nothing if not currently allocated -- matches
 *  ngb_treefree()'s own behavior exactly. This matters because run()
 *  calls this unconditionally on its very first pass, before anything
 *  has ever been allocated; ngb_treefree() handles that by simply
 *  having nothing to do, not by treating it as an error.
 */
void sidm_treefree(void)
{
  if(SidmTree_Nodes)
    {
      myfree_movable(SidmTree_Marker);
      myfree_movable(SidmTree_Father);
      myfree_movable(SidmTree_Nextnode);

      myfree_movable(SidmTree_Nodes + SidmTree_MaxPart);
      myfree_movable(SidmTree_DomainNodeIndex);

      SidmTree_Marker          = NULL;
      SidmTree_Father          = NULL;
      SidmTree_Nodes           = NULL;
      SidmTree_DomainNodeIndex = NULL;
      SidmTree_Nextnode        = NULL;
      SidmTree_MaxPart         = 0;
      SidmTree_MaxNodes        = 0;
    }
}

/*! \brief Converts double precision coordinate to unsigned long long int,
 *  for Peano/Morton-key based tree insertion. Identical to
 *  ngb_double_to_int / force_double_to_int -- same bit-trick, just not
 *  worth sharing across three near-identical one-line statics.
 */
static inline unsigned long long sidm_double_to_int(double d)
{
  union
  {
    double d;
    unsigned long long ull;
  } u;
  u.d = d;
  return (u.ull & 0xFFFFFFFFFFFFFllu);
}

/*! \brief Driver routine for constructing the SIDM DM-only tree.
 *
 *  Unlike ngb_treebuild(), there is no All.TotNumGas==0-style skip --
 *  this always builds when SIDM is active, over Type==1 particles only.
 *  Called independently of gravity's tree lifecycle (see sidm_tree.h),
 *  so cadence is entirely under SIDM's own control.
 */
int sidm_treebuild(void)
{
  TIMER_START(CPU_NGBTREEBUILD);

  mpi_printf("SIDM_TREE: construction.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  double t0 = second();

  int flag;
  do
    {
      int flag_single = sidm_treebuild_construct();

      MPI_Allreduce(&flag_single, &flag, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      if(flag == -1)
        {
          myfree(SidmTree_Node_Tmp_Sibling + SidmTree_MaxPart);
          sidm_treefree();

          All.NgbTreeAllocFactor *= 1.15;
          mpi_printf("SIDM_TREE: increasing NgbTreeAllocFactor, new value=%g\n", All.NgbTreeAllocFactor);

          sidm_treeallocate();
        }
    }
  while(flag == -1);

  int ntopleaves = DomainNLocalTopleave[ThisTask];
  int *list      = DomainListOfLocalTopleaves + DomainFirstLocTopleave[ThisTask];

  for(int i = 0; i < ntopleaves; i++)
    {
      int last = -1;
      int no   = SidmTree_DomainNodeIndex[list[i]];

      if(no < SidmTree_MaxPart || no >= SidmTree_MaxPart + SidmTree_MaxNodes)
        terminate("SIDM_TREE: i=%d no=%d  task=%d \n", i, no, DomainTask[list[i]]);

      sidm_update_node_recursive(no, SidmTree_Node_Tmp_Sibling[no], no, &last, 0);

      if(no == last)
        SidmTree_Nodes[no].u.d.nextnode = -1;

      SidmTree_Nodes[no].u.d.sibling = last;
    }

  sidm_exchange_topleafdata();

  for(int i = 0; i < NTopleaves; i++)
    {
      if(DomainTask[i] != ThisTask)
        {
          int index                       = SidmTree_DomainNodeIndex[i];
          SidmTree_Nodes[index].u.d.nextnode = SidmTree_MaxPart + SidmTree_MaxNodes + i;
        }
    }

  int last = -1;
  sidm_update_node_recursive(SidmTree_MaxPart, -1, -1, &last, 1);

  if(last >= SidmTree_MaxPart)
    {
      if(last >= SidmTree_MaxPart + SidmTree_MaxNodes)
        SidmTree_Nextnode[last - SidmTree_MaxNodes] = -1;
      else
        SidmTree_Nodes[last].u.d.nextnode = -1;
    }
  else
    SidmTree_Nextnode[last] = -1;

  double numnodes = SidmTree_NumNodes, tot_numnodes;
  MPI_Reduce(&numnodes, &tot_numnodes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  double t1 = second();
  mpi_printf("SIDM_TREE: construction done. took %g sec  <numnodes>=%g  NTopnodes=%d NTopleaves=%d\n", timediff(t0, t1),
             tot_numnodes / NTask, NTopnodes, NTopleaves);

  myfree(SidmTree_Node_Tmp_Sibling + SidmTree_MaxPart);

  TIMER_STOP(CPU_NGBTREEBUILD);

  return SidmTree_NumNodes;
}

/*! \brief Constructs the SIDM DM-only tree over Type==1 particles.
 *
 *  Mirrors ngb_treebuild_construct()'s Peano/Morton-key insertion
 *  algorithm, reusing the same domain top-level structure (TopNodes,
 *  DomainCorner, DomainInverseLen, DomainTask -- generic domain
 *  infrastructure, not gas-specific). Two deliberate deviations:
 *
 *  1. Iterates the FULL local particle range (0..NumPart-1) with an
 *     inline `if(P[i].Type != 1) continue` filter, rather than a
 *     contiguous npart=NumGas range -- DM particles are not
 *     contiguously indexed in P[] the way gas cells are by domain
 *     exchange convention, so there is no compact range to iterate.
 *     Tree leaves use the raw P[]-index directly (matches the gravity
 *     tree's convention, and what sidm_density_evaluate() already
 *     expects), so this needed no separate compaction pass.
 *
 *  2. Single-threaded (matches NUM_THREADS=1 in this build's
 *     generic_comm_helpers2.h usage) -- the original's thread-safe
 *     batched slot allocation (TAKE_NSLOTS_IN_ONE_GO, first_empty_slot
 *     bookkeeping) is simplified to direct SidmTree_NextFreeNode
 *     increments. Multi-threading can be reintroduced later by
 *     reinstating that batching if profiling shows it's needed.
 */
static int sidm_treebuild_construct(void)
{
  SidmTree_NextFreeNode = SidmTree_MaxPart;

  for(int i = 0; i < 8; i++)
    SidmTree_Nodes[SidmTree_NextFreeNode].u.suns[i] = -1;

  SidmTree_NumNodes = 1;
  SidmTree_NextFreeNode++;

  if(sidm_create_empty_nodes(SidmTree_MaxPart, 0, 1, 0, 0, 0) < 0)
    return -1;

  SidmTree_FirstNonTopLevelNode = SidmTree_NextFreeNode;

  SidmTree_Node_Tmp_Sibling = (int *)mymalloc("SidmTree_Node_Tmp_Sibling", (SidmTree_MaxNodes + 1) * sizeof(int));
  SidmTree_Node_Tmp_Sibling -= SidmTree_MaxPart;

  sidm_record_topnode_siblings(SidmTree_MaxPart, -1);

  unsigned long long *sidmTree_IntPos_list =
      (unsigned long long *)mymalloc("sidmTree_IntPos_list", 3 * NumPart * sizeof(unsigned long long));

  int out_of_space = 0;

  for(int i = 0; i < NumPart && out_of_space == 0; i++)
    {
      if(P[i].Type != 1) /* DM only */
        continue;

      if(P[i].Ti_Current != All.Ti_Current)
        drift_particle(i, All.Ti_Current);

      unsigned long long xxb  = sidm_double_to_int(((P[i].Pos[0] - DomainCorner[0]) * DomainInverseLen) + 1.0);
      unsigned long long yyb  = sidm_double_to_int(((P[i].Pos[1] - DomainCorner[1]) * DomainInverseLen) + 1.0);
      unsigned long long zzb  = sidm_double_to_int(((P[i].Pos[2] - DomainCorner[2]) * DomainInverseLen) + 1.0);
      unsigned long long mask = ((unsigned long long)1) << (52 - 1);
      unsigned char shiftx    = (52 - 1);
      unsigned char shifty    = (52 - 2);
      unsigned char shiftz    = (52 - 3);
      unsigned char levels    = 0;

      sidmTree_IntPos_list[3 * i + 0] = xxb;
      sidmTree_IntPos_list[3 * i + 1] = yyb;
      sidmTree_IntPos_list[3 * i + 2] = zzb;

      int no = 0;
      while(TopNodes[no].Daughter >= 0)
        {
          unsigned char subnode = (((unsigned char)((xxb & mask) >> (shiftx--))) | ((unsigned char)((yyb & mask) >> (shifty--))) |
                                   ((unsigned char)((zzb & mask) >> (shiftz--))));

          mask >>= 1;
          levels++;

          no = TopNodes[no].Daughter + TopNodes[no].MortonToPeanoSubnode[subnode];
        }

      no = TopNodes[no].Leaf;

      if(DomainTask[no] != ThisTask)
        terminate("SIDM_TREE: STOP! ID=%lld of type=%d is inserted into task=%d, but should be on task=%d no=%d\n",
                  (long long)P[i].ID, P[i].Type, ThisTask, DomainTask[no], no);

      int th = SidmTree_DomainNodeIndex[no];

      signed long long centermask = (0xFFF0000000000000llu) >> levels;
      (void)centermask; /* tracked for parity with the gravity/ngb tree construction pattern; not otherwise used here */

      int parent            = -1;
      unsigned char subnode = 0;

      while(1)
        {
          if(th >= SidmTree_MaxPart) /* internal node */
            {
              subnode = (((unsigned char)((xxb & mask) >> (shiftx--))) | ((unsigned char)((yyb & mask) >> (shifty--))) |
                         ((unsigned char)((zzb & mask) >> (shiftz--))));

              mask >>= 1;
              levels++;

              if(levels > MAX_TREE_LEVEL)
                {
                  for(int j = 0; j < 8; j++)
                    {
                      if(SidmTree_Nodes[th].u.suns[subnode] < 0)
                        break;

                      subnode++;
                      if(subnode >= 8)
                        subnode = 7;
                    }
                }

              int nn = SidmTree_Nodes[th].u.suns[subnode];

              if(nn >= 0)
                {
                  parent = th;
                  th     = nn;
                }
              else
                {
                  SidmTree_Nodes[th].u.suns[subnode] = i;
                  break;
                }
            }
          else
            {
              int thold = th;

              th = SidmTree_NextFreeNode;
              SidmTree_NextFreeNode++;

              if(SidmTree_NextFreeNode - SidmTree_MaxPart >= SidmTree_MaxNodes)
                {
                  out_of_space = 1;
                  break;
                }

              SidmTree_Nodes[parent].u.suns[subnode] = th;
              struct SidmNODE *nfreep                = &SidmTree_Nodes[th];

              for(int j = 0; j < 8; j++)
                nfreep->u.suns[j] = -1;

              unsigned long long *intppos = &sidmTree_IntPos_list[3 * thold];

              subnode = (((unsigned char)((intppos[0] & mask) >> shiftx)) | ((unsigned char)((intppos[1] & mask) >> shifty)) |
                         ((unsigned char)((intppos[2] & mask) >> shiftz)));

              nfreep->u.suns[subnode] = thold;
            }
        }
    }

  myfree(sidmTree_IntPos_list);

  if((SidmTree_NumNodes = SidmTree_NextFreeNode - SidmTree_MaxPart) >= SidmTree_MaxNodes)
    {
      if(All.NgbTreeAllocFactor > MAX_TREE_ALLOC_FACTOR)
        {
          terminate("SIDM_TREE: task %d: out of space for SIDM tree.\n", ThisTask);
        }
      else
        return -1;
    }

  return 0;
}

/*! \brief Creates empty tree nodes matching the top-level domain grid.
 *
 *  Direct rename of ngb_create_empty_nodes() -- fully generic domain
 *  top-level walk, no gas-specific logic to strip.
 */
static int sidm_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z)
{
  if(TopNodes[topnode].Daughter >= 0)
    {
      for(int i = 0; i < 2; i++)
        for(int j = 0; j < 2; j++)
          for(int k = 0; k < 2; k++)
            {
              if(SidmTree_NumNodes >= SidmTree_MaxNodes)
                {
                  if(All.NgbTreeAllocFactor > MAX_TREE_ALLOC_FACTOR)
                    terminate("SIDM_TREE: task %d: looks like a serious problem (NTopnodes=%d).\n", ThisTask, NTopnodes);
                  return -1;
                }

              int sub = 7 & peano_hilbert_key((x << 1) + i, (y << 1) + j, (z << 1) + k, bits);

              int count = i + 2 * j + 4 * k;

              SidmTree_Nodes[no].u.suns[count] = SidmTree_NextFreeNode;

              for(int n = 0; n < 8; n++)
                SidmTree_Nodes[SidmTree_NextFreeNode].u.suns[n] = -1;

              if(TopNodes[TopNodes[topnode].Daughter + sub].Daughter == -1)
                SidmTree_DomainNodeIndex[TopNodes[TopNodes[topnode].Daughter + sub].Leaf] = SidmTree_NextFreeNode;

              SidmTree_NextFreeNode++;
              SidmTree_NumNodes++;

              if(sidm_create_empty_nodes(SidmTree_NextFreeNode - 1, TopNodes[topnode].Daughter + sub, bits + 1, 2 * x + i, 2 * y + j,
                                         2 * z + k) < 0)
                return -1;
            }
    }

  return 0;
}

/*! \brief Records sibling pointers for the top-level tree.
 *
 *  Direct rename of ngb_record_topnode_siblings() -- fully generic.
 */
static void sidm_record_topnode_siblings(int no, int sib)
{
  if(SidmTree_Nodes[no].u.suns[0] >= 0)
    {
      SidmTree_Node_Tmp_Sibling[no] = -2;

      for(int j = 0; j < 8; j++)
        {
          int p = SidmTree_Nodes[no].u.suns[j];
          int nextsib;

          if(j < 7)
            nextsib = SidmTree_Nodes[no].u.suns[j + 1];
          else
            nextsib = sib;

          sidm_record_topnode_siblings(p, nextsib);
        }
    }
  else
    SidmTree_Node_Tmp_Sibling[no] = sib;
}

/*! \brief Computes node bounding-box moments and links sibling/nextnode/
 *  father pointers, recursively.
 *
 *  Stripped-down relative to ngb_update_node_recursive(): no
 *  vertex_vmin/vmax (mesh-vertex velocity, meaningless for static DM
 *  particles), no TREE_BASED_TIMESTEPS max-soundspeed tracking, no
 *  STAR_RADIATION_ACTIVE volume/kappa accumulation -- none of those
 *  concepts apply to a DM-only tree. Only range_min/range_max survive,
 *  which is exactly what sidm_density_evaluate()'s neighbour-search
 *  overlap test needs.
 */
static void sidm_update_node_recursive(int no, int sib, int father, int *last, int mode)
{
  int j, jj, k, p, pp, nextsib, suns[8];
  MyNgbTreeFloat range_min[3];
  MyNgbTreeFloat range_max[3];

  if(no >= SidmTree_MaxPart && no < SidmTree_MaxPart + SidmTree_MaxNodes) /* internal node */
    {
      if(*last >= 0)
        {
          if(*last >= SidmTree_MaxPart)
            {
              if(*last == no)
                terminate("SIDM_TREE: last == no");

              if(*last >= SidmTree_MaxPart + SidmTree_MaxNodes)
                SidmTree_Nextnode[*last - SidmTree_MaxNodes] = no;
              else
                SidmTree_Nodes[*last].u.d.nextnode = no;
            }
          else
            SidmTree_Nextnode[*last] = no;
        }

      *last = no;

      int not_interal_top_level = 0;

      if(mode == 1)
        {
          if(!(no >= SidmTree_MaxPart && no < SidmTree_FirstNonTopLevelNode))
            terminate("SIDM_TREE: can't be");

          if(SidmTree_Node_Tmp_Sibling[no] != -2)
            not_interal_top_level = 1;
        }

      if(not_interal_top_level)
        {
          p = SidmTree_Nodes[no].u.d.nextnode;

          if(p >= SidmTree_MaxPart + SidmTree_MaxNodes && p < SidmTree_MaxPart + SidmTree_MaxNodes + NTopleaves)
            sidm_update_node_recursive(p, sib, no, last, mode);
          else
            *last = SidmTree_Nodes[no].u.d.sibling;

          if(SidmTree_Node_Tmp_Sibling[no] != sib)
            terminate("SIDM_TREE: SidmTree_Node_Tmp_Sibling[no] != sib");

          SidmTree_Nodes[no].u.d.sibling = sib;
          SidmTree_Nodes[no].father      = father;
        }
      else
        {
          for(j = 0; j < 8; j++)
            suns[j] = SidmTree_Nodes[no].u.suns[j];

          for(k = 0; k < 3; k++)
            {
              range_min[k] = MAX_NGBRANGE_NUMBER;
              range_max[k] = -MAX_NGBRANGE_NUMBER;
            }

          for(j = 0; j < 8; j++)
            {
              if((p = suns[j]) >= 0)
                {
                  for(jj = j + 1; jj < 8; jj++)
                    if((pp = suns[jj]) >= 0)
                      break;

                  if(jj < 8)
                    nextsib = pp;
                  else
                    nextsib = sib;

                  sidm_update_node_recursive(p, nextsib, no, last, mode);

                  if(p >= SidmTree_MaxPart) /* internal node or pseudo particle */
                    {
                      if(p < SidmTree_MaxPart + SidmTree_MaxNodes) /* not a pseudo particle -- those are still zero, fixed up later */
                        {
                          for(k = 0; k < 3; k++)
                            {
                              if(range_min[k] > SidmTree_Nodes[p].u.d.range_min[k])
                                range_min[k] = SidmTree_Nodes[p].u.d.range_min[k];

                              if(range_max[k] < SidmTree_Nodes[p].u.d.range_max[k])
                                range_max[k] = SidmTree_Nodes[p].u.d.range_max[k];
                            }
                        }
                    }
                  else /* a particle */
                    {
                      for(k = 0; k < 3; k++)
                        {
                          if(range_min[k] > P[p].Pos[k])
                            range_min[k] = P[p].Pos[k];

                          if(range_max[k] < P[p].Pos[k])
                            range_max[k] = P[p].Pos[k];
                        }
                    }
                }
            }

          for(k = 0; k < 3; k++)
            {
              SidmTree_Nodes[no].u.d.range_min[k] = range_min[k];
              SidmTree_Nodes[no].u.d.range_max[k] = range_max[k];
            }

          SidmTree_Nodes[no].u.d.sibling = sib;
          SidmTree_Nodes[no].father      = father;

          SidmTree_Nodes[no].Ti_Current = All.Ti_Current;
        }
    }
  else /* single particle or pseudo particle */
    {
      if(*last >= 0)
        {
          if(*last >= SidmTree_MaxPart)
            {
              if(*last >= SidmTree_MaxPart + SidmTree_MaxNodes)
                SidmTree_Nextnode[*last - SidmTree_MaxNodes] = no;
              else
                SidmTree_Nodes[*last].u.d.nextnode = no;
            }
          else
            SidmTree_Nextnode[*last] = no;
        }

      if(no < SidmTree_MaxPart)
        {
          if(father < SidmTree_MaxPart)
            terminate("SIDM_TREE: no=%d father=%d\n", no, father);

          SidmTree_Father[no] = father;
        }

      *last = no;
    }
}

/*! \brief Exchanges top-leaf node moments (bounding boxes) across tasks.
 *
 *  Stripped-down relative to ngb_exchange_topleafdata(): only
 *  range_min/range_max are exchanged, no vertex/soundspeed/RT extras.
 */
static void sidm_exchange_topleafdata(void)
{
  struct DomainNODE
  {
    MyNgbTreeFloat range_min[3];
    MyNgbTreeFloat range_max[3];
  };

  struct DomainNODE *DomainMoment = (struct DomainNODE *)mymalloc("SidmDomainMoment", NTopleaves * sizeof(struct DomainNODE));

  int *recvcounts = (int *)mymalloc("recvcounts", sizeof(int) * NTask);
  int *recvoffset = (int *)mymalloc("recvoffset", sizeof(int) * NTask);
  int *bytecounts = (int *)mymalloc("bytecounts", sizeof(int) * NTask);
  int *byteoffset = (int *)mymalloc("byteoffset", sizeof(int) * NTask);

  for(int task = 0; task < NTask; task++)
    recvcounts[task] = 0;

  for(int n = 0; n < NTopleaves; n++)
    recvcounts[DomainTask[n]]++;

  for(int task = 0; task < NTask; task++)
    bytecounts[task] = recvcounts[task] * sizeof(struct DomainNODE);

  recvoffset[0] = 0, byteoffset[0] = 0;
  for(int task = 1; task < NTask; task++)
    {
      recvoffset[task] = recvoffset[task - 1] + recvcounts[task - 1];
      byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];
    }

  struct DomainNODE *loc_DomainMoment =
      (struct DomainNODE *)mymalloc("loc_SidmDomainMoment", recvcounts[ThisTask] * sizeof(struct DomainNODE));

  int idx = 0;
  for(int n = 0; n < NTopleaves; n++)
    {
      if(DomainTask[n] == ThisTask)
        {
          int no = SidmTree_DomainNodeIndex[n];

          for(int k = 0; k < 3; k++)
            {
              loc_DomainMoment[idx].range_min[k] = SidmTree_Nodes[no].u.d.range_min[k];
              loc_DomainMoment[idx].range_max[k] = SidmTree_Nodes[no].u.d.range_max[k];
            }

          idx++;
        }
    }

  MPI_Allgatherv(loc_DomainMoment, bytecounts[ThisTask], MPI_BYTE, DomainMoment, bytecounts, byteoffset, MPI_BYTE, MPI_COMM_WORLD);

  for(int task = 0; task < NTask; task++)
    recvcounts[task] = 0;

  for(int n = 0; n < NTopleaves; n++)
    {
      int task = DomainTask[n];
      if(task != ThisTask)
        {
          int no  = SidmTree_DomainNodeIndex[n];
          int idx2 = recvoffset[task] + recvcounts[task]++;

          for(int k = 0; k < 3; k++)
            {
              SidmTree_Nodes[no].u.d.range_min[k] = DomainMoment[idx2].range_min[k];
              SidmTree_Nodes[no].u.d.range_max[k] = DomainMoment[idx2].range_max[k];
            }

          SidmTree_Nodes[no].Ti_Current = All.Ti_Current;
        }
    }

  myfree(loc_DomainMoment);
  myfree(byteoffset);
  myfree(bytecounts);
  myfree(recvoffset);
  myfree(recvcounts);
  myfree(DomainMoment);
}

/*! \brief Prepares export of a SIDM-tree pseudo-particle node to its
 *  owning task.
 *
 *  Direct mirror of ngb_treefind_export_node_threads(), using
 *  SidmTree_-prefixed globals. No EXTENDED_GHOST_SEARCH image_flag --
 *  not used anywhere in this build.
 */
int sidm_treefind_export_node_threads(int no, int target, int thread_id)
{
  int task = DomainTask[no - (SidmTree_MaxPart + SidmTree_MaxNodes)];

  if(Thread[thread_id].Exportflag[task] != target)
    {
      Thread[thread_id].Exportflag[task]     = target;
      int nexp                               = Thread[thread_id].Nexport++;
      Thread[thread_id].PartList[nexp].Task  = task;
      Thread[thread_id].PartList[nexp].Index = target;
      Thread[thread_id].ExportSpace -= Thread[thread_id].ItemSize;
    }

  int nexp                      = Thread[thread_id].NexportNodes++;
  nexp                          = -1 - nexp;
  struct datanodelist *nodelist = (struct datanodelist *)(((char *)Thread[thread_id].PartList) + Thread[thread_id].InitialSpace);
  nodelist[nexp].Task           = task;
  nodelist[nexp].Index          = target;
  nodelist[nexp].Node           = SidmTree_DomainNodeIndex[no - (SidmTree_MaxPart + SidmTree_MaxNodes)];
  Thread[thread_id].ExportSpace -= sizeof(struct datanodelist) + sizeof(int);
  return 0;
}

#endif /* #ifdef SIDM */
