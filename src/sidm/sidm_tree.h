#ifndef SIDM_TREE_H
#define SIDM_TREE_H

#include "../main/allvars.h"

#ifdef SIDM

/*! \brief Tree node for the dedicated SIDM (DM-only) neighbour tree.
 *
 *  Mirrors struct NgbNODE in src/ngbtree/ngbtree.h, but built over
 *  Type==1 (DM) particles only, at a cadence SIDM controls itself --
 *  independent of both the gravity tree (which under
 *  HIERARCHICAL_GRAVITY only contains the currently-active particle
 *  subset, not a complete DM population -- see Day-N notes on the
 *  FLAG_FULL_TREE non-convergence bug) and the existing ngbtree module
 *  (which is gas-only, sized off All.MaxPartSph, and not safely
 *  extensible to also carry DM without auditing every caller that
 *  assumes its neighbour indices are gas cells).
 *
 *  Deliberately does NOT carry vertex_vmin/vertex_vmax (mesh-vertex
 *  velocity extrapolation -- meaningless for static DM particles, only
 *  relevant to the moving Voronoi mesh), and there are no SIDM
 *  equivalents of ExtNgbNODE (hydro sound-speed timestep bound) or
 *  RtNgbNODE (radiative transfer) -- neither concept applies to DM.
 */
extern struct SidmNODE
{
  union
  {
    int suns[8]; /*!< temporary pointers to daughter nodes, used only during construction */
    struct
    {
      int sibling;
      int nextnode;
      MyNgbTreeFloat range_min[3];
      MyNgbTreeFloat range_max[3];
    } d;
  } u;

  int father;

  integertime Ti_Current; /*!< for drift-particle correctness during the walk */
} * SidmTree_Nodes;

extern int SidmTree_MaxPart;
extern int SidmTree_MaxNodes;
extern int SidmTree_NumNodes;
extern int SidmTree_FirstNonTopLevelNode;
extern int SidmTree_NextFreeNode;

extern int *SidmTree_Nextnode;
extern int *SidmTree_Father;
extern int *SidmTree_Marker;
extern int *SidmTree_DomainNodeIndex;

void sidm_treeallocate(void);
void sidm_treefree(void);
int sidm_treebuild(void);
int sidm_treefind_export_node_threads(int no, int target, int thread_id);

#endif /* #ifdef SIDM */

#endif /* #ifndef SIDM_TREE_H */
