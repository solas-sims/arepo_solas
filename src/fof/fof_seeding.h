#ifndef HALO_SEED_REGISTRY_H
#define HALO_SEED_REGISTRY_H

/* FOF seeding public API header
 * (requires allvars.h to be included first, for MyIDType) */
void fof_seeding_init(int RestartFlag);
void fof_seeding_io(int modus);

/*! \brief One seeding decision made during an on-the-fly FOF pass.
 *
 *  DonorTask/DonorIndex identify the densest gas cell of the halo; they are
 *  only valid during the same synchronisation point as the FOF pass that
 *  produced them (no domain decomposition may occur in between).
 */
typedef struct halo_seed_event
{
  MyIDType HaloMinID; /*!< MinID identifying the FOF group */
  MyIDType DonorID;   /*!< particle ID of the donor gas cell */
  double HaloMass;    /*!< total FOF mass of the group (code units) */
  int DonorTask;      /*!< MPI task owning the donor cell */
  int DonorIndex;     /*!< local index of the donor cell on DonorTask */
} HaloSeedEvent;

/*! Run FOF, apply seeding criteria, and return the global list of seed
 *  events (identical on all tasks). Marks halos in the seed registry. */
int fof_seeding_list(HaloSeedEvent *events, int max_events);

#ifdef BLACKHOLE_SEEDING
/*! Spawn one black hole per seed event (collective call). */
void seed_black_holes_from_events(HaloSeedEvent *events, int n_events);
#endif

typedef struct HaloSeedRegistry{
    MyIDType *ids;
    int n;
    int max;
} HaloSeedRegistry;

extern HaloSeedRegistry HaloSeeds;

/* lifecycle */
void halo_seed_registry_init(HaloSeedRegistry *r, int restart_flag);
void halo_seed_registry_free(HaloSeedRegistry *r);
void halo_seed_registry_grow(HaloSeedRegistry *r);

/* restart I/O hooks */
void halo_seed_registry_io(HaloSeedRegistry *r, int modus);

/* API */
int  halo_is_seeded(HaloSeedRegistry *r, MyIDType id);
void halo_mark_seeded(HaloSeedRegistry *r, MyIDType id);

#endif