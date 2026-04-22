#ifndef HALO_SEED_REGISTRY_H
#define HALO_SEED_REGISTRY_H

/* FOF seeding public API header */
void fof_seeding_init(int RestartFlag);
void fof_seeding_io(int modus);

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

/* packing/unpacking for restart I/O */
void halo_seed_registry_pack(HaloSeedRegistry *r, MyIDType **buf, int *n);
void halo_seed_registry_unpack(HaloSeedRegistry *r, MyIDType *buf, int n);

#endif