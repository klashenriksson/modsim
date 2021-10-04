#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ran.h"

#define FNAMESIZE 100


typedef struct Measured {
  double epot, pressure;
#ifdef VEL
  double ekin, etot;
#endif
} Measured;

// All the values that define the model and the properties
// of the run are stored in a struct variable of type "Par".
typedef struct Par {
  int nblock, nsamp, ntherm, seed, n, df;
  double L[D], vol, rho;
  double t;
#ifdef MC
  double b;
#else
  double zeta;
  double alpha;
  double deltat;
#endif
} Par;


// In sim.c:
extern char *add_strings(char *str1, char *str2);
extern char *get_filename(Par *par);
extern void measure(Par *par, double *atoms, Measured *val);

// In common.c:
extern double pair_interaction(double r2);
extern double force_magnitude(double r2);
extern void one_force(double *f, double r2, double *dist);
extern double inrange(double r, double max);
extern double distance(double L, double r1, double r2);
extern double dist2(double *L, double *p1, double *p2, double *dist);
extern void forces_from_pos(Par *par, double *pos, double *force);
#ifdef VEL
extern void langevin_forces(Par *par, double *vel, double *force);
extern void vel_from_forces(Par *par, double *vel, double *force);
extern void pos_from_vel(Par *par, double *pos, double *vel);
#endif
extern int step(Par *par, double *atoms, double *force);

// In config.c:
extern void init_vel(Par *par, double *vel);
extern void init_pos(Par *par, double *pos);
extern double *init_conf(Par *par);
extern double *read_conf(Par *par, double *atoms, char *fname);
extern int write_conf(Par *par, double *atoms, char *fname);

#ifdef MC
// In (new file) mc.c?
extern int mc_sweep(Par *par, double *pos);
#endif


