#include "define.h"
#include "sim.h"


#define ICUT6 (1.0 / (CUT * CUT * CUT * CUT * CUT * CUT))



// *** pair_interaction ***
// Lennard-Jones interaction
// Input: distance squared between two particles
// Return value: interaction energy
double pair_interaction(double r2)
{
  double ir2, ir6;

  if (r2 > CUT * CUT)
    return 0.0;

  ir2 = 1.0 / r2;
  ir6 = ir2 * ir2 * ir2;

  return 4 * (ir6 * (ir6 - 1.0) - ICUT6 * (ICUT6 - 1.0));
}


double force_magnitude(double r2)
{
  double r = sqrt(r2);
  double ir2 = 1.0 / r2;
  double ir6 = ir2 * ir2 * ir2;
  double dudr = -(1.0 / r) * ((48. *ir6*ir6) - (24. * ir6));
  return -dudr;
}

// Input:  r2 = distance squared between two particles.
//         dist = vector with distances
// Output: f = vector with the force
void one_force(double *f, double r2, double *dist)
{
  int d;
  double force = force_magnitude(r2);
  double r = sqrt(r2);
  for (d = 0; d < D; d++)
  {
    f[d] = r > 0 ? (dist[d]/r)*force : 0.;
  }
}


// Function 'inrange' helps to implement periodic bc:s.
// Input: r, which has to be in the interval [-max, 2 max)
// Return value: r, in the interval [0,max).
double inrange(double r, double max)
{    
  if (r < 0.0)
    r += max;
  else
    if (r > max)
      r -= max;
  return r;
}


// One dimensional distance (to the closest mirror point)
double distance(double L, double r1, double r2)
{
  double dist = r1 - r2;
  if (fabs(dist) > L / 2) {
    if (dist > 0.0)
      dist -= L;
    else 
      dist += L;
  }
  return dist;
}


// D-dimensional distance squared.
double dist2(double *L, double *p1, double *p2, double *dist)
{
  int d;
  double sumdist2 = 0.0;

  for (d = 0; d < D; d++) {
    dist[d] = distance(L[d], p1[d], p2[d]);
    sumdist2 += dist[d] * dist[d];
  }

  return sumdist2;
}



// Calculate the total force on each particle
void forces_from_pos(Par *par, double *pos, double *force)
{
  int i, j, d;
  double r2;			// Distance squared between a pair of particles
  double f[D], dist[D];		// Force and distance between a pair of particles
  double *ipos, *jpos;		// Pointers to positions of particles i and j
  
  // Initialize the "total force" array to zero
  for (i = 0; i < par->n; i++)
    for (d = 0; d < D; d++)
      force[D * i + d] = 0.0;

  // Loop over all pairs of particles
  for (i = 0; i < par->n - 1; i++) {
    ipos = pos + D * i;
    for (j = i+1; j < par->n; j++) {
      jpos = pos + D * j;
      
      r2 = dist2(par->L, ipos, jpos, dist);
      if (r2 < CUT * CUT) {
        // Each pair of interacting particles should come here
        // Calculate the force due do this interaction
        one_force(f, r2, dist);		// Calculate the vector f on the basis of r2 and dist
        for (d = 0; d < D; d++) {
          force[D * i + d] += f[d];
          force[D * j + d] -= f[d];
        }
      }
    }
  }
}


#ifdef VEL
// This function takes care of the Langevin terms:
// -alpha * velocity + noise
void langevin_forces(Par *par, double *vel, double *force)
{
  int i, d;
  
  for (i = 0; i < par->n; i++)
    for (d = 0; d < D; d++) {
      // Fix this (5). Use dran_sign() which returns a value between -1 and 1.
      // Fix this (5):   force[D * i + d] += 0;

      double r = dran_sign();
      double vel_term = -par->alpha * vel[D * i + d];
      double noise_term = sqrt(3 * par->zeta) * r;
      force[D * i + d] += vel_term + noise_term;
    }
}

void vel_from_force(Par *par, double *vel, double *force)
{
  int i, d;

  for (i = 0; i < par->n; i++) {
    for (d = 0; d < D; d++) {
      vel[D*i + d] += force[D*i + d] * par->deltat;
    }
  }
}



void pos_from_vel(Par *par, double *pos, double *vel)
{
  int i, d;
  
  for (i = 0; i < par->n; i++) {
    for (d = 0; d < D; d++) {
      double newpos = pos[D * i + d] + par->deltat * vel[D * i + d];
      pos[D * i + d] = inrange(newpos, par->L[d]);
    }
  }
}
#endif

  



int step(Par *par, double *atoms, double *force)
{
  double *pos = atoms;
  
  forces_from_pos(par, pos, force);

#ifdef VEL
  double *vel = atoms + D * par->n;
  // If Langevin simulations then add noise and damping to the forces
  if (par->alpha > 0.0)
    langevin_forces(par, vel, force);

  vel_from_force(par, vel, force);
  
  pos_from_vel(par, pos, vel);
#endif

#ifdef BROWN
  // Brownian dynamics goes here
  pos_from_force(par, pos, force);
#endif
  
  return 0;
}
