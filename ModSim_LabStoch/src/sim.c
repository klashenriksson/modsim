#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "define.h"
#include "sim.h"

#ifdef VCORR
#include "vcorr.h"
#endif


double *pos, *vel, *force = NULL;



// Construct file names based on the parameters of the run, e.g.
// 0064_r0.500_T1.000_alpha0.10_dt010
// Note that sprintf returns the number of number of characters written to the string
// When the format specifier "%4d" gives "  64" the format "%4.4d" instead
// gives "0064". (We don't want spaces in the file names.)
char *get_filename(Par *par)
{
  static char fname[FNAMESIZE];		// without "static" fname would be just a temporary variable
  char *f = fname;

  f += sprintf(f, "%4.4d", par->n);
  f += sprintf(f, "_r%5.3f", par->rho);
  f += sprintf(f, "_T%5.3f", par->t);
#ifdef MC		// no alpha and deltat in an MC simulation
#else
  f += sprintf(f, "_alpha%4.2f", par->alpha);
  f += sprintf(f, "_dt%3.3d", (int) rint(1000.0 * par->deltat));
#endif
  
  if (f - fname >= FNAMESIZE) {		// f - fname is the size of the string
    fprintf(stderr, "Error: too long file name: %s...\n", fname);
    exit(EXIT_FAILURE);
  }

  return fname;
}


double standard_error(double x, double x2, int nblock)
{
  if (nblock < 1)
  {
    return 0;
  }

  double numerator = fabs(x2-x*x);
  double denom = (double)nblock - 1;
  return sqrt(numerator/denom);
}


void print_standard_error(char *str, double x, double x2, int nblock)
{
  if (nblock == 0)
  {
    printf("%s%g +/- %g\n", str, x, 0.);
  }

  double s = standard_error(x / nblock, x2 / nblock, nblock);
  printf("%s%g +/- %g\n", str, x / nblock, s);
  return;
}


void size_from_rho(Par *par)
{
  int d;

  par->vol = par->n / par->rho;
  for (d = 0; d < D; d++)
    par->L[d] = exp(log(par->vol) / D);
}

// Concatenate two strings and return a pointer to this new string
char *add_strings(char *str1, char *str2)
{
  char *str = malloc(strlen(str1) + strlen(str2) + 1);
  strcpy(str, str1);
  strcat(str, str2);
  return str;
}


// Measure potential energy, kinetic energy, and pressure

void measure(Par *par, double *atoms, Measured *val)
{
  int i, j, d;
  double dist[D], f[D], r2;
  double *pos = atoms;
  double epot = 0.0, virial = 0.0;
#ifdef VEL
  double *vel = atoms + D * par->n;
  double ekin = 0.0;
#endif


  // Potential energy
  for (i = 0; i < par->n; i++) {
    for (j = i + 1; j < par->n; j++) {
      r2 = dist2(par->L, pos + D * i, pos + D * j, dist);
      if (r2 < CUT * CUT) {
        epot += pair_interaction(r2);
	      virial += sqrt(r2) * force_magnitude(r2);
      }
    }
  }

  // Potential energy
  val->epot = epot / par->n;

#ifdef VEL
  // Kinetic energy. We here use m=1.
  for (i = 0; i < par->n; i++)
    for (d = 0; d < D; d++)
      ekin += vel[D * i + d] * vel[D * i + d] / 2.0;

  val->ekin = ekin / par->n;
#endif

  // pressure = (N * T + virial / Dimensionality) / volume
  val->pressure = (par->n * par->t + virial / D) / par->vol;
  return;
}

// Write total energy to stream
void energy_print(FILE *stream, int isamp, Measured *val)
{
  if (stream)
#ifdef VEL
    fprintf(stream, "%d %g\n", isamp, val->epot + val->ekin);
#else
    fprintf(stream, "%d %g\n", isamp, val->epot);
#endif
}
 

void run_simulation(Par *par, double *atoms)
{
  int i, itherm, nstep, istep, isamp, iblock;
  char *filename;
  Measured *val = NULL, *vsum, *v1sum, *v2sum;
  FILE *estream = NULL;
  double epot, ekin, pressure;
  
#ifdef MC
  double naccept = 0;
  printf("Gas with %d particles at T = %g, rho = %g, L = %5.3f,  b = %g\n", par->n, par->t,
	   par->rho, par->L[0], par->b);
#else
    if (par->alpha > 0.0) {
      printf("\n--- Langevin dynamics of a Lennard-Jones gas ---\n\n");
      par->zeta = 2 * par->alpha * par->t / par->deltat;
    }
    else
      printf("\n--- Molecular dynamics of a Lennard-Jones gas ---\n\n");

    printf("Gas with %d particles at rho = %g, T = %g, alpha = %g, deltat = %g, L = %5.3f\n",
	   par->n, par->rho, par->t, par->alpha, par->deltat, par->L[0]);
#endif

    
  // Initialize
  init_ran(par->seed);
  if (!force)
    force = malloc(D * par->n * sizeof(double));

  pos = atoms;
#ifdef VEL
  vel = atoms + D * par->n;
#endif
#ifdef VCORR
  int max_corr_count = 50;
  vcorr_t vcorr = init_vcorr(par->n, max_corr_count);
  int curr_vcorr_count = 0;
  int delta_samps_vcorr_msr = 0.1/par->deltat;
  double t = 0.f;
#endif

  // To store measured values
  if (!val) {
    int nval = sizeof(Measured) / sizeof(double);
    val  = calloc(nval, sizeof(double));	// Values from a single measurement
    vsum  = calloc(nval, sizeof(double));	// Accumulate values from a block with nsamp measurements
    v1sum = calloc(nval, sizeof(double));	// Accumulate block averages
    v2sum = calloc(nval, sizeof(double));	// Accumulate block averages squared (for error estimates)
  }



  // Open file for writing energy results
  filename = get_filename(par);		  // Get file name from parameters
  estream = fopen(add_strings("efile/", filename), "w");

  #ifdef VCORR
  FILE* vcorr_output_file = fopen(add_strings("log/", add_strings("vcorr_", add_strings(filename, ".txt"))), "w");
  #endif

  measure(par, atoms, val);
  printf("Potential energy = %g\n", val->epot);
#ifdef VEL

  printf("Kinetic energy   = %g\n", val->ekin);
#endif
  // Equilibrate
#ifndef MC
  nstep = rint(1.0 / par->deltat);
#endif
  step(par, atoms, force);	// Fix this (3)
  //printf("x, y, vx, vy = %g  %g  %g  %g\n", pos[0], pos[1], vel[0], vel[1]);
  if (par->ntherm) {
    printf("Equilibrate: %d...", par->ntherm);
    fflush(stdout);
    for (itherm = 0; itherm < par->ntherm; itherm++) {
      #ifdef MC
            mc_sweep(par, pos);
      #else    
        for (istep = 0; istep < nstep; istep++) {
          step(par, atoms, force);
        }
      #endif
    }
    printf("done\n");
  }

  
  printf("\nSimulate %d blocks x %d samples each: ", par->nblock, par->nsamp);
  fflush(stdout);


  // Production part
  // This is the main loop
  for (iblock = 0; iblock < par->nblock; iblock++) {
    vsum->epot   = 0.0;
    vsum->pressure = 0.0;
#ifdef VEL
    vsum->ekin   = 0.0;
    vsum->etot   = 0.0;
#endif
    
    for (isamp = 0; isamp < par->nsamp; isamp++) {

#ifdef MC
      naccept += mc_sweep(par, pos);
#else
      // This advances one unit of time since nstep = 1/deltat
      for (istep = 0; istep < nstep; istep++) {
          #ifdef VCORR
          if (istep % delta_samps_vcorr_msr == 0 && curr_vcorr_count <= max_corr_count)
          {
            vcorr_step(vel, &vcorr);
            curr_vcorr_count += 1;

            double correlation = vcorr_calc_correlation(&vcorr, curr_vcorr_count);
            fprintf(vcorr_output_file, "%g %g %g\n", t, (double)(isamp*nstep + istep)/delta_samps_vcorr_msr, correlation);
          }

          t += par->deltat;
          #endif

	        step(par, atoms, force);
     }	// End of "for (isamp = ..."
#endif
      
     measure(par, atoms, val);

     if (isamp % 100 == 0)	// Write energy from a single measurement, but not too often
	    energy_print(estream, isamp + par->nsamp * iblock, val);

      vsum->epot += val->epot;
      vsum->pressure += val->pressure;
#ifdef VEL
      vsum->ekin += val->ekin;
      double etot = val->epot + val->ekin;
      vsum->etot  += etot;
#endif

    }

    vsum->epot /= par->nsamp;
    v1sum->epot += vsum->epot;
    v2sum->epot += vsum->epot * vsum->epot;

    vsum->pressure /= par->nsamp;
    v1sum->pressure += vsum->pressure;
    v2sum->pressure += vsum->pressure * vsum->pressure;

#ifdef VEL    
    vsum->ekin /= par->nsamp;
    v1sum->ekin += vsum->ekin;
    v2sum->ekin += vsum->ekin * vsum->ekin;
    vsum->etot  /= par->nsamp;
    v1sum->etot += vsum->etot;
    v2sum->etot += vsum->etot * vsum->etot;
#endif
    
    printf("%d ", iblock + 1);	fflush(stdout);
  }	// End of loop "for (iblock..."
  printf("\n");

  if (estream) fclose(estream);

  #ifdef VCORR
  fclose(vcorr_output_file);
  #endif

  // Write configuration to the named file in "conf/"
  write_conf(par, atoms, filename);

  // Print out results
  print_standard_error("Potential E: ", v1sum->epot, v2sum->epot, par->nblock);
#ifdef VEL
  print_standard_error("Kinetic  E : ", v1sum->ekin, v2sum->ekin, par->nblock);
  print_standard_error("Total  E   : ", v1sum->etot, v2sum->etot, par->nblock);
#endif
  print_standard_error("Pressure   : ", v1sum->pressure, v2sum->pressure, par->nblock);
  printf("Pressure for ideal gas: %g\n", par->n * par->t / par->vol);

#ifdef MC
  printf("Acceptance ratio = %g\n", naccept / par->n / par->nsamp / par->nblock);
#endif



  
}

// The arg string is used to assign values to variables in par.
// The strings look like "N=50", "T=1.5" and so on.
// Return value: 1 for success, 0 for failure.
int read_args(Par *par, char *arg)
{
  static double *atoms = NULL;
  char *s;

  s = strchr(arg, '=');
  if (s)		// Let s point at string after '=', if that char is found
    *s++ = '\0';

  if (!strcmp(arg, "read")) {
    atoms = read_conf(par, atoms, s);

    if (atoms)	// Successful?
      return 1;
    else
      return 0;
  }

  
  if (!strcmp(arg, "run")) {
    if (!atoms)
      atoms = init_conf(par);

    run_simulation(par, atoms);
    return 1;
  }

  
  if (!s) {
    fprintf(stderr, "Command '%s' not recognized, expected format: '<name>=<value>'\n", arg);
    return 0;
  }


  if (!strcmp(arg, "N")) {
    par->n = strtol(s, NULL, 10);
    return 1;
  }

  if (!strcmp(arg, "rho")) {
    double oldrho = par->rho;
    par->rho = strtod(s, NULL);
    size_from_rho(par);

    if (pos) {		// Rescale existing configuration
      double fact = exp(log(par->rho/oldrho) / D);
      int id;
      for (id = 0; id < D * par->n; id++)
	pos[id] /= fact;
    }
    return 1;
  }

  if (!strcmp(arg, "T")) {
    par->t = strtod(s, NULL);
    return 1;
  }

  if (!strcmp(arg, "nblock")) {
    par->nblock = strtol(s, NULL, 10);
    return 1;
  }

  if (!strcmp(arg, "ntherm")) {
    par->ntherm = strtol(s, NULL, 10);
    return 1;
  }

  if (!strcmp(arg, "nsamp")) {
    par->nsamp = strtol(s, NULL, 10);
    return 1;
  }

  if (!strcmp(arg, "seed")) {
    par->seed = strtol(s, NULL, 10);
    return 1;
  }

#ifdef MC
  if (!strcmp(arg, "b")) {
    par->b = strtod(s, NULL);
    return 1;
  }
#else
  if (!strcmp(arg, "alpha")) {
    par->alpha = strtod(s, NULL);
    return 1;
  }
  if (!strcmp(arg, "deltat")) {
    par->deltat = strtod(s, NULL);
    return 1;
  }
#endif

  fprintf(stderr, "No such variable name: '%s'.\n", arg);
  return 0;
}


int main(int argc, char *argv[])
{
  int iarg;
  Par par;

  // Number of degrees of freedom per particle pos and vel: 2D, only pos: D.
#ifdef VEL
  par.df = 2 * D;
#else
  par.df = D;
#endif
  par.n = 64;
  par.nblock = 10;
  par.seed = 0;
  par.ntherm = 1000;
  par.nsamp = 1000;
  par.rho = 0.6;
  size_from_rho(&par);	// To set par->L
#ifdef MC
  par.b = 0.4;
#else
  par.alpha = 0.0;
  par.deltat = 0.01;
#endif

  if (argc == 1) {
#ifdef MC
    printf("Usage: %s N=50 rho=0.5 T=0.1 alpha=1.0\n", argv[0]);
    printf("Optional arguments (with defaults) b=%g nblock=%d nsamp=%d ntherm=%d seed=%d\n",
	   par.b, par.nblock, par.nsamp, par.ntherm, par.seed);
#else
    printf("Usage: %s N=50 L=14.0 T=0.1 alpha=1.0\n", argv[0]);
    printf("Optional arguments (with defaults) deltat=%g nblock=%d nsamp=%d ntherm=%d seed=%d\n",
	   par.deltat, par.nblock, par.nsamp, par.ntherm, par.seed);
#endif
    exit(EXIT_SUCCESS);
  }

  // Interpret the commands given in the argument list.
  for (iarg = 1; iarg < argc; iarg++)
    if (!read_args(&par, argv[iarg]))
      exit(EXIT_FAILURE);

  exit(EXIT_SUCCESS);
}

