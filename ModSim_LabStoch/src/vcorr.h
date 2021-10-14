#include "define.h"
#ifndef VCORR_H
#define VCORR_H

typedef struct vcorr
{
    int n_time_slots;
    int n_particles;
    int curr_slot;
    double* circ_buff;
    double time_step_size;
} vcorr_t;

vcorr_t init_vcorr(int n_particles, double step_size, int n_steps);

void vcorr_step(double* vel, vcorr_t* vcorr);

void vcorr_calc_correlation(vcorr_t* vcorr, double* out_corr_buff);

#endif