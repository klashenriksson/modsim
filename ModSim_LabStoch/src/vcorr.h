#include "define.h"
#ifndef VCORR_H
#define VCORR_H

typedef struct vcorr
{
    int n_time_slots;
    int n_particles;
    int curr_slot;
    double* circ_buff;
} vcorr_t;

vcorr_t init_vcorr(int n_particles, int time_slots);

void vcorr_step(double* vel, vcorr_t* vcorr);

double vcorr_calc_correlation(vcorr_t* vcorr, int steps);

#endif