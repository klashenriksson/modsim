#include <stdlib.h>

#include "vcorr.h"

int wrap(int n, int n_max) {
    int m = n % n_max;
    if (m < 0)
    {
        m = n_max + m;
    }

    return m;
}

vcorr_t init_vcorr(int n_particles, int time_slots) {
    vcorr_t vcorr;
    vcorr.circ_buff = (double*)calloc(n_particles * D * time_slots, sizeof(double));
    vcorr.n_time_slots = time_slots;
    vcorr.n_particles = n_particles;
    vcorr.curr_slot = 0;

    return vcorr;
}

void vcorr_step(double* vel, vcorr_t* vcorr) {
    int index = vcorr->curr_slot % vcorr->n_time_slots;

    for (int i = 0; i < vcorr->n_particles; i++)
    {
        for (int d = 0; d < D; d++)
        {
            vcorr->circ_buff[index * D * i + d] = vel[i];
        }
    }

    vcorr->curr_slot += 1;
}

double vcorr_calc_correlation(vcorr_t* vcorr, int n_slots) {
    double accum_corr = 0.f;
    int base_index = wrap(vcorr->curr_slot-1, vcorr->n_time_slots);

    for (int i = 0; i < n_slots; i++) {
        int time_slot_index = wrap(base_index - i, vcorr->n_time_slots);

        for (int n = 0; n < vcorr->n_particles; n++) {
            double* base_vel = vcorr->circ_buff + base_index * D * n;
            double* new_time_vel = vcorr->circ_buff + time_slot_index * D * n;

            for (int d = 0; d < D; d++) {
                accum_corr += base_vel[d]*new_time_vel[d];
            }
        }
    }

    return accum_corr / n_slots / vcorr->n_particles;
}