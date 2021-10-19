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

vcorr_t init_vcorr(int n_particles, double step_size, int n_steps) {
    vcorr_t vcorr;
    double time_slots = step_size * n_steps;
    vcorr.circ_buff = (double*)calloc(n_particles * D * time_slots, sizeof(double));
    vcorr.n_time_slots = time_slots;
    vcorr.n_particles = n_particles;
    vcorr.curr_slot = 0;
    vcorr.time_step_size = step_size;

    return vcorr;
}

void vcorr_step(double* vel, vcorr_t* vcorr) {
    int index = vcorr->curr_slot % vcorr->n_time_slots;

    for (int i = 0; i < vcorr->n_particles; i++)
    {
        for (int d = 0; d < D; d++)
        {
            vcorr->circ_buff[index * D * i + d] = vel[D*i+d];
        }
    }

    vcorr->curr_slot += 1;
}

void vcorr_calc_correlation(vcorr_t* vcorr, double* out_corr_buff, int n_delays, double n_delay_step_size) {
    for (int i = 0; i < n_delays; i++)
    {
        int steps_per_deci = n_delay_step_size/vcorr->time_step_size;
        double delay_time = i * 0.1;
        int time_iters = vcorr->n_time_slots - (int)(delay_time/vcorr->time_step_size);

        double accum_corr = 0.f;
        for (int t_idx = 0; t_idx < time_iters; t_idx++) {

            double time = t_idx * vcorr->time_step_size;
            double delayed_time = time + delay_time;
            int delayed_time_index = (int)(delayed_time / vcorr->time_step_size);
            
            for (int n = 0; n < vcorr->n_particles; n++) {
                double* iter_vel = vcorr->circ_buff + t_idx * D * n;
                double* delay_iter_vel = vcorr->circ_buff + delayed_time_index * D * n;

                for (int d = 0; d < D; d++) {
                    accum_corr += iter_vel[d]*delay_iter_vel[d];
                }
            }
        }

        out_corr_buff[i] = accum_corr / vcorr->n_particles / time_iters;
    }
}