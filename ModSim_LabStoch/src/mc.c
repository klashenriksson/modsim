#include "define.h"
#include "sim.h"
#include <stdio.h>
#include <memory.h>

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

int mc_sweep(Par *par, double *pos)
{
    double newpos[D];
    double dist[D];
    int naccept = 0;

    for (int i = 0; i < par->n; i++)
    {
        for (int d = 0; d < D; d++)
        {
            double p = pos[i * D + d] + dran_sign() * par->b;
            newpos[d] = inrange(p, par->L[d]);
        }

        double delta_pot = 0.f;
        for (int j = 0; j < par->n; j++)
        {
            if (i == j)
            {
                continue;
            }
                
            double old_r2 = dist2(par->L, pos + D * i, pos + D * j, dist);
            double old_pot = 0.;
            if (old_r2 < CUT * CUT) 
            {
                old_pot = pair_interaction(old_r2);
            }

            double new_r2 = dist2(par->L, newpos, pos + D * j, dist);
            double new_pot = 0.;
            if (new_r2 < CUT * CUT)
            {
                new_pot = pair_interaction(new_r2);
            }
            delta_pot += new_pot - old_pot;
        }

        double acc_prob = MIN(exp(-delta_pot/par->t), 1.);
        double r = dran();
        if (r < acc_prob)
        {
            for (int d = 0; d < D; d++)
            {
                pos[D * i + d] = newpos[d];
            }
            naccept++;
        }
    }

    return naccept;
}