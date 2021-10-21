function [i,s] = stoch_sim_iteration(i,s,mu, beta, dt, opt)
    if ~exist('opt', 'var')
        opt = "";
    end

    use_beta_i = false;
    if strcmp(opt, 'Group')
        use_beta_i = true;
    end

    numnodes = numel(i);
    for n_id = 1:numnodes
        num_i = i(n_id);
        num_s = s(n_id);

        if use_beta_i
            tot_pop = num_i + num_s;
            if tot_pop ~= 0
                beta_i = beta./(tot_pop);
            else
                beta_i = 0;
            end
        else
            beta_i = beta;
        end
        
        k_i_to_s = mu*dt;
        k_s_to_i = 1 - (1-beta_i*dt).^num_i;

        k_i_to_s = min(k_i_to_s, 1);
        k_s_to_i = min(k_s_to_i, 1);
        
        for infected = 1:num_i
           r = rand();
           if r < k_i_to_s
              i(n_id) = i(n_id) - 1;
              s(n_id) = s(n_id) + 1;
           end
        end
        
        for susp = 1:num_s
            r = rand();
            if r < k_s_to_i
               s(n_id) = s(n_id) - 1;
               i(n_id) = i(n_id) + 1;
            end
        end
        
    end
end