function [i,s] = stoch_sim_iteration_group_spread(i,s,mu, beta, dt)
    numnodes = numel(i);
    for n_id = 1:numnodes
        num_i = i(n_id);
        num_s = s(n_id);

        beta_i = beta./(num_i+num_s);
        
        k_i_to_s = mu*dt;
        k_s_to_i = 1 - (1-beta_i*dt).^num_i;
        
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