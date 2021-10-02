function [i,s,r] = sir_sim_iteration(i, s, r,mu, beta, dt)
    numnodes = numel(i);
    for n_id = 1:numnodes
        num_i = i(n_id);
        num_s = s(n_id);
        
        k_s_to_i = 1 - (1-beta*dt).^num_i;
        k_i_to_r = mu*dt;

        if k_s_to_i > 1 || k_i_to_r > 1
            disp("Error!! Probabilities must never exceed 1.");
        end
        
        for infected = 1:num_i
           prob = rand();
           if prob < k_i_to_r
              i(n_id) = i(n_id) - 1;
              r(n_id) = r(n_id) + 1;
           end
        end
        
        for susp = 1:num_s
            prob = rand();
            if prob < k_s_to_i
               s(n_id) = s(n_id) - 1;
               i(n_id) = i(n_id) + 1;
            end
        end
        
    end
end