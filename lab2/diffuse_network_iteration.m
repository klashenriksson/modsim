function [I,S] = diffuse_network_iteration(i,s,nbor_list, nbor_len_list, opt)
    if ~exist('opt', 'var')
        opt = "";
    end

    if strcmp(opt, 'Lockdown')
        I = i;
        S = s;
        return;
    end

    p = i + s;
    numnodes = numel(p);
    I = zeros(numnodes, 1);
    S = zeros(numnodes, 1);

    is_no_infect = false;
    if strcmp(opt, 'NoInfect')
        I = i;
        is_no_infect = true;
    end

    is_no_susp = false;
    if strcmp(opt, 'NoSusp')
        S = s;
        is_no_susp = true;
    end

    for n = 1:numnodes
        num_nbors = nbor_len_list(n);
        curr_neighbors = nbor_list(n,1:num_nbors);
        
        transfered_i = 0;
        transfered_s = 0;
        for individual = 1:p(n)
            neighbor_index = randi(num_nbors);
            nbor = curr_neighbors(neighbor_index);
            
            if transfered_i ~= i(n) && ~is_no_infect
                I(nbor) = I(nbor) + 1;
                transfered_i = transfered_i + 1;
            elseif transfered_s ~= s(n) && ~is_no_susp
                S(nbor) = S(nbor) + 1;
                transfered_s = transfered_s + 1; 
            end
        end
    end
end