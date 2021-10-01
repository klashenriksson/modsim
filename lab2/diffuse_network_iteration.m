function [I,S] = diffuse_network_iteration(i,s,nbor_list, nbor_len_list)
    p = i + s;
    numnodes = numel(p);
    I = zeros(numnodes, 1);
    S = zeros(numnodes, 1);
    for n = 1:numnodes
        num_nbors = nbor_len_list(n);
        curr_neighbors = nbor_list(n,1:num_nbors);
        
        transfered_i = 0;
        transfered_s = 0;
        for individual = 1:p(n)
            neighbor_index = randi(num_nbors);
            nbor = curr_neighbors(neighbor_index);
            
            if transfered_i ~= i(n)
                I(nbor) = I(nbor) + 1;
                transfered_i = transfered_i + 1;
            elseif transfered_s ~= s(n)
                S(nbor) = S(nbor) + 1;
                transfered_s = transfered_s + 1; 
            end
        end
    end
end