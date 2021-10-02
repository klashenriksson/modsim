function [I,S] = diffuse_network_iteration_no_infect(i,s,nbor_list, nbor_len_list)
    numnodes = numel(i);
    I = i;
    S = zeros(numnodes, 1);
    for n = 1:numnodes
        num_nbors = nbor_len_list(n);
        curr_neighbors = nbor_list(n,1:num_nbors);
        
        for individual = 1:s(n)
            neighbor_index = randi(num_nbors);
            nbor = curr_neighbors(neighbor_index);
            S(nbor) = S(nbor) + 1;
        end
    end
end