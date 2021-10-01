function [I,S] = diffuse_network_iteration_no_susp(i,s,nbor_list, nbor_len_list)
    numnodes = numel(i);
    I = zeros(numnodes, 1);
    S = s;
    for n = 1:numnodes
        num_nbors = nbor_len_list(n);
        curr_neighbors = nbor_list(n,1:num_nbors);
        
        for individual = 1:i(n)
            neighbor_index = randi(num_nbors);
            nbor = curr_neighbors(neighbor_index);
            I(nbor) = I(nbor) + 1;
        end
    end
end