function [I,S] = diffuse_network_iteration_no_susp(i,s,adj_mat)
    numnodes = numel(i);
    I = zeros(numnodes, 1);
    S = s;
    for n = 1:numnodes
        curr_neighbors = find(adj_mat(n,:) == 1);
        num_neighbors = numel(curr_neighbors);
        
        for individual = 1:i(n)
            neighbor_index = randi(num_neighbors);
            nbor = curr_neighbors(neighbor_index);
            I(nbor) = I(nbor) + 1;
        end
    end
end