function [I,S] = diffuse_network_iteration(i,s,adj_mat)
    p = i + s;
    numnodes = numel(p);
    I = zeros(numnodes, 1);
    S = zeros(numnodes, 1);
    for n = 1:numnodes
        curr_neighbors = find(adj_mat(n,:) == 1);
        num_neighbors = numel(curr_neighbors);
        
        transfered_i = 0;
        transfered_s = 0;
        for individual = 1:p(n)
            neighbor_index = randi(num_neighbors);
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