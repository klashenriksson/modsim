function adj_mat = gen_network(num_nodes, gamma)

    node_degrees = zeros(num_nodes,1);
    x_min = 2;
    x_max = sqrt(num_nodes);
    g1 = gamma+1;
    a = x_max.^(g1);
    b = x_min.^(g1);
    for i = 1:num_nodes
        r = rand(); 
        x_pl = ((a-b)*r + b).^(1./g1);

        node_degrees(i) = round(x_pl);
    end

    k_c = zeros(sum(node_degrees),1);
    curr_idx = 1;
    for n = 1:num_nodes
        dups = node_degrees(n);
        k_c(curr_idx:curr_idx+dups) = n;
        curr_idx = curr_idx + dups;
    end

%     degs = 1:max(node_degrees);
%     nodes_per_degs = zeros(max(degs), 1);
%     for i = 1:num_nodes
%         nodes_per_degs(node_degrees(i)) = nodes_per_degs(node_degrees(i)) + 1;
%     end
%     plot(degs, nodes_per_degs);

    adj_mat = zeros(num_nodes, num_nodes);
    while ~isempty(k_c)
        n = numel(k_c);
        if numel(unique(k_c)) == 1
            break
        end

        while true
            n1 = randi(n);
            n2 = randi(n);

            node_1 = k_c(n1);
            node_2 = k_c(n2);

            if adj_mat(node_1,node_2) ~= 1 && node_1 ~= node_2
                break;
            end
        end
        adj_mat(node_1,node_2) = 1;
        adj_mat(node_2, node_1) = 1;
        k_c(n1) = [];

        if n1 > n2
            k_c(n2) = [];
        else
            k_c(n2-1) = [];
        end
    end

    adj_mat = sparse(adj_mat);
end