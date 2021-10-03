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
        k_c(curr_idx:curr_idx+dups-1) = n;
        curr_idx = curr_idx + dups;
    end

%     degs = 1:max(node_degrees);
%     nodes_per_degs = zeros(max(degs), 1);
%     for i = 1:num_nodes
%         nodes_per_degs(node_degrees(i)) = nodes_per_degs(node_degrees(i)) + 1;
%     end
%     plot(degs, nodes_per_degs);

    adj_mat = zeros(num_nodes, num_nodes);
    max_fail_count = 20;

    num_elems = numel(k_c);
    while num_elems > 0
        if numel(unique(k_c)) == 1
            break
        end

        if num_elems == 2 && adj_mat(k_c(1), k_c(2)) == 1
            break
        end

        fail_count = 0;
        while fail_count < max_fail_count
            n1 = randi(num_elems);
            n2 = randi(num_elems);

            node_1 = k_c(n1);
            node_2 = k_c(n2);

            if adj_mat(node_1,node_2) ~= 1 && node_1 ~= node_2
                break;
            end

            fail_count = fail_count + 1;
        end

        if fail_count == max_fail_count
            break;
        end

        adj_mat(node_1,node_2) = 1;
        adj_mat(node_2, node_1) = 1;
        k_c(n1) = [];
        num_elems = num_elems -1;
        if num_elems == 0
            break;
        end

        if n1 > n2
            k_c(n2) = [];
        else
            k_c(n2-1) = [];
        end

        num_elems = num_elems - 1;
    end

    adj_mat = sparse(adj_mat);
end