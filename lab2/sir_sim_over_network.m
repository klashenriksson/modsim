function final_net = sir_sim_over_network(net, mu,beta, t_start, t_end, iterations, diffuse_ans, spread_ans)

    if ~exist('diffuse_ans', 'var')
        diffuse_ans = '';
    end
    if ~exist('spread_ans', 'var')
        spread_ans = '';
    end

    dt = (t_end-t_start)./iterations;
    t = zeros(1,iterations);
    t(1) = t_start;
    adj_mat = net.adjacency;
    numnodes = net.numnodes;

    neighbor_list = zeros(numnodes, max(degree(net)));
    neighbor_length_list = zeros(numnodes, 1);
    for i = 1:numnodes
        nbors = find(adj_mat(i,:) == 1);
        num_nbors = numel(nbors);

        neighbor_list(i,1:num_nbors) = nbors;
        neighbor_length_list(i) = num_nbors;
    end

    i = net.Nodes.I;
    s = net.Nodes.S;
    r = net.Nodes.R;

    hold on
    plot(t_start, sum(i), 'r*');
    plot(t_start, sum(s), 'c*');
    plot(t_start, sum(r), 'b*');
    hold on

    for iter = 2:iterations
        if strcmp(spread_ans, 'Group')
            [i,s,r] = sir_sim_iteration_group_spread(i,s,r,mu,beta,dt);
        else
            [i,s,r] = sir_sim_iteration(i,s,r,mu,beta,dt);
        end

        if strcmp(diffuse_ans, "NoSusp")
            error("TODO!!!!");
            [i,s,r] = diffuse_network_iteration_no_susp(i,s,neighbor_list, neighbor_length_list);
        elseif strcmp(diffuse_ans, "NoDiffuse")
        else
            [i,s,r] = sir_diffuse_network_iteration(i, s, r, neighbor_list, neighbor_length_list);
        end
        
        t(iter) = t(iter-1) + dt;
        
        plot(t(iter), sum(i), 'r*');
        plot(t(iter), sum(s), 'c*');
        plot(t(iter), sum(r), 'b*');
        hold on
    end

    legend("I", "S", "R")

    final_net = net;
    final_net.Nodes.Population = i+s;
    final_net.Nodes.I = i;
    final_net.Nodes.S = s;
    final_net.Nodes.R = r;
end