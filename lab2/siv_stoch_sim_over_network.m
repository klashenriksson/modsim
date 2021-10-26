function [final_net, sol] = siv_stoch_sim_over_network(net, mu,beta, gamma, t_start, t_end, iterations, diffuse_ans, spread_ans)

    if ~exist('diffuse_ans', 'var')
        diffuse_ans = '';
    end
    if ~exist('spread_ans', 'var')
        spread_ans = '';
    end

    dt = (t_end-t_start)./iterations;
    t = zeros(iterations, 1);
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
    v = net.Nodes.V;

    i_sol = zeros(iterations, numel(i));
    s_sol = zeros(iterations, numel(s));
    v_sol = zeros(iterations, numel(v));

    i_sol(1, :) = i;
    s_sol(1,:) = s;
    v_sol(1,:) = v;

    for iter = 2:iterations
        [i,s,v] = siv_stoch_sim_iteration(i,s,v,mu,beta,gamma, dt,spread_ans);
        [i,s,v] = siv_diffuse_network_iteration(i, s, v,neighbor_list, neighbor_length_list, diffuse_ans);

        i_sol(iter,:) = i;
        s_sol(iter,:) = s;
        v_sol(iter,:) = v;
        
        t(iter) = t(iter-1) + dt;
    end

    final_net = net;
    final_net.Nodes.Population = i+s+v;
    final_net.Nodes.I = i;
    final_net.Nodes.S = s;
    final_net.Nodes.V = v;

    %todo do per row stuff..
    sol = [t, i_sol, s_sol, v_sol];
end