function final_net = siv_stoch_sim_over_network_world(net, mu,beta, gamma, t_start, t_end, iterations, diffuse_ans, spread_ans)

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
    v = net.Nodes.V;
    ivc = net.Nodes.IVC;

    for iter = 2:iterations
        [i,s,v] = siv_stoch_sim_iteration_world(i,s,v,ivc,mu,beta,gamma, dt,spread_ans);
        [i,s,v] = siv_diffuse_network_iteration(i, s, v, neighbor_list, neighbor_length_list, diffuse_ans);
        
        t(iter) = t(iter-1) + dt;
    end

    final_net = net;
    final_net.Nodes.Population = i+s+v;
    final_net.Nodes.I = i;
    final_net.Nodes.S = s;
    final_net.Nodes.V = v;
end