function final_net = stoch_sim_over_network(net, mu,beta, t_start, t_end, iterations, diffuse_ans, spread_ans)

    if ~exist('diffus_ans', 'var')
        diffuse_ans = '';
    end
    if ~exist('spread_ans', 'var')
        spread_ans = '';
    end

    dt = (t_end-t_start)./iterations;
    t = zeros(1,iterations);
    t(1) = t_start;
    adj_mat = net.adjacency;

    i = net.Nodes.I;
    s = net.Nodes.S;
    for iter = 2:iterations
        if strcmp(spread_ans, 'Group')
            [i,s] = stoch_sim_iteration_group_spread(i,s,mu,beta,dt);
        else
            [i,s] = stoch_sim_iteration(i,s,mu,beta,dt);
        end

        if strcmp(diffuse_ans, "NoSusp")
            [i,s] = diffuse_network_iteration_no_susp(i,s,adj_mat);
        else
            [i,s] = diffuse_network_iteration(i, s, adj_mat);
        end
        
        t(iter) = t(iter-1) + dt;
        
       % if mod(iter, 10) == 0
       %     plot(t(iter), sum(i)/sum(i+s),'*');
       %     hold on
       % end
    end

    final_net = net;
    final_net.Nodes.Population = i+s;
    final_net.Nodes.I = i;
    final_net.Nodes.S = s;
end