clear
W = load("nets/pownet_N100_gamma-2.5.mat");
W = W.W;
W(W ~= 0) = 1;

g = init_graph(graph(W), 2, 0.1);

mu = 2;
beta = 1;

final_net = stoch_sim_over_network(g, mu, beta, 0, 5, 2000);

many_infected_nodes = [];
for i = 1:final_net.numnodes
   if final_net.Nodes.Population(i) == 0
       continue;
   end
   
   if final_net.Nodes.I(i)/final_net.Nodes.Population(i) > 0.3
       many_infected_nodes(end+1) = i;
   end
end
figure;
p = plot(final_net);
highlight(p, many_infected_nodes);

%%

function g = init_graph(g,rho, infected_fraction)
    pop = zeros(g.numnodes,1);
    total_population = rho*g.numnodes;
    curr_pop = 0;
    while curr_pop < total_population
        node = randi(g.numnodes);
    
        %pop_left = total_population-curr_pop;
        pop(node) = pop(node) + 1;
        curr_pop = curr_pop + 1;
    end
    g.Nodes.Population = pop;
        
    g.Nodes.S = pop;
    g.Nodes.I = zeros(size(pop,1),1);
    
    num_infected = 0;
    no_susp_node = [];
    intially_infected = round(infected_fraction*sum(pop));
    
    nonzero_node_ids = [];
    for i = 1:g.numnodes
        if g.Nodes.Population(i) > 0
            nonzero_node_ids(end+1) = i;
        end
    end

    while num_infected ~= intially_infected
        node = nonzero_node_ids(randi(numel(nonzero_node_ids)));
        if find(no_susp_node == node) ~= 0
           continue; 
        end
    
        infects_left = intially_infected - num_infected;
        infect_count = 1;%min(min(g.Nodes.S(node), randi(rho)), infects_left);
        
        g.Nodes.I(node) = g.Nodes.I(node) + infect_count;
        g.Nodes.S(node) = g.Nodes.S(node) - infect_count;
        
        if g.Nodes.S(node) == 0
            no_susp_node(end+1) = node;
        end
        
        num_infected = num_infected + infect_count;
    end
end

function final_net = stoch_sim_over_network(net, mu,beta, t_start, t_end, iterations)
    dt = (t_end-t_start)./iterations;
    t = zeros(1,iterations);
    t(1) = t_start;
    adj_mat = net.adjacency;

    i = net.Nodes.I;
    s = net.Nodes.S;
    for iter = 2:iterations
        [i,s] = stoch_sim_iteration(i,s,mu,beta,dt);
        [i,s] = diffuse_network_iteration(i, s, adj_mat);
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

function [i,s] = stoch_sim_iteration(i,s,mu, beta, dt)
    numnodes = numel(i);
    for n_id = 1:numnodes
        num_i = i(n_id);
        num_s = s(n_id);
        
        k_i_to_s = mu*dt;
        k_s_to_i = 1 - (1-beta*dt).^num_i;
        
        for infected = 1:num_i
           r = rand();
           if r < k_i_to_s
              i(n_id) = i(n_id) - 1;
              s(n_id) = s(n_id) + 1;
           end
        end
        
        for susp = 1:num_s
            r = rand();
            if r < k_s_to_i
               s(n_id) = s(n_id) - 1;
               i(n_id) = i(n_id) + 1;
            end
        end
        
    end
end