clear
W = load("nets/pownet_N100_gamma-2.5.mat");
W = W.W;
W(W ~= 0) = 1;
g = graph(W);

% avg population per node
rho = 2;
rho = max(rho,1.5);
%%
pop = zeros(g.numnodes,1);
for i = 1:g.numnodes
    population = rand() * rho + rho/2;
    pop(i) = round(population);
end
g.Nodes.Population = pop;


num_nodes_of_deg = zeros(max(degree(g)),1);
for j = 1:g.numnodes
    deg = degree(g,j); 
    num_nodes_of_deg(deg) = num_nodes_of_deg(deg) + 1;
end

old_pop = g.Nodes.Population;
for i = 1:4000
    new_pop = zeros(g.numnodes,1);
    
    for n = 1:g.numnodes
        curr_neighbors = find(W(n,:) == 1);
        for p = 1:old_pop(n)
            neighbor_index = randi(size(curr_neighbors,2));
            new_pop(curr_neighbors(neighbor_index)) = new_pop(curr_neighbors(neighbor_index)) + 1;
        end
    end
    
    old_pop = new_pop;
end
g.Nodes.Population = new_pop;

%% task 4
pop_per_deg = zeros(max(degree(g)),1);
rho_k = zeros(1,size(pop_per_deg,1));
for j = 1:g.numnodes
   deg = degree(g,j); 
   pop_per_deg(deg) = pop_per_deg(deg) + g.Nodes.Population(j);
end

expected_rho_k = zeros(max(degree(g)), 1);
for j = 1:size(rho_k,2)
   expected_rho_k(j) = (j./mean(degree(g)))*rho;
end

for j = 1:size(rho_k,2)
    if num_nodes_of_deg(j) == 0
        rho_k(j) = 0;
    else
        rho_k(j) = pop_per_deg(j)./num_nodes_of_deg(j);
    end
end

figure;
for i = 1:size(rho_k,2)
   if rho_k(i) ~= 0
       plot(i, rho_k(i), '*');
       hold on
   end
end

fplot(@(x) (x/mean(degree(g)))*rho, [0 10]);
hold on
for i = 1:size(expected_rho_k,1)
   if expected_rho_k(i) ~= 0
       plot(i, expected_rho_k(i), 's');
       hold on
   end
end

%% task 5
close all
pop = zeros(g.numnodes,1);
for i = 1:g.numnodes
    population = rand() * rho + rho/2;
    pop(i) = round(population);
end

initial_infected_fraction = 0.2;

g.Nodes.S = pop;
g.Nodes.I = zeros(size(pop,1),1);

num_infected = 0;
no_susp_node = [];
while num_infected ~= round(initial_infected_fraction*sum(pop))
    node = randi(g.numnodes);
    if find(no_susp_node == node) ~= 0
       continue; 
    end
    
    g.Nodes.I(node) = g.Nodes.I(node) + 1;
    g.Nodes.S(node) = g.Nodes.S(node) - 1;
    
    if g.Nodes.S(node) == 0
        no_susp_node(end+1) = node;
    end
    
    num_infected = num_infected + 1;
end

g.Nodes.Population = pop;

mu = 2;
beta = 1;

final_net = stoch_sim_over_network(g, W, mu, beta, 0, 50, 2000);
figure;

many_infected_nodes = [];
for i = 1:final_net.numnodes
   if final_net.Nodes.Population(i) == 0
       continue;
   end
   
   if final_net.Nodes.I(i)/final_net.Nodes.Population(i) > 0.3
       many_infected_nodes(end+1) = i;
   end
end
p = plot(final_net);
highlight(p, many_infected_nodes);

function final_net = stoch_sim_over_network(net, adj_mat, mu,beta, t_start, t_end, iterations)
    figure;
    dt = (t_start-t_end)./iterations;
    
    t = zeros(1,iterations);
    
    old_net = net;
    t(1) = t_start;
    for iter = 2:iterations
        final_net = stoch_sim_iteration(old_net, mu, beta, dt);
        final_net = diffuse_network_iteration(final_net, adj_mat);
        t(iter) = t(iter-1) + dt;
        
        if mod(iter, 10) == 0
            plot(t(iter), sum(final_net.Nodes.I)/sum(final_net.Nodes.Population),'*');
            hold on
        end
    end
end

function final_net = diffuse_network_iteration(net, adj_mat)
    p = zeros(size(net.Nodes.Population, 1), 1);
    new_i = p;
    new_s = p;
    for n = 1:net.numnodes
        curr_neighbors = find(adj_mat(n,:) == 1);
        
        transfered_i = 0;
        transfered_s = 0;
        for individual = 1:net.Nodes.Population(n)
            neighbor_index = randi(numel(curr_neighbors));
            nbor = curr_neighbors(neighbor_index);
            p(nbor) = p(nbor) + 1;
            
            if transfered_i ~= net.Nodes.I(n)
                new_i(nbor) = new_i(nbor) + 1;
                transfered_i = transfered_i + 1; 
            elseif transfered_s ~= net.Nodes.S(n)
                new_s(nbor) = new_s(nbor) + 1;
                transfered_s = transfered_s + 1; 
            end
        end
    end
    
    final_net = net;
    final_net.Nodes.Population = p;
    final_net.Nodes.I = new_i;
    final_net.Nodes.S = new_s;
end

function final_net = stoch_sim_iteration(net, mu, beta, dt)
    final_net = net;
    for n_id = 1:net.numnodes
        i = net.Nodes.I(n_id);
        s = net.Nodes.S(n_id);
        
        k_i_to_s = mu*dt;
        k_s_to_i = 1 - (1-beta*dt).^i;
        
        for infected = 1:i
           r = rand();
           if r < k_i_to_s
              final_net.Nodes.I(n_id) = final_net.Nodes.I(n_id) - 1;
              final_net.Nodes.S(n_id) = final_net.Nodes.S(n_id) + 1;
           end
        end
        
        for susp = 1:s
            r = rand();
            if r < k_s_to_i
               final_net.Nodes.S(n_id) = final_net.Nodes.S(n_id) - 1;
               final_net.Nodes.I(n_id) = final_net.Nodes.S(n_id) + 1;
            end
        end
        
    end
end