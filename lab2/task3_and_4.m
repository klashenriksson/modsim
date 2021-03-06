%clear
W = load("nets/pownet_N1000_gamma-2.5.mat");
W = W.W;
W(W ~= 0) = 1;
g = graph(W);

% avg population per node
rho = 2;
%%
g = init_graph(g, rho, 0);


num_nodes_of_deg = zeros(max(degree(g)),1);
numnodes = g.numnodes;
for j = 1:numnodes
    deg = degree(g,j); 
    num_nodes_of_deg(deg) = num_nodes_of_deg(deg) + 1;
end

g = stoch_sim_over_network(g, 0,0,0,50,10000,"","");


%% task 3
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

rho_deg_pairs = [];
%figure;
for i = 1:size(rho_k,2)
   if rho_k(i) ~= 0
       rho_deg_pairs(end+1,:) = [i,rho_k(i)];
   end
end

plot(rho_deg_pairs(:,1), rho_deg_pairs(:,2), 'c*', "DisplayName", "Experimental \rho_k, \gamma = -2.5");
hold on
fplot(@(x) (x/mean(degree(g)))*rho, [0 35], "c-", "DisplayName", "Expected \rho_k, \gamma = -2.5");
legend("FontSize", 12);

xlabel("Node Degree", "FontSize", 18);
ylabel("\rho", "FontSize", 18);

linfit = polyfit(rho_deg_pairs(:,1), rho_deg_pairs(:,2), 1);
disp("SLOP NUMERIC: " + linfit(1) + " EXPECTED: " + rho./mean(degree(g)) + " REL ERROR: " + (abs(linfit(1)-rho./mean(degree(g)))./(rho./mean(degree(g)))));

%% task 4
close all
pop = zeros(g.numnodes,1);
for i = 1:g.numnodes
    population = rand() * rho + rho/2;
    pop(i) = round(population);
end

initial_infected_fraction = 0.1;

g.Nodes.S = pop;
g.Nodes.I = zeros(size(pop,1),1);

num_infected = 0;
no_susp_node = [];
intially_infected = round(initial_infected_fraction*sum(pop));
while num_infected ~= intially_infected
    node = randi(g.numnodes);
    if find(no_susp_node == node) ~= 0
       continue; 
    end

    infects_left = intially_infected - num_infected;
    infect_count = min(min(g.Nodes.S(node), randi(rho)), infects_left);
    
    g.Nodes.I(node) = g.Nodes.I(node) + infect_count;
    g.Nodes.S(node) = g.Nodes.S(node) - infect_count;
    
    if g.Nodes.S(node) == 0
        no_susp_node(end+1) = node;
    end
    
    num_infected = num_infected + infect_count;
end

g.Nodes.Population = pop;

mu = 2;
beta = 1;

final_net = stoch_sim_over_network(g, mu, beta, 0, 50, 2000);
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