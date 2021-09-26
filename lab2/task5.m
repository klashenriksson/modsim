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

%% setup many sims
W = [load("nets/pownet_N100_gamma-2.5.mat"); 
    load("nets/pownet_N1000_gamma-2.5.mat");
    load("nets/pownet_N10000_gamma-2.5.mat");
    ];

adj_mats = {};
for i = 1:numel(W)
    mat = W(i).W;
    mat(mat ~= 0) = 1;
    adj_mats{end+1} = mat;
end

g10_2 = graph(cell2mat(adj_mats(1)));
g10_3 = graph(cell2mat(adj_mats(2)));
g10_4 = graph(cell2mat(adj_mats(3)));

%% many sims

rho = 2;
mu = 2;
beta = 1;
iters = 2000;
t_end = 500;
initial_infect_fract = 0.3;
close all

g10_2 = init_graph(g10_2,rho,initial_infect_fract);
g10_3 = init_graph(g10_3,rho,initial_infect_fract);
g10_4 = init_graph(g10_4,rho,initial_infect_fract);

g10_2 = stoch_sim_over_network(g10_2, mu,beta, 0, t_end, iters);
g10_3 = stoch_sim_over_network(g10_3, mu,beta, 0, t_end, iters/10);
g10_4 = stoch_sim_over_network(g10_4, mu,beta, 0, t_end, iters/100);

plot(sum(g10_2.Nodes.Population), sum(g10_2.Nodes.I)./sum(g10_2.Nodes.Population),'*');
hold on
plot(sum(g10_3.Nodes.Population), sum(g10_3.Nodes.I)./sum(g10_3.Nodes.Population),'*');
hold on
plot(sum(g10_4.Nodes.Population), sum(g10_4.Nodes.I)./sum(g10_4.Nodes.Population),'*');

%% many sims no susp transfer
rho = 2;
mu = 2;
beta = 1;
iters = 2000;
t_end = 500;
initial_infect_fract = 0.3;
close all

g10_2 = init_graph(g10_2,rho,initial_infect_fract);
g10_3 = init_graph(g10_3,rho,initial_infect_fract);
g10_4 = init_graph(g10_4,rho,initial_infect_fract);

g10_2 = stoch_sim_over_network(g10_2, mu,beta, 0, t_end, iters, "NoSusp");
g10_3 = stoch_sim_over_network(g10_3, mu,beta, 0, t_end, iters/10, "NoSusp");
g10_4 = stoch_sim_over_network(g10_4, mu,beta, 0, t_end, iters/100, "NoSusp");

plot(sum(g10_2.Nodes.Population), sum(g10_2.Nodes.I)./sum(g10_2.Nodes.Population),'*');
hold on
plot(sum(g10_3.Nodes.Population), sum(g10_3.Nodes.I)./sum(g10_3.Nodes.Population),'*');
hold on
plot(sum(g10_4.Nodes.Population), sum(g10_4.Nodes.I)./sum(g10_4.Nodes.Population),'*');

%% many sims different rho

rho_interval = [0.5,2];
rho_iters = 5;
rhos = linspace(rho_interval(1), rho_interval(2), rho_iters);

mu = 2;
beta = 1;
initial_infect_fract = 0.1;
t_end = 500;
iters = 2000;

for i = 1:rho_iters
    rho = rhos(i);

    g10_2 = init_graph(g10_2,rho,initial_infect_fract);
    g10_3 = init_graph(g10_3,rho,initial_infect_fract);
    g10_4 = init_graph(g10_4,rho,initial_infect_fract);
    
    g10_2 = stoch_sim_over_network(g10_2, mu,beta, 0, t_end, iters);
    g10_3 = stoch_sim_over_network(g10_3, mu,beta, 0, t_end, iters/10);
    g10_4 = stoch_sim_over_network(g10_4, mu,beta, 0, t_end, iters/100);
    
    plot(sum(g10_2.Nodes.Population), sum(g10_2.Nodes.I)./sum(g10_2.Nodes.Population),'r*');
    hold on
    plot(sum(g10_3.Nodes.Population), sum(g10_3.Nodes.I)./sum(g10_3.Nodes.Population),'g*');
    hold on
    plot(sum(g10_4.Nodes.Population), sum(g10_4.Nodes.I)./sum(g10_4.Nodes.Population),'b*');
end
