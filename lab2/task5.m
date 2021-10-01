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

%% many sims (i.e task 4)

rho_interval = [0.5,2];
rho_iters = 6;
rhos = linspace(rho_interval(1), rho_interval(2), rho_iters);

mu = 2;
beta = 1;
iters = 2000;
t_end = 100;
initial_infect_fract = 0.1;
close all

V10_2 = zeros(rho_iters, 2);
V10_3 = zeros(rho_iters, 2);
V10_4 = zeros(rho_iters, 2);

for i = 1:rho_iters
    rho = rhos(i);

    g10_2 = init_graph(g10_2,rho,initial_infect_fract);
    g10_3 = init_graph(g10_3,rho,initial_infect_fract);
    g10_4 = init_graph(g10_4,rho,initial_infect_fract);
    
    g10_2 = stoch_sim_over_network(g10_2, mu,beta, 0, t_end, iters);
    g10_3 = stoch_sim_over_network(g10_3, mu,beta, 0, t_end, iters);
    g10_4 = stoch_sim_over_network(g10_4, mu,beta, 0, t_end, iters);
    
    V10_2(i,:) = [sum(g10_2.Nodes.Population), sum(g10_2.Nodes.I)./sum(g10_2.Nodes.Population)];
    V10_3(i,:) = [sum(g10_3.Nodes.Population), sum(g10_3.Nodes.I)./sum(g10_3.Nodes.Population)];
    V10_4(i,:) = [sum(g10_4.Nodes.Population), sum(g10_4.Nodes.I)./sum(g10_4.Nodes.Population)];

%     plot(sum(g10_2.Nodes.Population), sum(g10_2.Nodes.I)./sum(g10_2.Nodes.Population),'r*');
%     hold on
%     plot(sum(g10_3.Nodes.Population), sum(g10_3.Nodes.I)./sum(g10_3.Nodes.Population),'g*');
%     hold on
%     plot(sum(g10_4.Nodes.Population), sum(g10_4.Nodes.I)./sum(g10_4.Nodes.Population),'b*');
end
%%
figure;
plot3(1:6, V10_2(:,1), V10_2(:,2));
hold on
plot3(1:6, V10_3(:,1), V10_3(:,2));
hold on
plot3(1:6, V10_4(:,1), V10_4(:,2));

figure;
plot(V10_2(:,1), V10_2(:,2), 'r*', "DisplayName", "V = 10^2");
hold on
text(V10_2(:,1), V10_2(:,2), "\rho = " + rhos, "FontSize", 11);
hold on
plot(V10_3(:,1), V10_3(:,2), 'g*', "DisplayName", "V = 10^3");
hold on
text(V10_3(:,1), V10_3(:,2), "\rho = " + rhos, "FontSize", 11);
hold on
plot(V10_4(:,1), V10_4(:,2), 'b*', "DisplayName", "V = 10^4");
hold on
text(V10_4(:,1), V10_4(:,2), "\rho = " + rhos, "FontSize", 11);
legend("FontSize", 12);
xlabel("Initial population size", "FontSize", 18);
ylabel("Rel. infected population", "FontSize", 18);


%% many sims no susp transfer
rho_interval = [1.5,3];
rho_iters = 6;
rhos = linspace(rho_interval(1), rho_interval(2), rho_iters);

mu = 2;
beta = 1;
iters = 2000;
t_end = 100;
initial_infect_fract = 0.1;

V10_2 = zeros(rho_iters, 2);
V10_3 = zeros(rho_iters, 2);
V10_4 = zeros(rho_iters, 2);

for i = 1:rho_iters
    rho = rhos(i);

    g10_2 = init_graph(g10_2,rho,initial_infect_fract);
    g10_3 = init_graph(g10_3,rho,initial_infect_fract);
    g10_4 = init_graph(g10_4,rho,initial_infect_fract);
    
    g10_2 = stoch_sim_over_network(g10_2, mu,beta, 0, t_end, iters, "NoSusp");
    g10_3 = stoch_sim_over_network(g10_3, mu,beta, 0, t_end, iters, "NoSusp");
    g10_4 = stoch_sim_over_network(g10_4, mu,beta, 0, t_end, iters, "NoSusp");
    
    V10_2(i,:) = [sum(g10_2.Nodes.Population), sum(g10_2.Nodes.I)./sum(g10_2.Nodes.Population)];
    V10_3(i,:) = [sum(g10_3.Nodes.Population), sum(g10_3.Nodes.I)./sum(g10_3.Nodes.Population)];
    V10_4(i,:) = [sum(g10_4.Nodes.Population), sum(g10_4.Nodes.I)./sum(g10_4.Nodes.Population)];

%     plot(sum(g10_2.Nodes.Population), sum(g10_2.Nodes.I)./sum(g10_2.Nodes.Population),'r*');
%     hold on
%     plot(sum(g10_3.Nodes.Population), sum(g10_3.Nodes.I)./sum(g10_3.Nodes.Population),'g*');
%     hold on
%     plot(sum(g10_4.Nodes.Population), sum(g10_4.Nodes.I)./sum(g10_4.Nodes.Population),'b*');
end
%% many sims GROUP SPREAD
rho = 2;
mu = 0.5;
beta = 1;
iters = 200;
t_end = 500;
initial_infect_fract = 0.2;
close all

g10_2 = init_graph(g10_2,rho,initial_infect_fract);
g10_3 = init_graph(g10_3,rho,initial_infect_fract);
g10_4 = init_graph(g10_4,rho,initial_infect_fract);

g10_2 = stoch_sim_over_network(g10_2, mu,beta, 0, t_end, iters, "", "Group");
g10_3 = stoch_sim_over_network(g10_3, mu,beta, 0, t_end, iters, "", "Group");
g10_4 = stoch_sim_over_network(g10_4, mu,beta, 0, t_end, iters, "", "Group");

plot(sum(g10_2.Nodes.Population), sum(g10_2.Nodes.I)./sum(g10_2.Nodes.Population),'*');
hold on
plot(sum(g10_3.Nodes.Population), sum(g10_3.Nodes.I)./sum(g10_3.Nodes.Population),'*');
hold on
plot(sum(g10_4.Nodes.Population), sum(g10_4.Nodes.I)./sum(g10_4.Nodes.Population),'*');

%% many sims different rho
rho_interval = [1.5,3];
rho_iters = 6;
rhos = linspace(rho_interval(1), rho_interval(2), rho_iters);

mu = 0.5;
beta = 1;
iters = 2000;
t_end = 100;
initial_infect_fract = 0.1;

V10_2 = zeros(rho_iters, 2);
V10_3 = zeros(rho_iters, 2);
V10_4 = zeros(rho_iters, 2);

for i = 1:rho_iters
    rho = rhos(i);

    g10_2 = init_graph(g10_2,rho,initial_infect_fract);
    g10_3 = init_graph(g10_3,rho,initial_infect_fract);
    g10_4 = init_graph(g10_4,rho,initial_infect_fract);
    
    g10_2 = stoch_sim_over_network(g10_2, mu,beta, 0, t_end, iters, "NoSusp", "Group");
    g10_3 = stoch_sim_over_network(g10_3, mu,beta, 0, t_end, iters, "NoSusp", "Group");
    g10_4 = stoch_sim_over_network(g10_4, mu,beta, 0, t_end, iters, "NoSusp", "Group");
    
    V10_2(i,:) = [sum(g10_2.Nodes.Population), sum(g10_2.Nodes.I)./sum(g10_2.Nodes.Population)];
    V10_3(i,:) = [sum(g10_3.Nodes.Population), sum(g10_3.Nodes.I)./sum(g10_3.Nodes.Population)];
    V10_4(i,:) = [sum(g10_4.Nodes.Population), sum(g10_4.Nodes.I)./sum(g10_4.Nodes.Population)];
end