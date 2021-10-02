close all
clear

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
%%
rho_interval = [0.5, 2];
rho_iters = 6;
rhos = linspace(rho_interval(1), rho_interval(2), rho_iters);

mu = 0.5;
gamma = 2;
beta = 1;
iters = 2000;
t_end = 1;
initial_infect_fract = 0.05;

V10_2 = zeros(rho_iters, 2);
V10_3 = zeros(rho_iters, 2);
V10_4 = zeros(rho_iters, 2);

for i = 1:rho_iters
    rho = rhos(i);

    g10_2_i = siv_init_graph(g10_2,rho,initial_infect_fract);
    g10_3_i = siv_init_graph(g10_3,rho,initial_infect_fract);
    g10_4_i = siv_init_graph(g10_4,rho,initial_infect_fract);
    
    g10_2_i = siv_stoch_sim_over_network(g10_2_i, mu,beta, gamma, 0, t_end, iters*10, "Lockdown", "Group");
    g10_3_i = siv_stoch_sim_over_network(g10_3_i, mu,beta, gamma, 0, t_end, iters, "Lockdown", "Group");
    g10_4_i = siv_stoch_sim_over_network(g10_4_i, mu,beta, gamma, 0, t_end, iters, "Lockdown", "Group");
    
    V10_2(i,:) = [sum(g10_2_i.Nodes.Population), sum(g10_2_i.Nodes.I)./sum(g10_2_i.Nodes.Population)];
    V10_3(i,:) = [sum(g10_3_i.Nodes.Population), sum(g10_3_i.Nodes.I)./sum(g10_3_i.Nodes.Population)];
    V10_4(i,:) = [sum(g10_4_i.Nodes.Population), sum(g10_4_i.Nodes.I)./sum(g10_4_i.Nodes.Population)];
end

%% plotting
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
