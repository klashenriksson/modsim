num_countries = 100;
world_graph = gen_world(num_countries, -2.5, -2.5);
rho = 200;
mu = 2;
beta = 1;
t_end = 1;
iters = 200;
initial_infect_fract = 0.05;

world_graph = init_graph(world_graph,rho,initial_infect_fract);
world_graph = stoch_sim_over_network(world_graph, mu, beta, 0, t_end, iters);
%%
plot_world_graph(world_graph);