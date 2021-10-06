num_countries = 100;
[world_graph, country_data] = gen_world(num_countries, -2.5, -2.5);
rho = 200;
mu = 2;
gamma = 2;
beta = 1;
t_end = 1;
iters = 200;
initial_infect_fract = 0.05;

world_graph = siv_init_graph(world_graph,rho,initial_infect_fract);

is_vax_country = zeros(world_graph.numnodes,1);
max_cities = max(country_data.city_count);

offset = 0;
for i = 1:country_data.n_countries
    if rand() < 2
        for n = 1:country_data.city_count
            n_id = n + offset;
            is_vax_country(n_id) = 1;
        end
    end

    offset = offset + country_data.city_count(i);
end
world_graph.Nodes.IVC = is_vax_country;
%%
world_graph = siv_stoch_sim_over_network_world(world_graph, mu, beta, gamma, 0, t_end, iters, "", "");
%%
plot_world_graph(world_graph);

%% viz
num_countries = 10;
world_graph = gen_world(num_countries, -2.5, 2.5);

world_graph = init_graph(world_graph, 200, 0);
plot_world_graph(world_graph);