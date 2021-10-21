function [world_graph, world] = gen_world(num_countries, gamma_countries, gamma_cities)
    country_adj_mat = gen_network(num_countries, gamma_countries);
    country_graph = graph(country_adj_mat);
    total_connections = numedges(country_graph);
    
    world.cntry_graph = country_graph;
    world.city_nets = cell(num_countries, 1);

    deg_factor = 50;

    num_cntries = country_graph.numnodes;
    degs = degree(country_graph);

    total_nodes = 0;
    max_city_count = 0;

    city_count = zeros(num_countries, 1);
    for n = 1:num_cntries
        num_cities = degs(n)*deg_factor;
        max_city_count = max(num_cities, max_city_count);
        city_count(n) = num_cities;

        city_net = gen_network(num_cities, gamma_cities);
        city_graph = graph(city_net);
        total_connections = total_connections + numedges(city_graph);
        total_nodes = total_nodes + num_cities;

        world.city_nets(n) = {city_graph};
    end

    %recombine
    world_edges = zeros(total_connections,2);
    country_offset = 0;
    country_offsets = zeros(num_cntries, 1);
    edge_index = 1;
    for n = 1:num_cntries
        city_net = world.city_nets(n);
        city_net = city_net{1};
        edges = city_net.Edges;
        edge_count = size(edges, 1);
        num_nodes = city_net.numnodes;
        for edge = 1:edge_count
            node1 = edges.EndNodes(edge,1);
            node2 = edges.EndNodes(edge,2);

            n1_id = node1 + country_offset;
            n2_id = node2 + country_offset;

            world_edges(edge_index,:) = [n1_id,n2_id];
            edge_index = edge_index + 1;
        end
        country_offsets(n) = country_offset;
        country_offset = country_offset + num_nodes;
        
    end

    country_data.offsets = country_offsets;
    country_data.city_count = city_count;
    country_data.n_countries = num_countries;

    % connect the largest cities
    cntry_edges = world.cntry_graph.Edges.EndNodes;
    edge_count = size(cntry_edges, 1);
    for edge = 1:edge_count
        cntry_1_id = cntry_edges(edge, 1);
        cntry_2_id = cntry_edges(edge, 2);
        city_net_1 = world.city_nets(cntry_1_id);
        city_net_1 = city_net_1{1};
        city_net_2 = world.city_nets(cntry_2_id);
        city_net_2 = city_net_2{1};

        [~, biggest_node_city_1] = max(degree(city_net_1));
        [~, biggest_node_city_2] = max(degree(city_net_2));

        city_1_node_id = biggest_node_city_1 + country_offsets(cntry_1_id);
        city_2_node_id = biggest_node_city_2 + country_offsets(cntry_2_id);

        world_edges(edge_index, :) = [city_1_node_id, city_2_node_id];
        edge_index = edge_index + 1;
    end

    weights = ones(size(world_edges,1),1);

    edge_table = table(world_edges, weights, 'VariableNames', {'EndNodes', 'Weight'});
    world_graph = graph(edge_table);

end
