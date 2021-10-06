function plot_world_graph(world_graph)
    total_pop = sum(world_graph.Nodes.Population);

    [~, max_infect_idx] = max(world_graph.Nodes.I);
    infect_fract = world_graph.Nodes.I./total_pop;

    figure;
    h = plot(world_graph, 'Layout', 'force', 'NodeCData', infect_fract);
    colorbar;
    %h = plot(world_graph, 'Layout', 'force');


    degs = degree(world_graph);
    max_deg = max(degs);
    min_deg = min(degs);
    for deg = 1:max_deg
        n_ids = find(degs == deg);

        scale_factor = (deg-min_deg)/(max_deg-min_deg);
        highlight(h,n_ids, 'MarkerSize', 2 + (3*scale_factor).^2);
    end
end