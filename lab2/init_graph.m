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