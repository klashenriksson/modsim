clear
W = load("nets/pownet_N100_gamma-2.5.mat");
W = W.W;
W(W ~= 0) = 1;
g = graph(W);
plot(g,'Layout','force');

% avg population per node
rho = 1;

pop = zeros(g.numnodes,1);
for i = 1:g.numnodes
    population = rand() * rho + rho/2;
    pop(i) = population;
end
g.Nodes.Population = pop;

figure;
for i = 1:500
    new_pop = zeros(g.numnodes,1);
    
    for n = 1:g.numnodes
        curr_neighbors = neighbors(g,n);
        for p = 1:g.Nodes.Population(n)
            neighbor_index = randi(size(curr_neighbors,1));
            new_pop(curr_neighbors(neighbor_index)) = new_pop(curr_neighbors(neighbor_index)) + 1;
        end
    end
    
    g.Nodes.Population = new_pop;
    
    pop_per_deg = zeros(max(degree(g)),1);
    rho_k = zeros(1,size(pop_per_deg,1));
    for j = 1:g.numnodes
       deg = degree(g,j); 
       pop_per_deg(deg) = pop_per_deg(deg) + g.Nodes.Population(j);
    end

    for j = 1:size(rho_k,2)
       rho_k(j) = (pop_per_deg(j)./mean(degree(g)))*rho;
    end
    
    if mod(i,10) == 0
        hold on
        plot(i,rho_k(2),'*');
        %disp("I: " + i + " RHO_K: " + rho_k(3));
    end
end

pop_per_deg = zeros(max(degree(g)),1);
rho_k = zeros(1,size(pop_per_deg,1));
for i = 1:g.numnodes
   deg = degree(g,i); 
   pop_per_deg(deg) = pop_per_deg(deg) + g.Nodes.Population(i);
end

for i = 1:size(rho_k,2)
   rho_k(i) = (pop_per_deg(i)./mean(degree(g)))*rho;
end