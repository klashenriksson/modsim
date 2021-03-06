%% task 1 a
close all
x0 = [600;4];
a = 0.4;
b = 0.001;
c = 0.001;
d = 0.9;
ode = @(t,x) [
    a*x(1) - b*x(1)*x(2);
    c*x(1)*x(2) - d*x(2);
    ];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-4);
disp("ODE45");
[t,x] = ode45(ode, [0, 1000], x0, options);
disp("ODE23S");
SOL23s = ode23s(ode, [0, 1000], x0, options);

plot(t,x);
legend("Rabbits", "Foxes");
xlabel("Time", 'FontSize', 18);
ylabel("Population Size", "FontSize", 18);

figure;
plot(x(:,1), x(:,2));
xlabel("Rabbits", "FontSize", 18);
ylabel("Foxes", "FontSize", 18);

%% task 1 b
x0 = [600;400];
a = 0.4;
b = 0.001;
c = 0.001;
d = 0.9;
alpha = d/a;

x0 = [x0(1)*c/d;x0(2)*b/a];
ode = @(t,x) [
    x(1)*(1-x(2));
    alpha*x(2)*x(1) - alpha*x(2);
    ];
tspan = 1:0.2:100;
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-4);
[t,x] = ode45(ode, tspan, x0, options);

plot(t,x);
legend("Rabbits", "Foxes");

figure;
plot(x(:,1), x(:,2));
xlabel("Rabbits");
ylabel("Foxes"),

%% tasb 1 c
close all
tspan = 1:0.2:10000;
x0 = [600;400];
a = 0.4;
b = 0.001;
c = 0.001;
d = 0.9;
r_grass = 5355*c/d;
r_sat = 11070*c/d;
alpha = d/a;

x0 = [x0(1)*c/d;x0(2)*b/a];
x0 = [1.0885;0.8894];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-4);
ode = @(t,x) [
    x(1)*(1-x(1)./r_grass) - (x(1)*x(2))/(1+(x(1)./r_sat));
    alpha*x(1)*x(2)/(1+(x(1)./r_sat)) - alpha*x(2);
    ];
[t,x] = ode45(ode, tspan, x0, options);

plot(t,x);
legend("Rabbits", "Foxes");

figure;
plot(x(:,1), x(:,2));
xlabel("Rabbits", "FontSize", 18);
ylabel("Foxes", "FontSize", 18);

figure;
for i = 1:20
    r_grass = 5.8 + i*0.05;
    r_sat = 10.8 + i*0.5;
    ode = @(t,x) [
    x(1)*(1-x(1)./r_grass) - (x(1)*x(2))/(1+(x(1)./r_sat));
    alpha*x(1)*x(2)/(1+(x(1)./r_sat)) - alpha*x(2);
    ];
[t,x] = ode45(ode, tspan, x0, options);
plot(x(:,1), x(:,2), 'DisplayName', "R_{grass} = " + r_grass + " R_{sat} = " + r_sat);
hold on
end

%% task 1 d
close all
x0 = [600;400;20];
a = 0.4;
b = 0.001;
c = 0.001;
d = 0.9;
e = c*5.5;
f = c/3;
g = 0.7;
h = 0.1;
ode = @(t,x) [
    a*x(1) - b*x(1)*x(2) - e*x(1)*x(3);
    c*x(1)*x(2) - d*x(2) - f*x(2)*x(3);
    e*x(1)*x(3) + f*x(2)*x(3) - g*x(3) - h*x(3)*x(3);
    ];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-4);
[t,x] = ode45(ode, [0, 1000], x0, options);

plot(t,x);
legend("Rabbits", "Foxes", "Wolves");
xlabel("Time", 'FontSize', 18);
ylabel("Population Size", "FontSize", 18);

figure;
plot3(x(:,1), x(:,2), x(:,3));
xlabel("Rabbits", "FontSize", 18);
ylabel("Foxes", "FontSize", 18);
zlabel("Wolves", "FontSize", 18);

%% bonus 1
a = 0.4;
b = 0.001;
d = 0.9;

ode = @(t,x) [
    a*x(1) - b*x(1)*x(2);
    b*x(1)*x(2) - d*x(2);
    ];

iters = 400000;
x0 = [600;400];
[t,x_sol] = ode45(ode, [0, 400], x0);

x = zeros(iters+1,2);
x(1,:) = x0;
time = zeros(iters+1,2);
time(1) = 0;
for i = 1:iters
    index = i+1;
    k = [a*x(i,1),b*x(i,1)*x(i,2),d*x(i,2)];
    r0 = rand();
    dt = -(1/sum(k))*log(r0);
    r1 = rand();
    
    time(index) = time(i) + dt;
    x(index,:) = x(i,:);
    
    if r1 < k(1)/sum(k)
        x(index,1) = x(i,1) + 1;
    elseif r1 < (k(1) + k(2))/sum(k)
        x(index,2) = x(i,2) + 1;
        x(index,1) = x(i,1) - 1;
    else
       x(index,2) = x(i,2) - 1; 
    end
end
plot(time, x(:,1));
hold on
plot(time, x(:,2));
hold on
plot(t,x_sol, 'MarkerSize', 3);
legend("Stochastic Rabbits", "Stochastic Foxes", "ODE Rabbits", "ODE Foxes");
figure;
plot(x(:,1),x(:,2));
hold on
plot(x_sol(:,1), x_sol(:,2));
legend("Stochastic", "ODE");

%% bonus 2 - extinct if b*x(1,1) < 1??
a = 0.4;
b = 0.001;
d = 0.9;

iters = 500000;
x0 = [600;400];

x_avg = zeros(iters+1,2);
avgs = 20;
for avg = 1:avgs
    x = zeros(iters+1,2);
    x(1,:) = x0;
    for i = 1:iters
        index = i+1;
        k = [a*x(i,1),b*x(i,1)*x(i,2),d*x(i,2)];
        r0 = rand();
        dt = -(1/sum(k))*log(r0);
        r1 = rand();

        x(index,:) = x(i,:);

        if r1 < k(1)/sum(k)
            x(index,1) = x(i,1) + 1;
        elseif r1 < (k(1) + k(2))/sum(k)
            x(index,2) = x(i,2) + 1;
            x(index,1) = x(i,1) - 1;
        else
           x(index,2) = x(i,2) - 1; 
        end
    end
    x_avg = x_avg + x;
end
x_avg = x_avg ./ avgs;
plot(1:iters+1, x_avg);
legend("Rabbits", "Foxes");
figure;
plot(x_avg(:,1),x_avg(:,2));

%% heatmap
a = 0.4;
b = 0.001;
d = 0.9;

iters = 500000;
start_rabbits = 100;
max_rabbits = 1000;
start_foxes = 100;
max_foxes = 1000;
delta_rabbit = 50;
delta_fox = 50;

rabbit_iter_count = ceil((max_rabbits-start_rabbits)/delta_rabbit);
fox_iter_count = ceil((max_foxes-start_foxes)/delta_fox);
x0 = zeros(rabbit_iter_count, fox_iter_count, 2);
xvalues = cell(rabbit_iter_count,1);
yvalues = cell(fox_iter_count,1);
for r = 1:rabbit_iter_count
    for c = 1:fox_iter_count
        rabbit = start_rabbits + delta_rabbit.*(r-1);
        foxes = start_foxes + delta_fox*(c-1);
        x0(r,c,:) = [
            rabbit,
            foxes;
        ];
    end
end

for r = 1:rabbit_iter_count
    xvalues{r} = start_rabbits + delta_rabbit.*(r-1);
end

for f = 1:fox_iter_count
    yvalues{f} = start_foxes + delta_fox.*(f-1);
end

size_x0 = size(x0);
fox_die_matrix = zeros(size_x0(1), size_x0(2));

avgs = 10;

for r = 1:size_x0(1)
    for c = 1:size_x0(2)
        for avg = 1:avgs
            x = zeros(iters+1,2);
            x(1,:) = x0(r,c,:);
            
            end_index = 0;
            for i = 1:iters
                index = i+1;
                k = [a*x(i,1),b*x(i,1)*x(i,2),d*x(i,2)];
                r0 = rand();
                dt = -(1/sum(k))*log(r0);
                r1 = rand();

                x(index,:) = x(i,:);

                if r1 < k(1)/sum(k)
                    x(index,1) = x(i,1) + 1;
                elseif r1 < (k(1) + k(2))/sum(k)
                    x(index,2) = x(i,2) + 1;
                    x(index,1) = x(i,1) - 1;
                else
                   x(index,2) = x(i,2) - 1; 
                end
                
                end_index = index;
                if x(index,2) < 1e-7
                    break;
                end
            end
            
            if x(end_index,2) < 1e-7
               fox_die_matrix(r,c) = fox_die_matrix(r,c) + 1; 
            end
        end
        
        fox_die_matrix(r,c) = fox_die_matrix(r,c) / avgs;
    end
end

h = heatmap(xvalues, yvalues, fox_die_matrix);
h.Title = "Fox Death Heatmap";
h.XLabel = "Rabbits";
h.YLabel = "Foxes";

%% task 1b - what is alpha?
x0 = [600;400];
a = 0.4;
b = 0.001;
c = 0.001;
d = 0.9;
alpha = 1;%d/a;

x0 = [x0(1)*c/d;x0(2)*b/a];
x0 = [0.6667; 1.000];
ode = @(t,x) [
    x(1)*(1-x(2));
    alpha*x(2)*x(1) - alpha*x(2);
    ];
tspan = 1:0.2:100;
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-4);
[t,x] = ode45(ode, tspan, x0, options);

plot(t,x);
legend("Rabbits", "Foxes");

figure;
plot(x(:,1), x(:,2));
xlabel("Rabbits");
ylabel("Foxes"),