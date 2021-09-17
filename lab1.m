%% task 1 a
x0 = [600;400];
a = 0.4;
b = 0.001;
c = 0.001;
d = 0.9;
ode = @(t,x) [
    a*x(1) - b*x(1)*x(2);
    c*x(1)*x(2) - d*x(2);
    ];
tspan = 0:0.2:100;

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-4);
[t,x] = ode45(ode, tspan, x0, options);

plot(t,x);
legend("Rabbits", "Foxes");

figure;
plot(x(:,1), x(:,2));
xlabel("Rabbits");
ylabel("Foxes");

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
r_grass = 1.8;
r_sat = 27.5;
tspan = 1:0.2:10000;
x0 = [600;300];
a = 0.4;
b = 0.001;
c = 0.001;
d = 0.9;
alpha = d/a;

x0 = [x0(1)*c/d;x0(2)*b/a];

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
xlabel("Rabbits");
ylabel("Foxes");

%% task 1 d
close all
r_grass = 1.8;
r_sat = 27.5;
tspan = 1:0.2:200;
x0 = [600;300];
a = 0.4;
b = 0.001;
c = 0.001;
d = 0.9;
alpha = d/a;

x0 = [x0(1)*c/d;x0(2)*b/a;0.05];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-4);
ode = @(t,x) [
    x(1)*(1-x(1)./r_grass) - (x(1)*x(2))/(1+(x(1)./r_sat) - alpha*6*x(1)*x(3));
    alpha*x(1)*x(2)/(1+(x(1)./r_sat)) - alpha*x(2) - alpha*0.3*x(2)*x(3);
    6*alpha*x(1)*x(3) + alpha*0.3*x(2)*x(3) - 2*alpha*x(3);
    ];
[t,x] = ode45(ode, tspan, x0, options);

plot(t,x);
legend("Rabbits", "Foxes", "Wolves");
figure;

plot3(x(:,1),x(:,2),x(:,3));
xlabel("Rabbit");
ylabel("Foxes");
zlabel("Wolves");

%% bonus 1
a = 0.4;
b = 0.001;
d = 0.9;

iters = 696969;
x0 = [600;400];

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
legend("Rabbits", "Foxes");
figure;
plot(x(:,1),x(:,2));

%% bonus 2 - extinct if b*x(1,1) < 1??
a = 0.4;
b = 0.001;
d = 0.9;

iters = 500000;
x0 = [600;400];

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
plot(1:iters+1, x);
legend("Rabbits", "Foxes");
figure;
plot(x(:,1),x(:,2));

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