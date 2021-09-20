%FP = (S,I) = (my/beta, 0)
% x = (S,I)

mu = 0.8;
beta = 0.001;

ode = @(t,x) [
    mu*x(2) - beta*x(2)*x(1);
    -mu*x(2) + beta*x(1)*x(2);
];
x0 = [900;100];
tot_pop = sum(x0);
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-4, 'Stats', 'on');

SOL = ode45(ode, [0 100], x0, options);

S = SOL.y(1,:);
I = SOL.y(2,:);
rho_s = S./(S+I);
rho_i = I./(I+S);

time = SOL.x;

plot(time, [S;I]);
legend("S", "I");
xlabel("Time", "FontSize", 18);

figure;
plot(time, [rho_s;rho_i]);
legend("S/total\_pop", "I/total\_pop");
xlabel("Time", "FontSize", 18);
ylabel("Ratio", "FontSize", 18);
hold on

ones_time = ones(1,size(time,2)) ./ tot_pop;
plot(time, [ones_time.*mu./beta;ones_time.*(tot_pop-(mu/beta))]);

figure;
plot(S,I);