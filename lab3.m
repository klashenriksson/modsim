clear
close all

alpha = 77.27;
beta = 8.375e-6;
gamma = 1.161;
x0 = [1;2;3];
ode = @(t,x) [
    alpha.*(x(2) + x(1).*(1-beta.*x(1)-x(2)));
    (1./alpha).*(x(3)-(1+x(1)).*x(2));
    gamma.*(x(1)-x(3))
];
jacobian = @(t,x) [
    alpha+2.*alpha.*beta.*x(1)-alpha.*beta.*x(2),alpha-alpha.*x(1), 0;
    -x(2)./alpha, (-1 + x(1))./alpha, 1./alpha;
    gamma, 0, -gamma
];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-4, 'Jacobian', jacobian, 'Stats', 'on');
tspan = 1:0.002:500;
disp("Ode 45");
%[t45,x45] = ode45(ode, [0 500], x0, options);
disp("Ode 113");
%[t113,x113] = ode113(ode, [0 500], x0, options);
disp("Ode 23s");
%[t23s,x23s] = ode23s(ode, [0, 500], x0, options);
SOL23s = ode23s(ode, [0, 500], x0, options);
disp("Ode 15s");
%[t15s, x15s] = ode15s(ode, [0, 500], x0, options);
%%
time = SOL23s.x;
timesteps = diff(time);
x1 = SOL23s.y(1,:);
x2 = SOL23s.y(2,:);
x3 = SOL23s.y(3,:);
plot3(x1,x2,x3);
figure;
plot(time(1:end-1), timesteps./max(timesteps));
hold on
plot(time, x1./max(x1));
hold on
plot(time, x2./max(x2));
hold on
plot(time, x3./max(x3));
legend("Time step size", "X1 conc.", "X2 conc.", "X3 conc.");

condNum = zeros(size(time,2),1);
for i = 1:size(time,2)
    jac = jacobian(1, [x1(i), x2(i), x3(i)]);
    condNum(i) = cond(jac);
end