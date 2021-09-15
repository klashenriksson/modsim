clear

alpha = 77.27;
beta = 8.375e-6;
gamma = 1.161;
x0 = [1;2;3];
ode = @(t,x) [
    alpha.*(x(2) + x(1).*(1-beta.*x(1)-x(2)));
    (1./alpha).*(x(3)-(1+x(1)).*x(2));
    gamma.*(x(1)-x(3))
];
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-4);
tspan = 1:0.002:500;
[t,x] = ode15s(ode, tspan, x0, options);
plot(t,x);