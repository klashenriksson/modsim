close all
a = 1;
b = 2;
x0 = [0.7;2.5];
FP = [1;b/a];

tspan = 0:0.2:50000;
u_fun = @(x,y) 1 - (b+1).*x + a.*x.*x.*y;
v_fun = @(x,y) b.*x - a.*x.*x.*y;
ode = @(t,x) [
    u_fun(x(1),x(2));
    v_fun(x(1),x(2));
    ];
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-4);
[t,x] = ode45(ode,[0 400],x0,options);
plot(t,x);
xlabel("Time", "FontSize", 18);
ylabel("Quantity of chemical", "FontSize", 18);
figure;
plot(x(:,1),x(:,2), 'LineWidth', 1.2);

hold on
arrowcount = 20;
end_coord = 2*(b/a);
[xspan,yspan] = meshgrid(0:end_coord/arrowcount:end_coord,0:end_coord/arrowcount:end_coord);
u = u_fun(xspan,yspan);
v = v_fun(xspan,yspan);
velo_plot = quiver(xspan,yspan,u./(sqrt(u.^2 + v.^2)),v./(sqrt(u.^2 + v.^2)));

nullcline_1 = @(x) (b+1)./(a.*x) - 1 ./ (a.*x.*x);
nullcline_2 = @(x) b./(a.*x);

hold on
fplot(nullcline_1, [0.01 3], 'LineWidth', 2);
hold on
fplot(nullcline_2, [0.01 3], 'LineWidth', 2);
axis([0 3 -4 4]);
xlabel("Quantity of X chemical", "FontSize", 18);
ylabel("Quantity of Y chemical", "FontSize", 18);

T = @(a,b) 4*pi/sqrt(abs(1-2*a-2*b+a*a-2*a*b+b*b));
% check this T for a = 1 b = 2, compare to oscillation plots