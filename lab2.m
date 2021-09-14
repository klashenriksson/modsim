a = 1;
b = 2;
x0 = [1;1];
FP = [1;b/a];

tspan = 0:0.2:500;
u_fun = @(x,y) 1 - (b+1).*x + a.*x.*x.*y;
v_fun = @(x,y) b.*x - a.*x.*x.*y;
ode = @(t,x) [
    u_fun(x(1),x(2));
    v_fun(x(1),x(2));
    ];
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-4);
[t,x] = ode45(ode,tspan,x0,options);
plot(t,x);
figure;
plot(x(:,1),x(:,2));

figure;
arrowcount = 10;
end_coord = 2*(b/a);
[xspan,yspan] = meshgrid(0:end_coord/arrowcount:end_coord,0:end_coord/arrowcount:end_coord);
u = u_fun(xspan,yspan);
v = v_fun(xspan,yspan);
quiver(xspan,yspan,u,v);