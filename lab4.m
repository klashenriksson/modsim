clear
close all

x0 = [-8;8;27];
beta = 8/3;
rho = 50;
sigma = 10;
end_t = 20;

ode = @(t,x) [
    sigma.*(x(2)-x(1));
    x(1).*(rho-x(3)) - x(2);
    x(1).*x(2) - beta.*x(3);
];

jacobian = @(t,x) [
    -sigma, sigma, 0;
    rho - x(3), -1, -x(1);
    x(2), x(1), -beta
];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-14, 'Jacobian', jacobian, 'Stats', 'on');
SOL45 = ode45(ode,[0 end_t], x0, options);
SOL23s = ode23s(ode, [0 end_t], x0, options);
SOL113 = ode113(ode, [0 end_t], x0, options);

plot_shit(SOL45, jacobian);
plot_shit(SOL23s, jacobian);
plot_shit(SOL113, jacobian);

figure;
x1_45 = SOL45.y(1,:);
plot(SOL45.x, x1_45);
hold on
x1_23s = SOL23s.y(1,:);
plot(SOL23s.x, x1_23s);
hold on
x1_113 = SOL113.y(1,:);
plot(SOL113.x, x1_113);

%% part 2
close all
clear
delta = 1e-15;
x0_1 = [3;pi;pi];
x0_2 = x0_1 + [delta; 0; 0];
beta = 8/3;
rho = 50;
sigma = 10;
end_t = 500;
dt = 0.0002;

ode = @(t,x) [
    sigma.*(x(2)-x(1));
    x(1).*(rho-x(3)) - x(2);
    x(1).*x(2) - beta.*x(3);
];

jacobian = @(t,x) [
    -sigma, sigma, 0;
    rho - x(3), -1, -x(1);
    x(2), x(1), -beta
];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-14, 'Jacobian', jacobian, 'Stats', 'on');
SOL1 = ode45(ode,0:dt:end_t, x0_1, options);
SOL11 = ode45(ode,0:dt:end_t, [0.5,-6,14], options);
SOL2 = ode45(ode,0:dt:end_t, x0_2, options);

x1 = SOL1.y;
x2 = SOL2.y;

x1 = x1(:,1:min(size(x1,2),size(x2,2)));
x2 = x2(:,1:min(size(x1,2),size(x2,2)));

time = SOL1.x;
time = time(1:min(size(time,2), size(x1,2)));

xdiff = x1-x2;
length_vec = dot(xdiff,xdiff);
norm_len = sqrt(dot(x1-x2,x1-x2));
log_len = log(norm_len);

figure;
plot(time, norm_len);
figure;
plot(time, log(norm_len));

figure;
plot3(SOL1.y(1,:),SOL1.y(2,:),SOL1.y(3,:));
hold on
plot3(SOL2.y(1,:),SOL2.y(2,:),SOL2.y(3,:));

figure;
plot3(SOL1.y(1,:),SOL1.y(2,:),SOL1.y(3,:));
hold on
plot3(SOL11.y(1,:),SOL11.y(2,:),SOL11.y(3,:));

figure;
x3 = SOL1.y(3,:);
pks = findpeaks(x3);
plot(pks(1:end-1),pks(2:end),'*');

figure;
plot(pks);

function plot_shit(SOL, jacobian)
time = SOL.x;
timesteps = diff(time);
x1 = SOL.y(1,:);
x2 = SOL.y(2,:);
x3 = SOL.y(3,:);

figure;
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

end
