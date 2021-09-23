clear
close all

x0 = [-8;8;27];
beta = 8/3;
rho = 50;
sigma = 10;
start_t = 0;
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
SOL45 = ode45(ode,[start_t end_t], x0, options);
SOL23s = ode23s(ode, [start_t end_t], x0, options);
SOL113 = ode113(ode, [start_t end_t], x0, options);

plot_shit(SOL45, jacobian);
%plot_shit(SOL23s, jacobian);
%plot_shit(SOL113, jacobian);

figure;
x1_45 = SOL45.y(1,:);
plot(SOL45.x, x1_45);
hold on
x1_23s = SOL23s.y(1,:);
plot(SOL23s.x, x1_23s);
hold on
x1_113 = SOL113.y(1,:);
plot(SOL113.x, x1_113);
legend("ode45", "ode23s", "ode113");
xlabel("Time", "FontSize", 18);
ylabel("X_1", "FontSize", 18);

%% part 2
close all
clear
delta = 1e-15;
x0_1 = [3;pi;pi];
x0_2 = x0_1 + [delta; 0; 0];
beta = 8/3;
rho = 50;
sigma = 10;
end_t = 20;
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
xlabel("Time", "FontSize", 18);
ylabel("s.u", "FontSize", 18);
legend("\delta");
figure;
plot(time, log(norm_len));
xlabel("Time", "FontSize", 18);
ylabel("s.u", "FontSize", 18);
legend("ln(\delta)");

figure;
plot3(SOL1.y(1,:),SOL1.y(2,:),SOL1.y(3,:));
hold on
plot3(SOL2.y(1,:),SOL2.y(2,:),SOL2.y(3,:));
legend("Init. cond. x_0", "Init. cond. x_0 + \delta");
xlabel("x_1", "FontSize", 18);
ylabel("x_2", "FontSize", 18);
zlabel("x_3", "FontSize", 18);

%figure;
%plot3(SOL1.y(1,:),SOL1.y(2,:),SOL1.y(3,:));
%hold on
%plot3(SOL11.y(1,:),SOL11.y(2,:),SOL11.y(3,:));
curve_1 = @(x) 0.129.*x.^2 - 13.29.*x + 397.7;
curve_2 = @(x) 0.09815.*x.^2 - 15.65.*x + 677.5;

x3 = SOL1.y(3,:);

[pks, locs] = findpeaks(x3);
pks_1 = pks(1:end-1);
pks_offsetted = pks(2:end);
figure;
plot(pks(1:end-1),pks(2:end),'*');
hold on
fplot(curve_1, [53 65], "LineWidth", 2);
hold on
fplot(curve_2, [64, 80], "LineWidth", 2);
xlabel("x_3", "FontSize", 18);
ylabel("x_3", "FontSize", 18);
legend("(peak x_{n+1}, peak x_n)", "Quad. approx", "Quad. approx");

approx_peaks = zeros(1, size(pks,2));

approx_peaks(1) = pks(1);
for i = 2:size(approx_peaks, 2)
   old_peak = approx_peaks(i-1);
   if old_peak < 65
       approx_peaks(i) = curve_1(old_peak);
   elseif old_peak >= 65 && old_peak < 89
       approx_peaks(i) = curve_2(old_peak);
   end
end
figure;
plot(time(locs), pks);
hold on
plot(time(locs), approx_peaks);
xlabel("Time", "FontSize", 18);
ylabel("X_3", "FontSize", 18);
legend("Real local peaks", "Approximated local peaks");

figure;
plot(pks);

%% FP test
clear
close all

beta = 8/3;
rho = 50;
sigma = 10;
end_t = 1.1;

x0_3 = rho - 1;
x0_1 = sqrt(beta*(x0_3));
x0 = [x0_1; x0_1; x0_3] + [0.01;0.01;0.01]; %[-8; 8; 27];

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
plot_shit(SOL45, jacobian);

function plot_shit(SOL, jacobian)
time = SOL.x;
timesteps = diff(time);
x1 = SOL.y(1,:);
x2 = SOL.y(2,:);
x3 = SOL.y(3,:);

figure;
plot3(x1,x2,x3);
xlabel("X_1", "FontSize", 18);
ylabel("X_2", "FontSize", 18);
zlabel("X_3", "FontSize", 18);
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
