clear
close all
alpha = 10.0;

% potential energy, standard error
data_mc = [-0.813144, 0.00095619];
data_brown_dt_0_001 = [-0.812391, 0.00230368];
data_brown_dt_0_002 = [-0.807627, 0.00229925];
data_brown_dt_0_003 = [-0.805504, 0.00258021];
data_brown_dt_0_004 = [-0.801663, 0.00231858];

dt = [0.001, 0.002, 0.003, 0.004];

xdata = [0, dt./alpha];
ydata = [data_mc(1), data_brown_dt_0_001(1), data_brown_dt_0_002(1), ...
    data_brown_dt_0_003(1), data_brown_dt_0_004(1)];
edata = [data_mc(2),data_brown_dt_0_001(2), data_brown_dt_0_002(2), ...
    data_brown_dt_0_003(2), data_brown_dt_0_004(2)];

figure(1);
% use halve edata as the total width of error bar is 2*err originally
errorbar(xdata,ydata,edata./2);
xlabel("\Delta_t/\alpha", "FontSize", 18);
ylabel("Energy", "FontSize", 18);
legend("E_{pot} w/ standard error as errorbar", "FontSize", 14);
