clear

data_dt_002 = [-0.299671, 5.67969e-05];
data_dt_005 = [-0.296918, 0.000160269];
data_dt_010 = [-0.299772, 0.000192871];
data_dt_014 = [-0.292456, 0.00105202];
data_dt_020 = [4.13966e+14, 3.1706e+06];

N = 64;
NBlocks = 50;

% note we discard last as that is very divergent
xdata = [0.002, 0.005, 0.010, 0.014];
ydata = [data_dt_002(1), data_dt_005(1), data_dt_010(1), data_dt_014(1)];
edata_standard_error = [data_dt_002(2), data_dt_005(2), data_dt_010(2), data_dt_014(2)];
edata_std = edata_standard_error * sqrt(NBlocks-1);


figure(1);
% use halve edata as the total width of error bar is 2*err originally
errorbar(xdata,ydata,edata_standard_error./2);
xlabel("Time", "FontSize", 18);
ylabel("Energy", "FontSize", 18);
legend("E_{tot} w/ stderr as errorbar", "FontSize", 14);

figure(2);
errorbar(xdata,ydata,edata_std./2);
xlabel("Time", "FontSize", 18);
ylabel("Energy", "FontSize", 18);
legend("E_{tot} w/ std as errorbar", "FontSize", 14);