clear

dt_002 = load("data/0064_r0.500_T1.000_alpha0.00_dt002");
dt_005 = load("data/0064_r0.500_T1.000_alpha0.00_dt005");
dt_010 = load("data/0064_r0.500_T1.000_alpha0.00_dt010");
dt_014 = load("data/0064_r0.500_T1.000_alpha0.00_dt014");
dt_016 = load("data/0064_r0.500_T1.000_alpha0.00_dt016");

N = 64;

% note we discard last as that is very divergent
xdata = [0.002, 0.005, 0.010, 0.014];
ydata = [mean(dt_002(:,2)), mean(dt_005(:,2)), mean(dt_010(:,2)), mean(dt_014(:,2))];
edata_std = [std(dt_002(:,2)), std(dt_005(:,2)), std(dt_010(:,2)), std(dt_014(:,2))];
edata_standard_error = edata_std ./sqrt(N);


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