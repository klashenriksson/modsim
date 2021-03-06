%rho, T, pressure
close all
clear

data = load("data/exc6_3_data.txt");

data_T_1_00 = data(data(:,2) == 1.0,:);
data_T_0_90 = data(data(:,2) == 0.9,:);
data_T_0_80 = data(data(:,2) == 0.8,:);
data_T_0_70 = data(data(:,2) == 0.7,:);
data_T_0_60 = data(data(:,2) == 0.6,:);
data_T_0_50 = data(data(:,2) == 0.5,:);

data = [data_T_0_50, data_T_0_60, data_T_0_70, data_T_0_80, data_T_0_90, data_T_1_00];
for i = 1:6
    rho_idx = (i-1)*3 + 1;
    t_idx = rho_idx + 1;
    p_idx = rho_idx + 2;

    rhos = data(:,rho_idx);
    pressures = data(:,p_idx);

    plot(1./rhos, pressures, 'DisplayName', "T = " + data(1,t_idx));
    hold on
end

rhos = data_T_0_80(:,1);
p_ideal_T_1_0 = 1.0*rhos;
p_ideal_T_0_5 = 0.5*rhos;

plot(1./rhos, p_ideal_T_1_0, '--', 'DisplayName', "Ideal T = 1.0");
hold on
plot(1./rhos, p_ideal_T_0_5, "--", "DisplayName", "Ideal T = 0.5");
legend("FontSize", 12);

xlabel("1/\alpha", "FontSize", 18);
ylabel("Pressure", "FontSize", 18);