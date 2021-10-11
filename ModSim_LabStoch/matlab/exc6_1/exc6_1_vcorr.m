close all

%time N vcorr
data_alpha_0_01 = load("data/vcorr_0064_r0.200_T1.000_alpha0.01_dt010.txt");
data_alpha_0_10 = load("data/vcorr_0064_r0.200_T1.000_alpha0.10_dt010.txt");
data_alpha_1_00 = load("data/vcorr_0064_r0.200_T1.000_alpha1.00_dt010.txt");

plot(data_alpha_0_01(:,1), data_alpha_0_01(:,3), "DisplayName", "\alpha = 0.01");
hold on
plot(data_alpha_0_10(:,1), data_alpha_0_10(:,3), "DisplayName", "\alpha = 0.10");
hold on
plot(data_alpha_1_00(:,1), data_alpha_1_00(:,3), "DisplayName", "\alpha = 1.00");
xlabel("Time", "FontSize", 18);
ylabel("Velocity Correlation", "FontSize", 18);

legend();