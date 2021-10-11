clear
close all

alpha_0_01 = load("data/0064_r0.600_T1.000_alpha0.01_dt010");
alpha_0_10 = load("data/0064_r0.600_T1.000_alpha0.10_dt010");
alpha_1_00 = load("data/0064_r0.600_T1.000_alpha1.00_dt010");

std_var = [std(alpha_0_01(:,2)), std(alpha_0_10(:,2)), std(alpha_1_00(:,2))];
alpha = [0.01, 0.10, 1.00];

figure(1);
plot(alpha, std_var);
xlabel("\alpha", "FontSize", 18);
ylabel("Std. deviation", "FontSize", 18);

%relate energy fluctuations to m/alpha..?
figure(2);
plot(alpha_0_01(:,1), alpha_0_01(:,2), "DisplayName", "\alpha = 0.01");
hold on
plot(alpha_0_10(:,1), alpha_0_10(:,2), "DisplayName", "\alpha = 0.10");
hold on
plot(alpha_1_00(:,1), alpha_1_00(:,2), "DisplayName", "\alpha = 1.00");
xlabel("Time", "FontSize", 18);
ylabel("Energy", "FontSize", 18);
legend("FontSize", 14);