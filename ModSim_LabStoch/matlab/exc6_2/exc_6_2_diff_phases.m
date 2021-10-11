clear
close all

data = [];
for i = 1:10
    suffix = "data/0200_r0.250_T";
    T = 1.000 - (i-1)*0.1;
    T_str = sprintf("%1.3f", T);
    filename = suffix + T_str;

    data = [data, load(filename)];
end

for i = 1:10
    x_idx = (i-1)*2 + 1;
    x_cords = data(:,x_idx);
    y_cords = data(:,x_idx+1);
    T = 1.000 - (i-1)*0.1;

    figure;
    plot(x_cords, y_cords,'*');
    xlabel("X", "FontSize", 18);
    ylabel("Y", "FontSize", 18);
end