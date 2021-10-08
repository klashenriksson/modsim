import os

data = [];
path = "log/exc6_3";

for file in os.listdir(path):
    curr = os.path.join(path, file);
    if os.path.isfile(curr):
        rho = -1;
        T = -1;
        p = -1;
        with open(curr) as f:
            lines = f.readlines();

            for line in lines:
                T_idx = line.find("T = ");
                if T_idx >= 0:
                    T_str = line[T_idx+4:line.find(",",T_idx)];
                    T = float(T_str);

                rho_idx = line.find("rho = ");
                if rho_idx >= 0:
                    rho_str = line[rho_idx+6:line.find(",",rho_idx)];
                    rho = float(rho_str);

                pressure_idx = line.find("Pressure   : ");
                if pressure_idx >= 0:
                    pressure_str = line[rho_idx+13:line.find(" +",pressure_idx)];
                    p = float(pressure_str);
        if rho == -1:
            print("Failed to find rho");
        if T == -1:
            print("Failed to find T");
        if p == -1:
            print("Failed to find pressure");

        if rho == -1 or T == -1 or p == -1:
            exit();
        data.append((rho, T, p));
        print("Read rho={rho}, T={t}, pressure={p}".format(rho=rho, t=T, p=p));

sdata = sorted(data, reverse=True);
with open("log/exc6_3_data.txt", "w") as ofile:
    ofile.write("# rho, T, pressure\n");
    for datapoint in sdata:
        rho, t, p = datapoint;
        ofile.write("{rho} {t} {p}\n".format(rho=rho, t=t, p=p));