echo "Simulating w/ T = 1.0";
./sim N=200 rho=0.25 nblock=10 T=1.0 run > log/ex6_2_mc_T_1.0.txt;

for r in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1;
do
o=$(python -c "print(\"{:.3f}\").format($r+0.1)");
echo "Simulating w/ T = ${r}, using conf from where T = ${o}";
./sim N=200 rho=0.25 nblock=10 T=${r} read=0200_r0.250_T${o} run > log/ex6_2_mc_T${r}.txt;
done