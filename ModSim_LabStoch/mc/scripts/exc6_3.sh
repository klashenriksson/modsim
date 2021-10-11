for r in 0.15 0.20 0.25 0.30 0.40 0.50 0.60 0.70 0.80;
do
    echo "Simulating.. rho=${r}, T=1.0..";
    ./sim N=64 rho=${r} ntherm=10000 nblock=100 T=1.0 run > log/ec6_3_mc_r${r}_T1.txt
done
echo "Done!";

for r in 0.15 0.20 0.25 0.30 0.40 0.50 0.60 0.70 0.80;
do
    for t in 0.9 0.8 0.7 0.6 0.5;
    do
        echo "Simulating.. rho=${r}, T=${t}..";
        o=$(python -c "print(\"{:.3f}\").format($t+0.1)");
        ./sim N=64 rho=${r} ntherm=10000 nblock=100 T=${o} read T=${t} run > log/ec6_3_mc_r${r}_T${t}.txt;
    done
done

echo "All done!";