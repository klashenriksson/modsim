for r in 0.01 0.1 1.0; do
    echo "Simulating w/ alpha=${r}";
    ./sim N=64 rho=0.2 T=1.0 alpha=${r} > log/exc6_1_alpha${r}.txt run &
done
