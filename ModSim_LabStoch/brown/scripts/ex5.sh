for r in 0.001 0.002 0.003 0.004; do
    echo "Simulating w/ deltat=${r}";
    ./sim N=64 rho=0.3 T=1.0 alpha=10 deltat=$r > log/exc5_deltat${r}.txt run &
done
