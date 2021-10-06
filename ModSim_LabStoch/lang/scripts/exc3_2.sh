for r in 0.01 0.1 1.0;
do 
    echo "Simulating w/ alpha=${r}";
    ./sim N=64 rho=0.6 nblock=10 T=1 deltat=0.01 alpha=$r run > log/exc3_2_alpha_${r}.txt & 
done
echo "Done";
