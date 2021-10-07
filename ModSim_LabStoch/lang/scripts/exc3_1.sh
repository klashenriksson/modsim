for r in 0.002 0.005 0.01 0.014 0.020;
do 
    echo "Simulating w/ deltat=${r}";
    ./sim N=64 rho=0.5 nblock=50 T=1 deltat=$r read=0064_r0.500_T1.000_start run > log/exc3_1_deltat_${r}.txt & done
    echo "Done";
