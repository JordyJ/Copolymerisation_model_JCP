#!/bin/sh
maindirectory=$(pwd inputlist.csv)

echo "Making input files from inputlist.csv"
#### Generate the input files from $1
awk -F, -vOFS="," 'NR>1{ 
print "seed,length_limit,time_limit,max_pols,k0,k,kact,ConcEff,Gact,Gbb,Conc0,Conc1,G0,G1,Ggen,Gend,outputfilename,outputtrajfilename,print_flag,print_valid_trans_flag,write_output_flag,write_traj_flag,traj_write_interval,breakable_backbone_flag,off_rate_discrim_flag,t_activation_flag,local_updates_flag,use_test_rates_flag" >> "input_"(NR-1)".csv";
print $0 >> "input_"(NR-1)".csv"
}' inputlist.csv

END="$(wc -l < inputlist.csv)"

echo "$END input files generated"

echo Running sim over $END input files

for ((i=1;i<$END;i++)); do

mkdir "$maindirectory/$i/"    

mv "input_${i}.csv" $maindirectory/${i}/input.csv

cd $i
echo  
echo Running input $i
echo  
../sim input.csv
echo  
echo Analysing and plotting output.
echo  

python ../plot_output.py output.csv

cd ../
done

