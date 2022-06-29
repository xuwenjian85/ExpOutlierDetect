
# before use OUTRIDER to control cofounders, we need to know the best encoding dimension q0 
# python处理数据集，输出配置文件：q值的范围，文件夹、文件名之类的
date
find_enc_dim=/media/eys/xwj/expression_outlier/find_enc_dim.R
inj_mask=/media/eys/xwj/expression_outlier/inject_mask_w_seed.R
inj_outlier=/media/eys/xwj/expression_outlier/inject_zscore_outlier_w_seed.R
run_OUTRDIER=/media/eys/xwj/expression_outlier/run_outrider.R

id=1
dname=blood # gtex_skin and blood
N_sample=601
N_simu=100
ZScore=2,3
cts='/media/eys/xwj/RNAseq/public_normal/df_cts_HC1157_corrupt_fc2_ngene100_nrep3.txt'
outdir=/media/eys/xwj/expression_outlier/output

mkdir -p $outdir/logs
mkdir -p $outdir/${dname}

#######--------1.find the best encode dimensions q0------------#####

# manually set 10 values around (1/5)*N_sample. 
# if N_sample = 1000, then q_array = [..., 180, 190, 200, 210, 220,...]
Ns=$((N_sample/5/10*10)) ## N_sample rounded to the nearest tens
q_array=$(seq $((Ns-50)) 10 $((Ns+50)))
echo $id,$dname,$N_sample,$N_simu, $ZScore, $outdir
echo $Ns, $q_array

source /public/home/test1/soft/anaconda3/etc/profile.d/conda.sh
conda activate R4.1_OUTRIDER
# loop through a array of encode dimensions q with find_enc_dim.R
i=0 # I assign jobs to many cpus
# for q in $q_array; 
# do 
#     echo "$dname: encdim $q -> cpu $i"
#     out1=$outdir/${dname}/find_enc_dim_$q.encDimqTable
#     out2=$outdir/${dname}/find_enc_dim_$q.RDS
#     log=$outdir/logs/${dname}_find_enc_dim_$q.log
#     # nohup command > t.log 2>&1 &
#     taskset -c $i Rscript $find_enc_dim $dname $cts $q $out1 $out2 > $log 2>&1 &
#     i=$((i+1))
# done
wait

q0=$(sort -nrk 3 $outdir/${dname}/find_enc_dim_*Table | head -n1 | awk '{print $2}')
echo $q0
date
### TODO plot q0 on q curve in python

#######--------2.simulate datasets------------#####
#  by introducing random corruptions: 
# "inj_mask" random mask (gene,sample) in matrix using the string i as random seed
# "inj_outlier" zscore * mask to simulate count matrix
Freq=1e-4
for j in $(seq -w 1 $N_simu); do
    
    in_ods=$outdir/$dname/find_enc_dim_$q0.RDS
    out_ods_true_mask=$outdir/$dname/ods_true_mask_$j.RDS
    # Rscript $inj_mask $in_ods $j $Freq $out_ods_true_mask
    
    for Z in ${ZScore//,/ }; do
        out_ods_simu=$outdir/$dname/ods_simu_$Z"_"$j.RDS
        echo $j,$Z, $out_ods_simu
        # Rscript $inj_outlier $out_ods_true_mask $Z $out_ods_simu
    done
done

#######--------3.runOUTRDIER first & control confounders-----------#####
# 
i=0
for Z in ${ZScore//,/ }; do
    for j in $(seq -f "%03g" 0 10); do
        ods_simu=$outdir/$dname/ods_simu_$Z"_"$j.RDS
        ods_simu_outrider=$outdir/$dname/ods_simu_outrider_$Z"_"$j.RDS
        log=$outdir/logs/${dname}_run_OUTRDIER_q$q0"_z"$Z"_"$j.log
        #taskset xx Rscript run_OUTRDIER.R simu.RDS q0 OUTRIDER_res.RDS
        echo $j,$q0, $ods_simu,'->', $i
        taskset -c $i Rscript $run_OUTRDIER $ods_simu $q0 $ods_simu_outrider > $log 2>&1 &
        i=$((i+1))
    done
    wait
done
exit 0


for i in N_simu:
    
    
# convert RDS to txt
for i in N_simu:
    Rscript rds2tsv.R OUTRIDER_res.RDS cts_simu.txt OUTRIDER_res.txt df_corrupt.txt

# python 
makefeature--> X_y.pkl (for lopp N_simu)
model(lengfei)
compare_models --> PRcurve + 95%CI
