#!/usr/bin/bash


########################################################################################################################
## 版本: 1.1.1
## 作者: 李伟宁 liwn@cau.edu.cn
## 日期: 2023-07-05
## 简介: 用于基于BLUP/Bayes交叉验证准确性计算
##
## 使用: ./single_breed.sh --label "breedA" --phef /path/to/your/phenotype ...(详细参数请通过--help查看)
## 依赖软件/环境: 
##  1. R
##  2. plink/1.9
##  3. gmatrix
##  4. mbBayesAB
##  5. 其他R语言和Bash脚本
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


###################  参数处理  #####################
####################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## 参数名
TEMP=$(getopt -o 4hp:b:v:m:o: --long rmNeg,intercept,debug,label:,phef:,bfile:,varf:,pedf:,gmat:,DIR:,rep:,fold:,add_rf:,gen:,year:,iyse:,tbvf:,tbv_col:,phereal:,binf:,seed:,all_eff:,iter:,burnin:,bin:,thin:,report:,ran_eff:,method:,thread:,miss:,gidf:,add_sol:,invA:,gebvf:,num_int:,code:,alpha:,out:,overWri,valphe,dmu4,dense,help \
  -n 'javawrap' -- "$@")
if [ $? != 0 ]; then
  echo "Terminating..." >&2
  exit 1
fi
eval set -- "$TEMP"
## 解析参数
while true; do
  case "$1" in
  -p |--phef)    phef="$2";       shift 2 ;; ## 表型文件
  -b |--bfile)   bfile="$2";      shift 2 ;; ## plink二进制文件前缀
      --label)   label="$2";      shift 2 ;; ## 品种标志符
  -v |--varf)    varf="$2";       shift 2 ;; ## 方差组分文件(如没提供则用系谱/SNP信息估计)
      --pedf)    pedf="$2";       shift 2 ;; ## 系谱文件(不提供则运行GBLUP)
  -m |--gmat)    gmat="$2";       shift 2 ;; ## 用户提供关系矩阵或逆矩阵(id id value)
      --gidf)    gidf="$2";       shift 2 ;; ## 基因型个体id，与用户指定的G阵文件中个体id一致
      --DIR)     DIR="$2";        shift 2 ;; ## 参数卡文件前缀 [within]
      --rep)     rep="$2";        shift 2 ;; ## 交叉验证重复数 [1]
      --fold)    fold="$2";       shift 2 ;; ## 交叉验证折数 [1]
      --add_rf)  add_rf="$2";     shift 2 ;; ## 加性效应所在分组，默认在第1个随机效应组 [1]
      --gen)     gen="$2";        shift 2 ;; ## 有表型个体中倒数gen个世代个体为验证群 [ALL]
      --year)    year="$2";       shift 2 ;; ## 有表型个体中在该年份之后出生为验证群
      --iyse)    iyse="$2";       shift 2 ;; ## 提供year时须提供，格式为idcol:yearcol:yeas_star:year_end
      --tbvf)    tbvf="$2";       shift 2 ;; ## 真实育种值文件，若提供tbv_col但tbvf缺失，则tbvf设为表型文件
      --tbv_col) tbv_col="$2";    shift 2 ;; ## 真实育种值在表型文件的第几列，当真实育种值即为表型文件中相应列表型时参数应为"same"
      --phereal) phereal="$2";    shift 2 ;; ## 表型在表型文件中实数列的位置 [1]
      --all_eff) all_eff="$2";    shift 2 ;; ## DIR中$MODEL第3行，前3位数不需要，只需要所有效应所在的列，如"2 3 1"
      --ran_eff) ran_eff="$2";    shift 2 ;; ## DIR中$MODEL第4行，第1位数不需要，只需要所有随机效应所在分组，如"1" [1]
      --method)  method="$2";     shift 2 ;; ## 育种值估计方法，可为PBLUP/GBLUP/ssGBLUP/BayesAS [GBLUP]
      --bin )    bin="$2";        shift 2 ;; ## 是否合并临近窗口，fix/frq/ld/ind [ind]
      --binf )   binf="$2";       shift 2 ;; ## 区间文件
      --seed )   seed="$2";       shift 2 ;; ## MCMC抽样及验证群划分时的随机种子 [40296]
      --dense )  dense=true;      shift   ;; ## 将DMU的ANALYSE中的method设为31，利用多线程计算方差组分
      --iter )   iter="$2";       shift 2 ;; ## MCMC 总循环次数
      --burnin ) burnin="$2";     shift 2 ;; ## MCMC burn-in循环数
      --thin )   thin="$2";       shift 2 ;; ## MCMC抽样步长
      --report ) report_sep="$2"; shift 2 ;; ## MCMC报告间隔
      --thread)  thread="$2";     shift 2 ;; ## dmu任务并行数
      --miss)    miss="$2";       shift 2 ;; ## 缺失表型表示符 [-99]
      --invA)    invA="$2";       shift 2 ;; ## A逆构建方式(1/2/3/4/6)，1考虑近交，2不考虑近交，其他见DMU说明书 [1]
      --gebvf )  gebvf="$2";      shift 2 ;; ## Bayes模型中gebvf输出文件夹的名前缀 [EBV_${bin}]
      --num_int) num_int="$2";    shift 2 ;; ## 整型列列数 [根据文件计算]
      --code)    code="$2";       shift 2 ;; ## 代码路径 [NULL]
      --alpha)   alpha="$2";      shift 2 ;; ## G阵校正系数(是否考虑近交) [0.05]
      --intercept) mean=true;     shift   ;; ## 在整数列最后一列后增加一列群体均值
      --rmNeg)    rmNeg=true;     shift   ;; ## 在整数列最后一列后增加一列群体均值
      --valphe)   valphe=true;    shift   ;; ## 验证群中个体只需要有表型，不需要含有基因型
      --overWri)  overwrite=true; shift   ;; ## 若已经存在将要生成的文件，进行覆盖
      --debug)    debug=true;     shift   ;; ## 不跑DMU和gmatrix
  -o | --out)     out="$2";       shift 2 ;; ## 输出校正表型文件名
  -4 | --dmu4)    dmu4=true;      shift   ;; ## 用DMU4模型估计育种值
  -h | --help)    grep ";; ##" $0 | grep -v help && exit 1 ;;
  -- ) shift; break ;;
  * ) break ;;
  esac
done

## 避免执行R脚本时的警告("ignoring environment value of R_HOME")
unset R_HOME

## 默认参数
gebvf=${gebvf:=EBV_${bin}}
seed=${seed:=40296}                           ## 随机数种子
nchr=${nchr:=30}                              ## 染色体数目
DIR=${DIR:=within}                            ## 参数卡名称
iter=${iter:=30000}                           ## MCMC 总循环次数
burnin=${burnin:=20000}                       ## MCMC burn-in循环数
thin=${thin:=10}                              ## MCMC抽样步长
report_sep=${report_sep:=100}                 ## MCMC报告间隔
alpha=${alpha:=0.05}                          ## G阵参数
thread=${thread:=1}                           ## dmu任务并行数
phereal=${phereal:=1}                         ## 表型列
add_rf=${add_rf:=1}                           ## 加性随机效应所在组
ran_eff=${ran_eff:=1}                         ## 随机效应
rep=${rep:=1}                                 ## 重复
fold=${fold:=1}                               ## 验证分组数
miss=${miss:=-99}                             ## 缺失表型表示
method=${method:=GBLUP}                       ## 育种值估计方法
out=${out:=accuracy.txt}                      ## 准确性输出文件名
invA=${invA:=2}                               ## A逆构建方式，默认考虑近交
[[ ${tbv_col} && ! ${tbvf} ]] && tbvf=${phef} ## 含有真实育种值的文件
[[ ${tbvf} && ! ${tbv_col} ]] && tbv_col=2    ## 真实育种值文件存在，默认育种值在文件第二列，id在第一列

## 主文件夹
workdir=$(pwd)

## 日志文件夹
logp=${workdir}/log
mkdir -p ${logp}

## 脚本所在文件夹
if [[ ${code} ]]; then
  [[ ! -d ${code} ]] && echo "${code} not exists! " && exit 5
else
  script_path=$(dirname "$(readlink -f "$0")")
  code=$(dirname "$script_path")
fi

## 脚本
phe_group=${code}/R/validation_population_define.R
accur_cal=${code}/R/accuracy_bias_calculation.R
keep_phe_gid=${code}/R/keep_pheno_geno_individuals.R
job_pool=${code}/shell/parallel_jobs_func.sh
func=${code}/shell/function.sh

## 载入自定义函数
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${job_pool}
source ${func}

## 将程序路径加到环境变量中
export PATH=${code}/bin:$PATH

## 检查必要参数是否提供
if [[ ! -s ${phef} ]]; then
  echo "phenotype file ${phef} not found! "
  exit 1
elif [[ ! -s ${bfile}.fam && ! -s ${bfile}.map && ! -s ${gidf} && ! -s ${pedf} ]]; then
  echo "plink file ${bfile}.fam, pedigree file ${pedf} or genotyped individuals id file ${gidf} not found! "
  exit 1
elif [[ -s ${gmat} && ${method} == "ssGBLUP" && ! -s ${gidf} ]]; then
  echo "genotyped individuals id file ${gidf} not found! "
  exit 1
elif [[ ! ${label} ]]; then
  echo "breed label parameter --label not found! "
  exit 1
fi

## 检查需要的程序是否在环境变量中能检索到并且可执行
check_command plink gmatrix mbBayesAB LD_mean_r2 run_dmu4 run_dmuai

## 检查需要的脚本文件是否存在且具有执行权限
check_command $phe_group $accur_cal $keep_phe_gid $job_pool $func


##################  解析命令行参数  ##############
################################################

## DMU多线程
if [[ ${method} == "GBLUP" || ${dense} ]]; then
  method_dmu="31"
else
  method_dmu="1"
fi

# ## 作业线程
# [[ ${SLURM_CPUS_ON_NODE} ]] && thread=${SLURM_CPUS_ON_NODE}

## 随机数种子
if [[ ! ${seed} ]]; then
  seed=$RANDOM
  echo "$seed" >MCMC_seed.txt
fi

## 判断是否已完成验证群划分等步骤
val_phe="${workdir}/val1/rep1/pheno.txt"
if [[ ! -s ${val_phe} || ${overwrite} ]]; then
  overwrite=true
else
  unset overwrite
fi

## 线程池初始化
job_pool_init ${thread} 0

## 复制一份表型文件到工作文件夹
cp ${phef} ${workdir}/pheno_within.txt
phef=${workdir}/pheno_within.txt

## 表型文件整型、实型变量列数
if [[ ! ${num_int} ]]; then
  ncol=$(awk 'END{print NF}' ${phef})
  for i in $(seq 1 ${ncol}); do
    dot=$(awk -vl=${i} '{print $l}' ${phef} | grep -c "\.")
    [[ ${dot} -gt 0 ]] && num_int=$((i - 1)) && break
  done
fi
num_real=$(($(awk 'END{print NF}' ${phef}) - num_int))

## 效应个数
nA=$(echo ${all_eff} | awk '{print NF}')
nR=$(echo ${ran_eff} | awk '{print NF}')

###################  基因型文件  ####################
#####################################################
if [[ ${bfile} && ! -s ${bfile}.bim ]]; then
  if [[ ! -s ${bfile}.map ]]; then
    echo "plink file ${bfile}.map not found! "
    exit 1
  else
    ## 提取指定染色体上的标记
    plink \
      --file ${bfile} \
      --chr-set ${nchr} \
      --make-bed --out ${label} >${logp}/plink_single_copy.log
    bfile=$(pwd)/${label}
  fi
elif [[ -s ${bfile}.bim ]]; then
  plink \
    --bfile ${bfile} \
    --chr-set ${nchr} \
    --make-bed --out ${label} >${logp}/plink_single_copy.log
  bfile=$(pwd)/${label}
fi

#####################  表型文件  ####################
####################################################
phe_col=$((phereal + num_int))
## 筛选出同时有基因型和表型和个体作为参考群
if [[ -s ${bfile}.fam ]]; then
  $keep_phe_gid \
    --famf "${bfile}.fam" \
    --phef ${phef} \
    --rm single \
    --num_int ${num_int} \
    --phec ${phe_col} \
    --rmOut gid_miss_phe.txt
fi
## 删除没有表型信息的基因型个体
if [[ -s gid_miss_phe.txt ]]; then
  echo "remove $(wc -l <gid_miss_phe.txt) individuals without phenotype"
  plink \
    --bfile ${bfile} \
    --chr-set ${nchr} \
    --remove gid_miss_phe.txt \
    --make-bed --out ${bfile} >${logp}/plink_single_rm_phemiss.log
fi
## 表型文件中整型列后添加1列截距项(整列设为"1")
if [[ ${mean} ]]; then
  ## 参数卡信息
  ((nA++))
  ((num_int++))
  all_eff="${num_int} ${all_eff}"

  ## 表型文件
  awk -v column="${num_int}" -v value="1" '
    BEGIN {
        FS = OFS = " ";
    }
    {
        for ( i = NF + 1; i > column; i-- ) {
            $i = $(i-1);
        }
        $i = value;
        print $0;
    }
    ' ${phef} >${phef}.tmp
  echo "add populations mean column (${num_int}) in the phenotype file."
  mv ${phef}.tmp ${phef}
fi
## 检查表型文件中是否还有个体
if [[ ! -s ${phef} ]]; then
  echo "no individuals in phenotype file ${phef}! "
  [[ ${workdir} =~ "mbGS" ]] && rm -r ${workdir}
  exit 1
fi

###################  dmu参数卡模板  ####################
########################################################
if [[ ${method} =~ "BLUP" ]]; then
  ## ANALYSE
  if [[ ${dmu4} ]]; then
    ANALYSE="11 9 0 0"
  else
    ANALYSE="1 ${method_dmu} 0 0"
  fi

  [[ -s ${DIR}.DIR ]] && echo "warn: ${DIR}.DIR will be overwrited! "
  {
    echo "\$COMMENT"
    echo "get EBV of individuals in validation population with reduced phenotypes"
    echo "\$ANALYSE ${ANALYSE}"
    echo "\$DATA  ASCII (${num_int}, ${num_real}, ${miss}) pheno.txt"
    echo -e "\$MODEL\n1\n0\n${phereal} 0 ${nA} ${all_eff}\n${nR} ${ran_eff}\n0\n0"
    echo '$VAR_STR %VAR_STR%'
    echo '$PRIOR %PRIOR%'
    echo '$SOLUTION'
  } >${DIR}.DIR
fi

#####################  方差结构  #####################
######################################################
if [[ ${method} == "GBLUP" ]]; then
  ## 方差协方差结构文件
  if [[ ${gmat} ]]; then
    if [[ ! -s ${gmat} ]]; then
      echo "${gmat} not found! "
      exit 1
    else
      echo "Use the user specified genome relationship matrix: ${gmat}"
    fi
  else
    echo "Read the plink bed file and Calculate the additive G matrix..."
    [[ ! ${debug} ]] && gmatrix --bfile ${bfile} --grm agrm --out full --inv
    echo "G matrix created."
    gmat=${workdir}/full.agiv.id_fmt
  fi
  ## 参数卡
  sed -i "s#%VAR_STR%#${add_rf} GREL ASCII ${gmat}#g" ${DIR}.DIR
  ## 加性效应在SOL结果文件中的代码
  add_sol=3
elif [[ ${method} == "PBLUP" ]]; then
  ## PBLUP
  sed -i "s#%VAR_STR%#${add_rf} PED ${invA} ASCII ${pedf}#g" ${DIR}.DIR
  add_sol=4
elif [[ ${method} == "ssGBLUP" ]]; then
  if [[ ${gmat} ]]; then
    [[ ! -s ${gmat} ]] && echo "${gmat} not found! " && exit 1
    [[ ! -s ${gidf} ]] && echo "${gidf} not found! " && exit 1
  else
    echo "Read the plink bed file and Calculate the additive G matrix..."
    [[ ! ${debug} ]] && gmatrix --bfile ${bfile} --grm agrm --out full
    echo "G matrix created."
    gmat=${workdir}/full.agrm.id_fmt
    gidf=${workdir}/full.id
  fi
  sed -i "s#%VAR_STR%#${add_rf} PGMIX ${invA} ASCII ${pedf} ${gidf} ${gmat} ${alpha} G-ADJUST#g" ${DIR}.DIR
  add_sol=4
fi

## 指定方差组分(初始值)
if [[ ${method} =~ "BLUP" ]]; then
  if [[ -s ${varf} ]]; then
    sed -i "s#%PRIOR%#${varf}#g" ${DIR}.DIR
  else
    if [[ ${dmu4} ]]; then
      echo "${varf} not found! "
      exit 1
    else
      [[ ${varf} ]] && echo "${varf} not found! "
      sed -i '/$PRIOR.*/d' ${DIR}.DIR
    fi
  fi
fi

###################  验证群划分  #####################
######################################################
## 参数设置
if [[ ${overwrite} ]]; then
  option=' '
  [[ ${gen} ]] && option="${option} --gen ${gen} --pedf ${pedf}"               ## 按世代划分
  [[ ${year} && ${iyse} ]] && option="${option} --year ${year} --iyse ${iyse}" ## 按出生年份划分
  [[ ! ${valphe} && -s ${bfile}.fam ]] && option="${option} --fam ${bfile}.fam"

  ## 划分出参考群和验证群
  $phe_group \
    --phef "${phef}" \
    --nonmiss "${all_eff}" \
    --rep ${rep} \
    --fold ${fold} \
    --seed ${seed} \
    --outvid "val.id" \
    --outdir "${workdir}/val#val#/rep#rep#" \
    --rminvail \
    --keepid "${workdir}/keep_fid_id.txt" \
    --pheCol ${phe_col} \
    ${option}

  ## 只保留可以作为参考群的个体基因型(有表型、基因型、固定效应水平不缺失)
  if [[ -s ${workdir}/keep_fid_id.txt ]]; then
    plink --bfile ${bfile} \
      --chr-set ${nchr} \
      --keep ${workdir}/keep_fid_id.txt \
      --make-bed --out ${bfile} >${logp}/plink_full_set.log

    rm ${bfile}.fam~ ${bfile}.bim~ ${bfile}.bed~
    echo "keep $(wc -l <${workdir}/keep_fid_id.txt) individuals in geneotype file."
  fi
fi

###################  验证群育种值估计  #####################
############################################################
for r in $(seq 1 ${rep}); do # r=1;f=1
  for f in $(seq 1 ${fold}); do
    ## 替换缺失值
    sed -i "s/NA/${miss}/Ig" ${workdir}/val${f}/rep${r}/pheno.txt
    check_alphabet ${workdir}/val${f}/rep${r}/pheno.txt

    ## 育种值估计
    if [[ ${method} == "BayesAS" ]]; then
      ## 如果没有提供binf文件，则执行BayesA
      if [[ ! -s ${binf} ]]; then
        binf=${workdir}/val${f}/rep${r}/fixed_bins_1.txt
        awk '{print "1"}' ${bfile}.bim >${binf}
        echo "Fit bayesA with only one snp in each bins"
      else
        ## 划分类型
        case ${binf} in
          *fix*) bin="fix" ;;
          *ld*) bin="ld" ;;
          *frq*) bin="frq" ;;
          *) bin="unknown" ;;
        esac
      fi

      ## 固定效应
      fix_eff=${all_eff%" ${ran_eff}"}

      ## 运行Bayes模型
      [[ $(find ${workdir}/val${f}/rep${r} -name "${gebvf}*.txt" 2>/dev/null | wc -l) -ge 1 ]] && continue # temp
      job_pool_run mbBayesAB \
        --bfile ${bfile} \
        --phef ${workdir}/val${f}/rep${r}/pheno.txt \
        --fix "${fix_eff}" \
        --phe ${phe_col} \
        --binf ${binf} \
        --iter ${iter} \
        --burnin ${burnin} \
        --thin ${thin} \
        --outdir ${workdir}/val${f}/rep${r} \
        --report ${report_sep} \
        --varOut var_${bin}.txt \
        --effOut effect_${bin}.txt \
        --gebvOut ${gebvf} \
        --mcmcOut MCMC_process_${bin}.txt \
        --seed ${seed} \
        --logf ${bin}_gibs_${SLURM_JOB_ID}.log
    fi

    if [[ ${method}  =~ "BLUP" ]]; then
      ## 参数卡
      cp ${workdir}/${DIR}.DIR ${workdir}/val${f}/rep${r}/${DIR}.DIR

      ## 运行BLUP模型
      [[ $(find ${workdir}/val${f}/rep${r} -name "*.SOL" 2>/dev/null | wc -l) -ge 1 ]] && continue # temp
      if [[ ${dmu4} ]]; then
        [[ ! ${debug} ]] && job_pool_run run_dmu4 ${DIR} ${workdir}/val${f}/rep${r}
      else
        [[ ! ${debug} ]] && job_pool_run run_dmuai ${DIR} ${workdir}/val${f}/rep${r}
      fi
    fi
  done
done

## 等待后台进程结束
job_pool_wait
## 释放线程
job_pool_shutdown

####################  准确性计算  ##################
###################################################
## 参数
[[ ${rmNeg} ]] && rmNeg=" --rmNeg " ## 是否删除准确性计算结果中的负值（异常）
if [[ ${tbv_col} == "same" ]]; then
  ## 文件中的原表型列作为校正表型
  option=" --tbv_col ${phe_col}"
elif [[ ${tbv_col} ]]; then
  option=" --tbv_col ${tbv_col}"
fi
[[ ${tbvf} ]] && option="${option} --tbvf ${tbvf} "
if [[ ${method} =~ 'BLUP' ]]; then
  ## BLUP模型
  option="${option} --add_sol ${add_sol} --dir_val ${workdir}/val#val#/rep#rep#/${DIR}"
  ebv_col=1
else
  ## BayesAS模型
  ebv_col=2
  option="${option} --ebvf ${workdir}/val#val#/rep#rep#/${gebvf}_y1.txt"
fi
## 计算准确性
$accur_cal \
  --ebv_col ${ebv_col} \
  --val_idf ${workdir}/val#val#/rep#rep#/val.id \
  --rep ${rep} \
  --fold ${fold} \
  ${rmNeg} \
  ${option} \
  --out ${workdir}/${out}

###################  删除中间文件  #####################
########################################################
## 参数卡
[[ -s ${workdir}/${DIR}.DIR ]] && rm ${workdir}/${DIR}.DIR

[[ $? -ne 0 ]] && echo "Accuracy calculation completed, file output to: ${out}"
