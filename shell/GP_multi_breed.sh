#!/usr/bin/bash


########################################################################################################################
## 版本: 1.1.1
## 作者: 李伟宁 liwn@cau.edu.cn
## 日期: 2023-07-05
## 简介: 用于估计联合评估BLUP/Bayes的预测准确性，需要先完成群体内准确性计算(即不同群体目录下有val*/rep*文件夹)再运行此脚本
##
## 使用: ./multi_breed.sh --pops "breedA breedB" --bfileM /path/to/your/plinkFile ...(详细参数请通过--help查看)
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

###################  参数处理  ###################
#################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
TEMP=$(getopt -o h4 --long ref:,rg_local,keep_all,nbin:,binf:,summs:,dirPre:,noCov,debug,traceplot,bfileW:,bfileM:,pops:,pedf:,VaSori:,R0Sori:,VaSSori:,R0SSori:,priorCor:,priSaDmu:,VaSd:,R0Sd:,DIR:,vidf:,method:,rg:,re:,type:,GmatM:,gmat:,gidf:,invA:,all_eff:,ran_eff:,tbvf:,tbv_col:,phereal:,add_rf:,fold:,rep:,miss:,bincol:,num_int:,win:,nsnp_win:,phyNum:,r2_merge:,bin:,LD_maxD:,r2:,inter:,fix_snp:,seed:,iter:,burnin:,thin:,report_sep:,code:,gebvf:,thread:,alpha:,out:,nchr:,updatepri:,prefix:,dfva:,dfve:,varSamp,VarA,dmu4,overwri,suffix,dense \
  -n 'javawrap' -- "$@")
if [ $? != 0 ]; then
  echo "Terminating..." >&2
  exit 1
fi
eval set -- "$TEMP"

## 解析参数
while true; do
  case "$1" in
  --pops )      pops="$2";       shift 2 ;; ## 群体A标识符，如'YY DD'
  --bfileW )    bfileW="$2";     shift 2 ;; ## 每个品种的基因型文件，如'./data/YY ./data/DD'
  --bfileM )    bfileM="$2";     shift 2 ;; ## plink二进制文件前缀，群体合并后的
  --binf )      binf="$2";       shift 2 ;; ## 区间文件
  --bincol )    bincol="$2";     shift 2 ;; ## 区间文件中指示区间内标记数的列 [1]
  --pedf )      pedf="$2";       shift 2 ;; ## 系谱文件，可为合并的单个系谱，也可为以空格分隔的每个群体的系谱文件
  --summs )     summs="$2";      shift 2 ;; ## GEMMA结果文件前缀
  --priorCor )  priorCor="$2";   shift 2 ;; ## 相关系数先验值文件 [NULL]
  --VaSori )    VaSori="$2";     shift 2 ;; ## 方差组分先验信息来源，phe/iden/FILE [NULL]
  --R0Sori )    R0Sori="$2";     shift 2 ;; ## 方差组分先验信息来源，phe/iden/FILE [NULL]
  --VaSSori )   VaSSori="$2";    shift 2 ;; ## 方差组分先验分布的尺度矩阵信息来源，phe/iden [NULL]
  --R0SSori )   R0SSori="$2";    shift 2 ;; ## 方差组分先验分布的尺度矩阵信息来源，phe/iden [NULL]
  --priSaDmu )  priSaDmu="$2";   shift 2 ;; ## 方差组分先验信息来源，phe/iden/FILE [NULL]
  --VaSd )      VaSd="$2";       shift 2 ;; ## 加性效应先验分布尺度参数服从的分布 NULL/wishart/invgamma [NULL]
  --R0Sd )      R0Sd="$2";       shift 2 ;; ## 残差效应先验分布尺度参数服从的分布 NULL/wishart [NULL]  
  --h2 )        h2="$2";         shift 2 ;; ## 遗传力，用于生成初始方差组分，两个性状不同时空格隔开，如 '0.3 0.1' [0.5]
  --DIR )       DIR="$2";        shift 2 ;; ## 参数卡前缀 [type]
  --method )    method="$2";     shift 2 ;; ## dmu评估模型，PBLUP/GBLUP/ssGBLUP
  --type )      type="$2";       shift 2 ;; ## 参考群合并方式 blend/union/multi(两性状模型) [blend]
  --GmatM )     GmatM="$2";      shift 2 ;; ## 基因型矩阵构建方式[multi/single]
  --gmat )      gmat="$2";       shift 2 ;; ## 用户提供关系矩阵或逆矩阵(id id value)
  --gidf )      gidf="$2";       shift 2 ;; ## 基因型个体id，与用户指定的G阵文件中个体id一致
  --invA )      invA="$2";       shift 2 ;; ## A逆构建方式(1/2/3/4/6)，1为考虑近交，2为不考虑近交，其他见DMU说明书 [1]
  --all_eff  )  all_eff="$2";    shift 2 ;; ## DIR中$MODEL第3行(所有效应)，前3位不用，如"2 1"，品种间不同时用-分隔
  --ran_eff  )  ran_eff="$2";    shift 2 ;; ## DIR中$MODEL第4行(随机效应分)，第1位不用，如"1"，品种间不同时用-分隔
  --tbvf )      tbvf="$2";       shift 2 ;; ## 包含真实育种值的文件
  --tbv_col )   tbv_col="$2";    shift 2 ;; ## 真实育种值在表型文件的第几列，当真实育种值即为表型文件中相应列表型时参数应为"same"
  --phereal )   phereal="$2";    shift 2 ;; ## 表型在表型文件中实数列的位置
  --add_rf )    add_rf="$2";     shift 2 ;; ## 加性效应所在分组
  --fold )      fold="$2";       shift 2 ;; ## 交叉验证倍数
  --rep )       rep="$2";        shift 2 ;; ## 重复计算次数
  --miss )      miss="$2";       shift 2 ;; ## 缺失表型表示符
  --num_int )   num_int="$2";    shift 2 ;; ## 表型文件中整型变量的列数 [从A品种表型文件中获取]
  --nsnp_win )  nsnp_win="$2";   shift 2 ;; ## 合并前每个窗口SNP数
  --phyNum )    phyNum="$2";     shift 2 ;; ## 以固定SNP数还是物理位置position划分区间，num/pos [num]
  --win )       win="$2";        shift 2 ;; ## 计算指定SNP侧翼R2均值时的窗口大小 [50]
  --bin )       bin="$2";        shift 2 ;; ## 是否合并临近窗口，fix/frq/ld [fix]
  --nbin )      nbin="$2";       shift 2 ;; ## 当bin为fix时，划分区间数约束为与nbin相近
  --ref )       ref="$2";        shift 2 ;; ## Reference panel for dividing blocks M/1[index of breed]/2/... [M]
  --r2_merge )  r2_merge="$2";   shift 2 ;; ## 合并临近窗口时的LD阈值
  --LD_maxD )   LD_maxD="$2";    shift 2 ;; ## plink参数，计算LD时，SNP对的最大距离kb
  --r2 )        r2="$2";         shift 2 ;; ## plink参数，计算LD时，输出结果文件中R2的阈值
  --inter )     inter="$2";      shift 2 ;; ## plink参数，计算LD时，间隔位点数多于inter的标记对不计算LD值
  --seed )      seed="$2";       shift 2 ;; ## MCMC抽样时的随机种子
  --iter )      iter="$2";       shift 2 ;; ## MCMC 总循环次数
  --burnin )    burnin="$2";     shift 2 ;; ## MCMC burn-in循环数
  --thin )      thin="$2";       shift 2 ;; ## MCMC抽样步长
  --report )    report_sep="$2"; shift 2 ;; ## MCMC报告间隔
  --code )      code="$2";       shift 2 ;; ## 脚本路径
  --thread )    thread="$2";     shift 2 ;; ## dmu任务并行数
  --alpha )     alpha="$2";      shift 2 ;; ## ssGBLUP relationship matrix alpha
  --vidf )      vidf="$2";       shift 2 ;; ## 验证群个体id文件名
  --dfva )      dfva="$2";       shift 2 ;; ## 加性效应分布方差先验自由度 [3]
  --dfve )      dfve="$2";       shift 2 ;; ## 残差效应分布方差先验自由度 [$dfva]
  --out )       out="$2";        shift 2 ;; ## 输出准确性文件名前缀
  --gebvf )     gebvf="$2";      shift 2 ;; ## Bayes模型中gebvf输出文件夹的名前缀 [EBV]
  --dirPre )    dirPre="$2";     shift 2 ;; ## ebv输出文件夹增加的前缀
  --updatepri ) updatepri="$2";  shift 2 ;; ## 指定轮次更新prior方差尺度矩阵先验 [0]
  --prefix )    prefix="$2";     shift 2 ;; ## 在union/blend/multi等文件夹前添加指定前缀
  --nchr )      nchr="$2" ;      shift 2 ;; ## 染色体数目 [30]
  --dense )     dense=true;      shift   ;; ## 将DMU的ANALYSE中的method设为31，利用多线程计算方差组分
  --suffix )    suffix=true;     shift   ;; ## 在union/blend/multi等文件夹后添加品种名称后缀
  --debug )     debug=true;      shift   ;; ## 不跑DMU、gmatrix、bayes等时间长的步骤
  --rg_local )  rg_local=true;   shift   ;; ## 每个基因组分区的协方差不同，为局部ld、frq相关系数
  --noCov )     noCov=true;      shift   ;; ## 性状间的残差效应约束为0
  --keep_all )  keep_all=true;   shift   ;; ## 保留所有Bayes结果文件
  --traceplot ) traceplot=true;  shift   ;; ## MCMC burn-in期后的trace plot
  --varSamp )   varSamp=true;    shift   ;; ## 输出标记效应方差样本
  --overwri )   overwrite=true;  shift   ;; ## 若存在已生成的准备文件(表型、id文件等)，则进行覆写
  -4 | --dmu4 ) dmu4=true;       shift   ;; ## 使用dmu4进行评估，不估计方差组分
  -h | --help)  grep ";; ##" $0 | grep -v help && exit 1 ;;
  -- ) shift; break ;;
  * ) break ;;
  esac
done

## 调用软件需要的模块
# module load GCC/11.3.0

## 工作目录
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
  code=$(dirname "$script_path")
fi

## 脚本
pheMerge=${code}/R/multibreed_pheno.R
accur_cal=${code}/R/accuracy_bias_calculation.R
fix_frq_ld_bolck=${code}/R/fix_frq_ld_bolck.R
job_pool=${code}/shell/parallel_jobs_func.sh
func=${code}/shell/function.sh
multiG=${code}/R/multibreed_relationship_matrix.R
variance_prior=${code}/R/variance_prior_setting.R
MCMC_polt=${code}/R/MCMC_chain_plot.R
local_rg=${code}/R/local_block_rg.R
# keep_phe_gid=${code}/R/pheno_pre/pheno_all_genoid.R
# phe_order=${code}/R/pheno_pre/pheno_order.R

## 载入自定义函数
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${job_pool}
source ${func}

## 将程序路径加到环境变量中
export PATH=${code}/bin:$PATH

## 检查需要的程序是否在环境变量中能检索到并且可执行
check_command plink gmatrix mbBayesAB LD_mean_r2 run_dmu4 run_dmuai

## 检查需要的脚本文件是否存在且具有执行权限
check_command $pheMerge $accur_cal $fix_frq_ld_bolck $job_pool $func
check_command $multiG $variance_prior $MCMC_polt $local_rg

## 避免执行R脚本时的警告("ignoring environment value of R_HOME")
unset R_HOME

## suffix
if [[ "$suffix" = true ]]; then
  suffix="_$(echo "${pops}" | tr ' ' '_')"
else
  suffix=""
fi

## 存放结果的目录
tpath=${workdir}/${prefix}${type}${suffix}
mkdir -p ${tpath}

## 默认参数
gebvf=${gebvf:=EBV}
nchr=${nchr:=30}
updatepri=${updatepri:=0}
LD_maxD=${LD_maxD:=10000}
r2=${r2:=0}
inter=${inter:=99999}
nsnp_win=${nsnp_win:=100}
win=${win:=50}
bin=${bin:=fix}
ref=${ref:=M}
bincol=${bincol:=1}
fold=${fold:=1}
rep=${rep:=1}
miss=${miss:=-99}
type=${type:=blend}
phereal=${phereal:=1}
h2=${h2:=0.5}
rg=${rg:=0.001}
re=${re:=0.001}
iter=${iter:=30000}
burnin=${burnin:=20000}
thin=${thin:=10}
DIR=${DIR:=${type}}
report_sep=${report_sep:=100}
scale=${scale:="phe"}
GmatM=${GmatM:="single"}
invA=${invA:=1}
alpha=${alpha:=0.05}
method=${method:=GBLUP}
maf=${maf:=-0.01}
vidf=${vidf:=val.id}
r2_merge=${r2_merge:=0.1}
add_rf=${add_rf:=1}
dfva=${dfva:=3}
geno=${geno:=0.2}
mind=${mind:=0.2}
code=${code:=${HOME}/liwn/code}
[[ ${tbvf} && ! ${tbv_col} ]] && tbv_col=2 ## 真实育种值文件存在，默认育种值在文件第二列，id在第一列

## 参数初始化
# [[ ${noCov} ]] && noCov=" --constraint "
[[ ${varSamp} ]] && varSamp=" --varSamp "
[[ ${noCov} ]] && noCov=" --nocov "
[[ ${VaSori} ]] && VaSori=" --VaSori ${VaSori} "
[[ ${R0Sori} ]] && R0Sori=" --R0Sori ${R0Sori} "
[[ ${VaSSori} ]] && VaSSori=" --VaSSori ${VaSSori} "
[[ ${R0SSori} ]] && R0SSori=" --R0SSori ${R0SSori} "
[[ ${VaSd} ]] && VaSd=" --VaSd ${VaSd} "
[[ ${R0Sd} ]] && R0Sd=" --R0Sd ${R0Sd} "
[[ ${priorCor} ]] && priorCor=" --corf ${priorCor} "
[[ ${phyNum} ]] && phyNum=" --phyNum ${phyNum} "
[[ ${dfva} ]] && dfva=" --dfva ${dfve} "
[[ ${dfve} ]] && dfve=" --dfve ${dfve} "
# [[ ${all_samps} ]] && all_samps=" --all_samp_out "
[[ ${keep_all} ]] && keep_all=" --keep_all "
if [[ $(echo ${h2} | awk '{print NF}') -gt 1 ]]; then
  h2B=$(echo ${h2} | awk '{print $2}')
  h2=$(echo ${h2} | awk '{print $1}')
else
  h2B=${h2}
fi

## 随机数种子
if [[ ! ${seed} ]]; then
  seed=$RANDOM
  echo "$seed" >MCMC_seed.txt
fi

## 准确性输出文件名
if [[ ! ${out} ]]; then
  out=accur_GBLUP
  [[ ${type} == "multi" ]] && out=accur_bayes
fi

## 并行作业数设置
if [[ ! ${thread} ]]; then
  if [[ ${SLURM_CPUS_ON_NODE} ]]; then thread=${SLURM_CPUS_ON_NODE}; else thread=$((rep * fold)); fi
fi

## 作业数池
# shellcheck source=/dev/null
source "${job_pool}"

## 作业池初始化
[[ ! ${debug} ]] && job_pool_init ${thread} 0

##################  解析命令行参数  ##############
################################################
np=$(echo ${pops} | awk '{print NF}')
IFS=" " read -r -a popN <<<"$pops"
IFS=" " read -r -a bfileWa <<<"$bfileW"

## DMU多线程
if [[ ${method} == "GBLUP" && ${dense} ]]; then
  MKL_NUM_THREADS=1
  method_dmu="31"
else
  method_dmu="1"
fi

###############  检查品种内评估是否完成  ###########
#################################################
## 交叉验证和重复
rep=$(find ${workdir}/${popN[0]}/val1/rep* -type d | wc -l)
fold=$(find ${workdir}/${popN[0]}/val*/rep1 -type d | wc -l)

## 每个子集下文件是否存在
for r in $(seq 1 ${rep}); do
  for f in $(seq 1 ${fold}); do
    for b in "${popN[@]}"; do
      for file in pheno.txt ${vidf}; do
        [[ ! -s ${workdir}/${b}/val${f}/rep${r}/${file} ]] &&
          echo "${workdir}/${b}/val${f}/rep${r}/${file} not found! " &&
          exit 1

        ## 获取表型文件中的整型和实型变量数
        if [[ ! ${num_int} ]]; then
          phef_within=${workdir}/${b}/val${f}/rep${r}/pheno.txt
          ncol=$(awk 'END{print NF}' ${phef_within})
          for i in $(seq 1 ${ncol}); do
            dot=$(awk -vl=${i} '{print $l}' ${phef_within} | grep -c "\.")
            [[ ${dot} -gt 0 ]] && num_int=$((i - 1)) && break
          done
        fi
        [[ ! ${num_real} ]] && num_real=$(($(awk 'END{print NF}' ${phef_within}) - num_int))
      done
    done
  done
done

##################  获取品种内模型效应设定  #############
########################################################
## 从其中一个品种子集中获取效应设定(不同品种效应相同，以第一个品种为准)
firstB=$(echo ${pops} | cut -d ' ' -f 1)
within_DIR=$(find ${workdir}/${firstB}/val1/rep1 -name "*.DIR" | head -n 1)
[[ ! -s ${within_DIR} ]] && echo "${within_DIR} not found! " && exit 1
model_line=$(grep -n "MODEL" ${within_DIR} | head -n 1 | cut -d ':' -f 1)
all_eff=$(sed -n "$((model_line + 3))p" ${within_DIR} | cut -d ' ' -f '4-')
ran_eff=$(sed -n "$((model_line + 4))p" ${within_DIR} | cut -d ' ' -f '2-')
## 效应个数
nA=$(echo ${all_eff} | awk '{print NF}')
nR=$(echo ${ran_eff} | awk '{print NF}')
## 增加品种固定效应
((nA++))
all_eff="$((num_int + 1)) ${all_eff}"

################  基因型文件  #############
##########################################
## 品种内的基因型文件路径
# IFS=" " read -r -a bfiles <<<"${bfileM} ${bfileW[*]}"
if [[ ! ${debug} ]]; then
  if [[ ! ${bfileW} ]]; then
    IF="\n" mapfile -t bfiles < <(printf "%s\n" "${popN[@]}" | xargs -I {} echo "${workdir}/{}/{}")
  else
    IFS=" " read -r -a bfiles <<<"${bfileW}"
  fi
  ## 检查fam/bim/bed文件是否存在
  check_plink "${bfiles[@]}" ${nchr}
  ## 生成各品种合并后的plink文件
  if [[ ${type} != "multi" && ! -s ${gmat} ]] || [[ ${type} == "multi" ]]; then
    if [[ ${bfileM} ]]; then
      check_plink "${bfileM}" ${nchr}
    else
      ## 合并plink文件
      bfileM=${tpath}/merge
      merge_plink "${bfiles[*]}" ${bfileM}
    fi
  fi
fi

##################  参考群表型文件  ##############
#################################################
$pheMerge \
  --pops "${pops}" \
  --phef "${workdir}/#breed#/val#val#/rep#rep#/pheno.txt" \
  --rep ${rep} \
  --fold ${fold} \
  --nInt ${num_int} \
  --type ${type} \
  --overwri ${overwrite} \
  --pheCol ${phereal} \
  --fixPop \
  --out "${tpath}/val#val#/rep#rep#/pheno.txt"

#####################  系谱文件  #####################
#####################################################
if [[ ${pedf} ]]; then
  for pedi in ${pedf}; do
    if [[ ! -s ${pedi} ]]; then
      echo "${pedi} not found! " && exit 1
    else
      cat ${pedi} >pedi_merge.txt
    fi
  done
  echo "Number of individuals in pedigree: $(wc -l <pedi_merge.txt)"
else
  [[ -s pedi_merge.txt ]] && rm pedi_merge.txt
fi

####################  基因组区域确定  ##################
########################################################
cd ${tpath} || exit
## 基因组分区文件名
if [[ ! -s ${binf} ]]; then
  if [[ ${bin} == "fix" ]]; then
    bin_prefix=${bin}_${nsnp_win}
  elif [[ ${bin} == "frq" || ${bin} == "ld" ]]; then
    bin_prefix=${bin}_${r2_merge}_${nsnp_win}
  else
    echo "${bin} can only be fix, frq, ld! "
    exit 1
  fi
  binf=${tpath}/${bin_prefix}.txt
else
  # binf_base=$(basename "${binf}")
  # bin_prefix="${binf_base%.*}"
  bin_prefix=self_bin
  if [[ ${bincol} != 1 ]]; then
    awk -v col=${bincol} '{print $col}' ${binf} > bin_col1.txt
    binf=$(pwd)/bin_col1.txt
  fi
fi
## 生成文件
if [[ ${type} == 'multi' ]] && [[ ! -s ${binf} || ${overwrite} ]]; then
  echo 'Generating genome regions file...'

  ## 区间划分
  if [[ ${bin} == "fix" ]]; then
    if [[ ${nbin} ]]; then
      nsnp=$(wc -l <${bfileM}.bim)
      ## 约束最终划分的区间数与nbin相近
      nsnp_win=$((nsnp / nbin))
      echo "Number of SNPs per window set to: ${nsnp_win}"
    fi

    [[ ! ${debug} ]] && \
      $fix_frq_ld_bolck \
        --win ${nsnp_win} \
        --map ${bfileM}.bim \
        --bin_merge ${bin} \
        ${phyNum} \
        --out ${binf}
  elif [[ ${bin} == "frq" || ${bin} == "ld" ]]; then
    # ## 检查是否提供了每个品种的基因型文件 
    # [[ ! "${bfileW[*]}" ]] && echo "Required parameter 'bfileW' is missing! " && exit 1

    ## 品种共有标记
    awk '{print $2}' ${bfileM}.bim >SNP_share_id.txt

    ## LD计算(只计算共有标记)
    for i in $(seq 0 $((np - 1))); do
      plink --bfile ${bfiles[i]} \
        --chr-set ${nchr} \
        --extract SNP_share_id.txt \
        --r2 \
        --freq \
        --ld-window-kb ${LD_maxD} \
        --ld-window ${inter} \
        --ld-window-r2 ${r2} \
        --out ${popN[i]} >>${logp}/plink_ld_frq.log
    done

    ## 根据要求定义区间
    $fix_frq_ld_bolck \
      --prefixs "${popN[*]}" \
      --bin_merge ${bin} \
      --win ${nsnp_win} \
      --seg ${r2_merge} \
      --map ${bfileM}.bim \
      --out ${binf}~

    ## 评估软件所需区间文件
    awk '{print $3}' ${binf}~ >${binf}
  fi

  ## 检查是否出错
  if [[ ! -s ${binf} ]]; then
    echo 'error in creat bins file! '
    exit 1
  fi
fi

#####################  方差组分准备  #####################
#########################################################
## 方差组分文件名
if [[ ${priSaDmu} ]]; then
  parfA=$(basename "$(find ${workdir}/${A} -name "*PAROUT" | head -n 1)")
  parfB=$(basename "$(find ${workdir}/${B} -name "*PAROUT" | head -n 1)")

  ## 每个区间协方差先验不同
  [[ ${rg_local} ]] && rg_local=" --rg_local ${binf}~ "

  if [[ ${priSaDmu} == "pheno" ]]; then
    ## 根据表型方差推算遗传和残差方差
    varf_para1="${workdir}/${A}/val#val#/rep#rep#/pheno.txt"
    varf_para2="${workdir}/${B}/val#val#/rep#rep#/pheno.txt"
  elif [[ ${priSaDmu} == "predict" ]]; then
    ## 用品种内评估的方差推算多品种遗传和残差方差
    varf_para1="${workdir}/${A}/val#val#/rep#rep#/${parfA}"
    varf_para2="${workdir}/${B}/val#val#/rep#rep#/${parfB}"
  elif [[ ${priSaDmu} == "A" ]]; then
    ## 用A品种内估计的方差组分
    varf_para1="${workdir}/${A}/val#val#/rep#rep#/${parfA}"
    varf_para2="null"
  elif [[ ${priSaDmu} == "B" ]]; then
    ## 用B品种内估计的方差组分
    varf_para1="${workdir}/${B}/val#val#/rep#rep#/${parfB}"
    varf_para2="null"
  fi

  ## 生成方差组分(初值)文件
  $variance_prior \
    --filea ${varf_para1} \
    --fileb ${varf_para2} \
    --rep ${rep} \
    --fold ${fold} \
    --h2 ${h2} \
    --h2B ${h2B} \
    --rg ${rg} \
    --re ${re} \
    --type ${type} \
    --var ${priSaDmu} \
    --norec \
    --add_rf ${add_rf} \
    --overwri ${overwrite} \
    --pcol $((num_int + phereal)) \
    ${rg_local} \
    --out ${tpath}/val#val#/rep#rep#/${type}_${dirPre}${bin}_prior.txt
fi

#####################  关系矩阵准备  #####################
##########################################################
if [[ ${type} == 'blend' || ${type} == 'union' ]]; then
  ## 生成基因型关系矩阵G阵
  if [[ -s ${gmat} ]]; then
    if [[ ${method} == "ssGBLUP" && ! -s ${gidf} ]]; then
      echo "${gidf} not found! " && exit 1
    else
      echo "Use user provided relationship matrix: ${gmat}"
    fi
  elif [[ ${GmatM} == 'multi' ]]; then
    ## 多品种关系矩阵
    [[ ! "${bfileW[*]}" ]] && echo "Required parameter 'bfileW' is missing! " && exit 1

    :> all_breed.id
    for i in $(seq 0 $((np - 1))); do
      awk '{print $2}' ${bfileWa[i]}.fam >${bfileW[i]}.ids
      cat ${bfileW[i]}.ids >>all_breed.id
      idf="${bfileW[i]}.ids ${idf}"
    done

    ## 格式转换，顺便剔除存在缺失的位点(只存在某个群体中的位点)
    plink \
      --bfile ${bfileM} \
      --chr-set ${nchr} \
      --geno 0 \
      --recode A \
      --out merge >>${logp}/plink_recodeA.log

    ## 多品种关系矩阵（待修改）
    $multiG \
      --rawf merge.raw \
      --idf ${idf} \
      --out merge
    mv merge.grm merge.agrm.id_fmt
  elif [[ ${GmatM} == 'single' ]]; then
    ## 生成基因型关系矩阵G阵
    if [[ ${method} == "GBLUP" ]]; then
      gmat_inv="--inv"
      [[ ! ${gmat} ]] && gmat=${tpath}/merge.agiv.id_fmt
    elif [[ ${method} == "ssGBLUP" ]]; then
      gmat_inv=""
      [[ ! ${gmat} ]] && gmat=${tpath}/merge.agrm.id_fmt
      [[ ! ${gidf} ]] && gidf=${tpath}/merge.id
    fi

    ## 构建基因组关系矩阵
    echo "Read the plink bed file and Calculate the additive G matrix..."
    [[ ! ${debug} ]] && gmatrix --bfile ${bfileM} --grm agrm --out ${tpath}/merge ${gmat_inv}
    echo "G matrix created."

    ## 改名
    if [[ ${method} == "GBLUP" ]]; then
      ## 若指定了gmat文件名，则改名(同时可能移动文件)
      [[ ${gmat} && ! -s ${gmat} ]] && \
        mv merge.agiv.id_fmt ${gmat} && \
        echo "gamt created and has been rename to: ${gmat}"
    elif [[ ${method} == "ssGBLUP" ]]; then
      ## 若指定了gmat文件名，则改名(同时可能移动文件)
      [[ ${gmat} && ! -s ${gmat} ]] && \
        mv merge.agrm.id_fmt ${gmat} && \
        echo "gamt created and has been rename to: ${gmat}"
      [[ ${gidf} && ! -s ${gidf} ]] && mv ${tpath}/merge.id ${gidf}
    fi
  else
    echo "Gmat can only be multi or single! "
    exit 1
  fi
fi

###################  dmu参数卡  ####################
####################################################
cd ${tpath} || exit
if [[ ${type} == 'blend' || ${type} == 'union' ]]; then
  ## ANALYSE
  [[ ${dmu4} ]] && ANALYSE="11 9 0 0" || ANALYSE="1 ${method_dmu} 0 0"

  ## MODEL
  if [[ ${type} == 'blend' ]]; then
    num_real=1
    MODEL="1\n0\n1 0 ${nA} ${all_eff}\n${nR} ${ran_eff}\n0\n0"
  else
    ABSORB="0"
    MODELS="1 0 ${nA} ${all_eff}"
    RANDOMS="${nR} ${ran_eff}"
    REGRES="0"
    NOCOV=""
    nNOCOV=0
    for i in $(seq 2 ${np}); do
      ABSORB="${ABSORB}\n0"
      MODELS="${MODELS}\n${i} 0 ${nA} ${all_eff}"
      RANDOMS="${RANDOMS}\n${nR} ${ran_eff}" # RANDOM为shell关键字，故改名
      REGRES="${REGRES}\n0"
      for j in $(seq 1 $((i -1))); do
        NOCOV="${NOCOV}\n${j} ${i}"
        ((nNOCOV++))
      done
    done

    num_real=${np}
    MODEL="${np}\\n${ABSORB}\\n${MODELS}\\n${RANDOMS}\\n${REGRES}\\n${nNOCOV}${NOCOV}"
  fi

  ## 写出参数卡
  [[ -s ${DIR}.DIR ]] && echo "warn: ${DIR}.DIR has been overwrited! "
  {
    echo "\$COMMENT"
    echo "get EBV of individuals in validation population with reduced phenotypes"
    echo "\$ANALYSE ${ANALYSE}"
    echo "\$DATA  ASCII ($((num_int + 1)), ${num_real}, ${miss}) pheno.txt"
    echo -e "\$MODEL\n${MODEL}"
    echo "\$PRIOR %varf%"
  } >${DIR}.DIR

  ## 方差组分结构
  if [[ ${method} == "GBLUP" ]]; then
    echo "\$VAR_STR ${add_rf} GREL ASCII ${gmat}" >>${DIR}.DIR
    add_sol=3
  elif [[ ${method} == "ssGBLUP" ]]; then
    [[ ! -s pedi_merge.txt ]] && echo "pedi_merge.txt not found! " && exit 1
    echo "\$VAR_STR ${add_rf} PGMIX ${invA} ASCII pedi_merge.txt ${gidf} ${tpath}/merge.agrm.id_fmt ${alpha} G-ADJUST" >>${DIR}.DIR
    add_sol=4
  elif [[ ${method} == "PBLUP" ]]; then
    [[ ! -s pedi_merge.txt ]] && echo "pedi_merge.txt not found! " && exit 1
    echo "\$VAR_STR ${add_rf} PED ${invA} ASCII pedi_merge.txt" >>${DIR}.DIR
    add_sol=4
  else
    echo "method can only be PBLUP/GBLUP/ssGBLUP! "
    exit 1
  fi

  ## 输出效应值
  echo "\$SOLUTION" >>${DIR}.DIR
fi

#####################  子集处理和评估  ###################
##########################################################
for r in $(seq 1 ${rep}); do # r=1;f=1
  for f in $(seq 1 ${fold}); do
    ## 更改工作文件夹
    vali_path=${tpath}/val${f}/rep${r}
    cd ${vali_path} || exit

    ## 如果为调试模式则跳过估计育种值步骤
    [[ ${debug} ]] && continue

    #####  GBLUP模型  #####
    #######################
    if [[ ${type} == 'blend' || ${type} == 'union' ]]; then
      ## 拷贝DMU参数卡
      cp ${tpath}/${DIR}.DIR .

      ## 根据需要提供方差组分(或迭代初值)
      if [[ -s ${type}_${dirPre}${bin}_prior.txt ]]; then
        sed -i "s#%varf%#${type}_${dirPre}${bin}_prior.txt#" ${DIR}.DIR
      else
        sed -i '/$PRIOR/d' ${DIR}.DIR
      fi

      ## 育种值估计
      [[ $(find ${vali_path} -name "*.SOL" 2>/dev/null | wc -l) -ge 1 ]] && continue # temp
      if [[ ${dmu4} ]]; then
        job_pool_run run_dmu4 ${DIR} ${vali_path}
      else
        job_pool_run run_dmuai ${DIR} ${vali_path}
      fi
    fi

    ##### MT-BayesAS模型 #####
    ############################
    if [[ ${type} == 'multi' ]]; then
      ## 固定效应
      fix_eff=${all_eff%" ${ran_eff}"}

      [[ $(find ${vali_path} -name "${gebvf}*.txt" 2>/dev/null | wc -l) -ge 1 ]] && continue # temp
      # --binf ${binf} \ # debug
      job_pool_run mbBayesAB \
        --bfile ${bfileM} \
        --phef ${vali_path}/pheno.txt \
        --fix "${fix_eff}" \
        --binf ${binf} \
        --iter ${iter} \
        --burnin ${burnin} \
        --thin ${thin} \
        --outdir ${vali_path} \
        --report ${report_sep} \
        --varOut var_${bin}_${r2_merge}_${nsnp_win}.txt \
        --effOut effect_${bin}_${r2_merge}_${nsnp_win}.txt \
        --gebvOut ${gebvf} \
        --mcmcOut MCMC_process_${bin}.txt \
        --seed ${seed} \
        --logf ${bin}_${r2_merge}_${nsnp_win}_gibs_${SLURM_JOB_ID}.log \
        ${noCov} \
        ${dfva} \
        ${dfve} \
        ${VaSori} \
        ${R0Sori} \
        ${VaSSori} \
        ${R0SSori} \
        ${VaSd} \
        ${R0Sd} \
        ${varSamp} \
        ${priorCor}
    fi
  done
done

## 等待后台程序运行完毕
[[ ! ${debug} ]] && job_pool_wait
## 释放线程
[[ ! ${debug} ]] && job_pool_shutdown

####################  准确性计算  ##################
###################################################
## 更换工作路径
cd ${tpath} || exit
## 准确性计算参数
if [[ ${tbv_col} == "same" ]]; then
  option=" --tbv_col $((num_int + phereal))"
elif [[ ${tbv_col} ]]; then
  option=" --tbv_col ${tbv_col}"
fi
[[ ${tbvf} ]] && option="${option} --tbvf ${tbvf} "
## 其他参数
if [[ ${type} == 'multi' ]]; then
  ## MT-Bayes模型
  option="${option} --ebvf ${tpath}/val#val#/rep#rep#/${gebvf}_y%i%.txt"
  ebv_col=2
else
  ## MT-GBLUP模型
  option="${option} --add_sol ${add_sol} --dir_val ${tpath}/val#val#/rep#rep#/${DIR}"
  if [[ ${type} == 'blend' ]]; then
    ebv_col=1
  else
    ebv_col=%i%
  fi
fi
## 计算准确性
for i in $(seq 0 $((np - 1))); do
  ip1=$((i+1))
  ## 输出文件名
  if [[ ${type} == 'multi' ]]; then
    outf=${out}_${dirPre}_${bin}_${popN[${i}]}.txt
  else
    outf=${out}_${popN[${i}]}.txt
  fi

  ## 计算准确性
  $accur_cal \
    --ebv_col ${ebv_col/\%i\%/${ip1}} \
    --val_idf ${workdir}/${popN[${i}]}/val#val#/rep#rep#/${vidf} \
    --fid ${popN[${i}]} \
    --rep ${rep} \
    --fold ${fold} \
    ${option/\%i\%/${ip1}} \
    --out ${tpath}/${outf}
done

#################  BayesAS 模型收敛作图  #################
#########################################################
## BayesAS 模型 收敛作图
if [[ ${type} == 'multi' ]]; then
  if [[ ${traceplot} ]]; then
    for r in $(seq 1 ${rep}); do
      ## Bayes MCMC链作图
      for f in $(seq 1 ${fold}); do
        $MCMC_polt \
          --files ${vali_path}/${bin}_MCMC_process.txt \
          --start ${burnin} \
          --end ${iter} \
          --thin ${report_sep} \
          --names "μ1 μ2 alpha11 alpha12" \
          --out "${vali_path}/${bin}_MCMC_process"
      done
    done
  fi
fi
