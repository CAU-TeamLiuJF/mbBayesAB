#!/usr/bin/bash
#SBATCH --job-name=GP_CV

########################################################################################################################
## 版本: 1.1.1
## 作者: 李伟宁 liwn@cau.edu.cn
## 日期: 2023-09-09
## 
## 功能：
## 1.根据提供的plink格式的基因型文件进行表型模拟
## 2.根据提供的表型和基因型(系谱)信息进行"品种内"遗传评估并计算评估准确性
## 3.                              "多品种"
## 
## 使用: GP_cross_validation.sh --help
## 
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################



###################  参数处理  #####################
###################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## 参数名
TEMP=$(getopt -o h --long code:,proj:,type:,breeds:,thread:,varf:,traits:,trait:,h2:,vg:,rg_sim:,rg_pri:,rg_dist:,means:,method:,phef:,pedf:,binf:,nqtl:,nbin_cor:,nsnp_cor:,nsnp_win:,phyNum:,VaSori:,R0Sori:,VaSSori:,R0SSori:,priorCor:,VaSd:,R0Sd:,bin_sim:,tbv_col:,all_eff:,ran_eff:,iter:,burnin:,ref:,dirPre:,nbin:,min:,bfile:,seed:,fold:,rep:,gen:,nsnp:,nsnp_sim:,out:,sim_dir:,nginds:,seg_gens:,extentLDs:,last_males:,last_females:,founder_sel:,seg_sel:,last_sel:,last_litters:,geno_gen:,maf:,binDiv:,binThr:,nchr:,nmloc:,nqloci:,QMSim_h2:,gebvf:,dfva:,dfve:,varSamp,debug,suffix,dense,noCov,append,overlap,evenly,ped_var,dmu4,help \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
## 解析参数
while true; do
  case "$1" in
    --proj )     proj="$2";     shift 2 ;; ## 项目目录 [必需]
    --breeds )   breeds="$2";   shift 2 ;; ## 群体/品种标识符，如'YY DD' [必需]
    --traits )   traits="$2";   shift 2 ;; ## 表型文件中所有性状名，如'AGE BF' [必需]
    --type )     type="$2";     shift 2 ;; ## 分析类型，bin/adj/gsim/geno/psim/within/blend/union/multi/accur/var [必需]
    --trait )    trait="$2";    shift 2 ;; ## 需要进行分析的性状名，如'AGE' [traits]
    --bfile )    bfile="$2";    shift 2 ;; ## plink文件前缀，如"/public/home/merge" [NULL]
    --phef )     phef="$2";     shift 2 ;; ## 表型文件 [NULL]
    --pedf )     pedf="$2";     shift 2 ;; ## 系谱文件 [NULL]
    --varf )     varf="$2";     shift 2 ;; ## 方差组分文件(如没提供则用系谱/SNP信息估计)
    --dfva )     dfva="$2";     shift 2 ;; ## 加性效应分布方差先验自由度 [3]
    --dfve )     dfve="$2";     shift 2 ;; ## 残差效应分布方差先验自由度 [$dfva]
    --rg_dist )  rg_dist="$2";  shift 2 ;; ## 加性遗传相关大小服从的分布，uniform/normal/identical [identical]
    --rg_sim )   rg_sim="$2";   shift 2 ;; ## 性状模拟时群体间的加性遗传相关大小 [0.5]
    --rg_pri )   rg_pri="$2";   shift 2 ;; ## 多性状模型中群体间的加性遗传相关大小先验 [NULL]
    --h2 )       h2="$2";       shift 2 ;; ## 性状模拟时各品种性状的遗传力大小 [0.2]
    --vg )       vg="$2";       shift 2 ;; ## QTL效应方差 [1.0]
    --means )    means="$2";    shift 2 ;; ## 性状模拟时各品种群体均值 [NULL]
    --iter )     iter="$2";     shift 2 ;; ## MCMC 总循环次数 [30000]
    --burnin )   burnin="$2";   shift 2 ;; ## MCMC burn-in循环数 [20000]
    --ref )      ref="$2";      shift 2 ;; ## 用于区间划分的SNP面板，M/1/2/... [M]
    --dirPre )   dirPre="$2";   shift 2 ;; ## multi情形中在文件夹增加的前缀，如pre_multi_A_B [NULL]
    --method)    method="$2";   shift 2 ;; ## 育种值估计方法，可为PBLUP/GBLUP/ssGBLUP/BayesAS [GBLUP]
    --nsnp_cor ) nsnp_cor="$2"; shift 2 ;; ## 当设定指定数目区间内品种间存在遗传相关时，区间内的QTL数目 [10]
    --nbin_cor ) nbin_cor="$2"; shift 2 ;; ## 设定指定数目区间内品种间存在遗传相关 [10]
    --nsnp_sim ) nsnp_sim="$2"; shift 2 ;; ## 性状模拟时，以固定数目进行区间划分时，每个区间内所含的标记数 [60]
    --bin_sim )  bin_sim="$2";  shift 2 ;; ## 性状模拟时区间划分方法，binf路径/win/chr [win]
    --nsnp_win ) nsnp_win="$2"; shift 2 ;; ## 多品种评估以固定数目进行区间划分时，每个区间内所含的标记数 [100]
    --phyNum )   phyNum="$2";   shift 2 ;; ## 以固定SNP数还是物理位置position划分区间，num/pos [num]
    --nqtl )     nqtl="$2";     shift 2 ;; ## 性状模拟时的QTL数目 [300]
    --nbin )     nbin="$2";     shift 2 ;; ## 多品种评估以固定数目SNP为区间划分依据时，尽可能接近这个设置的区间数 [NULL]
    --min )      min="$2";      shift 2 ;; ## 保证挑选QTL的区间中SNP标记数的最小值 [NULL]
    --tbv_col )  tbv_col="$2";  shift 2 ;; ## 不进行校正表型计算，真实育种值在表型文件中tbv_col列 [NULL]
    --VaSori )   VaSori="$2";   shift 2 ;; ## 方差组分先验信息来源，phe/iden/FILE [phe]
    --R0Sori )   R0Sori="$2";   shift 2 ;; ## 方差组分先验信息来源，phe/iden/FILE [$VaSori]
    --VaSSori )  VaSSori="$2";  shift 2 ;; ## 方差组分先验分布的尺度矩阵信息来源，phe/iden [phe]
    --R0SSori )  R0SSori="$2";  shift 2 ;; ## 方差组分先验分布的尺度矩阵信息来源，phe/iden [$VaSSori]
    --priorCor ) priorCor="$2"; shift 2 ;; ## 相关系数先验文件 [NULL]
    --VaSd )     VaSd="$2";     shift 2 ;; ## 加性效应先验分布尺度参数服从的分布 NULL/wishart/invgamma [NULL]
    --R0Sd )     R0Sd="$2";     shift 2 ;; ## 残差效应先验分布尺度参数服从的分布 NULL/wishart [NULL]  
    --all_eff )  all_eff="$2";  shift 2 ;; ## 固定效应和随机效应所在列，如"2 1" ["2 1"]
    --ran_eff )  ran_eff="$2";  shift 2 ;; ## 随机效应所在列，如"1" [1]
    --binf )     binf="$2";     shift 2 ;; ## 区间划分文件 [NULL]
    --rep)       rep="$2";      shift 2 ;; ## 交叉验证重复数 [1]
    --fold)      fold="$2";     shift 2 ;; ## 交叉验证折数 [NULL]
    --gen)       gen="$2";      shift 2 ;; ## 有表型个体中倒数gen个世代个体为验证群 [ALL]
    --seed )     seed="$2";     shift 2 ;; ## MCMC抽样及验证群划分时的随机种子 [40296]
    --code )     code="$2";     shift 2 ;; ## 脚本文件所在目录，如/BIGDATA2/cau_jfliu_2/liwn/code [NULL]
    --thread )   thread="$2";   shift 2 ;; ## 线程数 [NULL]
    --gebvf )    gebvf="$2";    shift 2 ;; ## Bayes模型中gebvf输出文件夹的名前缀 [EBV]
    --sim_dir )      sim_dir="$2";      shift 2 ;; ## 模拟结果文件输出文件夹名，注意该文件夹不能已存在 [rep1]
    --nginds )       nginds="$2";       shift 2 ;; ## 各品种在输出基因型的群体中选择的基因型个体数 ["600 ..."]
    --seg_gens )     seg_gens="$2";     shift 2 ;; ## 各品种从历史群体中分离后，分离的世代数 ["40 10"]
    --extentLDs )    extentLDs="$2";    shift 2 ;; ## 各品种最后稳定LD经历的世代数 ["10 10"]
    --last_males )   last_males="$2";   shift 2 ;; ## 各品种最后一个阶段的群体中的雄性个体数 ["100 10"]
    --last_females ) last_females="$2"; shift 2 ;; ## 各品种最后一个阶段的群体中的雌性个体数 ["500 50"]
    --founder_sel )  founder_sel="$2";  shift 2 ;; ## 各品种从历史群体中选择个体的依据 ["tbv /h,tbv /l"]
    --seg_sel )      seg_sel="$2";      shift 2 ;; ## 各品种在世代选择阶段的个体选留依据 ["phen /h,phen /l"]
    --last_sel )     last_sel="$2";     shift 2 ;; ## 各品种在LD稳定阶段的个体选留依据 ["rnd,rnd"]
    --last_litters ) last_litters="$2"; shift 2 ;; ## 各品种在LD稳定阶段的每窝个体数 ["10 10"]
    --geno_gen )     geno_gen="$2";     shift 2 ;; ## 输出基因型个体的世代 [8-10]
    --nchr )         nchr="$2";         shift 2 ;; ## 染色体数目 [30]
    --nmloc )        nmloc="$2";        shift 2 ;; ## 每条染色体上的标记数 [300000]
    --nqloci )       nqloci="$2";       shift 2 ;; ## 每条染色体上的QTL数 [100]
    --QMSim_h2 )     QMSim_h2="$2";     shift 2 ;; ## 性状的广义遗传力大小 [0.3]
    --maf )          maf="$2";          shift 2 ;; ## 品种内模拟基因型标记选择标准 [0.01]
    --binDiv )       binDiv="$2";       shift 2 ;; ## 抽样SNP时区间划分的依据，pos/frq [pos]
    --binThr )       binThr="$2";       shift 2 ;; ## 抽样SNP时区间划分长度，物理位置cM或基因频率步长 [10]
    --out )          out="$2";          shift 2 ;; ## 输出文件名 [依分析类型而定]
    --nsnp )         nsnp="$2";         shift 2 ;; ## 模拟基因型中需要选择的标记个数 [50000]
    --dense )        dense=true;        shift   ;; ## 用DMUAI模块时，应用稠密矩阵算法（多线程）即method为31
    --rg_local )     rg_local=true;     shift   ;; ## 每个基因组分区的协方差不同，为局部ld、frq相关系数
    --all_comb )     all_comb=true;     shift   ;; ## breeds中所有可能的品种组合都进行评估
    --noCov )        noCov=true;        shift   ;; ## 性状间的残差效应约束为0
    --overlap )      overlap=true;      shift   ;; ## SNP标记中包含QTL
    --suffix )       suffix=true;       shift   ;; ## 在union/blend/multi等文件夹后添加品种名称后缀，如blend_YY_LL
    --evenly )       evenly=true;       shift   ;; ## 在一个区间内挑选SNP作为存在遗传相关的QTL时，使SNP在区间中均匀分布
    --append )       append=true;       shift   ;; ## 计算校正表型时，把各个群体的校正表型进行合并
    --varSamp )      varSamp=true;      shift   ;; ## 输出标记效应方差样本
    --debug )        debug=true;        shift   ;; ## 不跑DMU、gmatrix、bayes等时间长的步骤
  -4 | --dmu4 )     dmu4=true;    shift   ;;  # 用DMU4模型估计育种值
  -v | --ped_var )   ped_var=true;      shift   ;; ## 使用ssGBLUP时用系谱信息估计方差组分
  -h | --help)    grep ";; ##" $0 | grep -v help && exit 1 ;;
  -- ) shift; break ;;
  * ) shift; break ;;
  esac
done

## 保存当前的工作路径
workdir=$(pwd)

## 检查必要参数是否提供
if [[ ! ${proj} ]]; then
  echo "parameter proj must be provided! "
  exit 1
# elif [[ ! ${breeds} ]]; then
#   echo "parameter --breeds is reduired! "
#   exit 1
fi

[[ ! -d ${proj} ]] && mkdir -p ${proj}

## 日志文件夹
logp=${proj}/log
mkdir -p ${logp}

## 避免执行R脚本时的警告("ignoring environment value of R_HOME")
unset R_HOME

## 脚本所在文件夹
if [[ ${code} ]]; then
  [[ ! -d ${code} ]] && echo "${code} not exists! " && exit 5
else
  script_path=$(dirname "$(readlink -f "$0")")
  code=$(dirname "$script_path")
fi

## 脚本
phe_adj=${code}/shell/dmu_get_pheno_adj.sh
GP_single=${code}/shell/GP_single_breed.sh
GP_multi=${code}/shell/GP_multi_breed.sh
genomic_structure=${code}/shell/genomic_structure.sh
accuracy_summ=${code}/shell/accuracy_summary.sh
varComp_summ=${code}/shell/varcomp_summary.sh
func=${code}/shell/function.sh

## 将程序路径加到环境变量中
export PATH=${code}/bin:$PATH

## 加载自定义函数
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${func}

## 检查需要的程序是否在环境变量中能检索到并且可执行
check_command plink gmatrix mbBayesAB LD_mean_r2 run_dmu4 run_dmuai

## 检查需要的脚本文件是否存在且具有执行权限
check_command $GP_single $GP_multi $phe_adj

## 默认参数
gebvf=${gebvf:=EBV}
method=${method:="GBLUP"}
all_eff=${all_eff:="2 1"}
ran_eff=${ran_eff:="1"}
seed=${seed:="8123"}
h2=${h2:="0.2"}
iter=${iter:="30000"}
burnin=${burnin:="20000"}
ref=${ref:="M"}
nqtl=${nqtl:="300"}
nsnp_cor=${nsnp_cor:="10"}
nbin_cor=${nbin_cor:="10"}
nsnp_sim=${nsnp_sim:="60"}
nsnp=${nsnp:="50000"}
nsnp_win=${nsnp_win:="100"}
bin_sim=${bin_sim:="win"}
type=${type:="sim"}
trait=${trait:="${traits}"}
QMSim_rep=${QMSim_rep:="1"}
nchr=${nchr:="30"}
nmloc=${nmloc:="300000"}
nqloci=${nqloci:="100"}
QMSim_h2=${QMSim_h2:="0.3"}
QMSim_qtlh2=${QMSim_qtlh2:="${QMSim_h2}"}
bottleneck=${bottleneck:="250"}
prmpath=${prmpath:="${code}/prm"}
sim_dir=${sim_dir:="rep1"}
geno_gen=${geno_gen:="8-10"}
SLURM_JOB_ID=${SLURM_JOB_ID:="${RANDOM}"}
binDiv=${binDiv:="pos"}
binThr=${binThr:="10"}
min=${min:="${nsnp_cor}"}
maf=${maf:=-0.01}
dfva=${dfva:=3}
dfve=${dfve:="${dfva}"}
VaSori=${VaSori:=phe}
VaSSori=${VaSSori:=phe}
R0Sori=${R0Sori:="${VaSori}"}
R0SSori=${R0SSori:="${VaSSori}"}
R0Sd=${R0Sd:="${VaSd}"}

## 不同类型分析下的默认参数
if [[ ${type} == "var" || ${type} == "accur" ]]; then
  ## 默认参数
  rg_sim=${rg_sim:=/}
  vg=${vg:=/}
  rgv=${rgv:=/}
  ref=${ref:=/}
  rep=${rep:=/}
  rg_dist=${rg_dist:=/}
  [[ ${out} ]] && out=" --out ${out} "
elif [[ ${type} == "within" ]]; then
  if [[ ${method} == "GBLUP" ]]; then
    out=${out:="accur_GBLUP.txt"}
  # else
    # out=${out:="accur_BayesAS.txt"}
  fi
fi

## 其他默认参数
vg=${vg:="1.0"}
rgv=${rgv:="1"}
rg_sim=${rg_sim:="0.5"}
rg_dist=${rg_dist:=identical}

## 从参数获取信息
read -ra breeds_array <<<"$breeds"
read -ra traits_array <<<"$traits"
read -ra trait_array <<<"$trait"
read -ra sim_dirs <<<"$sim_dir"
nbreed=${#breeds_array[@]}
ntrait=${#trait_array[@]}
[[ ${varSamp} ]] && varSamp=" --varSamp "
[[ ${overlap} ]] && overlap=" --overlap "
[[ ${evenly} ]] && evenly=" --evenly "
[[ ${dense} ]] && dense=" --dense "
[[ ${dmu4} ]] && dmu4=" --dmu4 "
[[ ${debug} ]] && debug=" --debug "
[[ ${suffix} ]] && suffix=" --suffix "
[[ ${append} ]] && append=" --append "
[[ ${noCov} ]] && noCov=" --res_const "
[[ ${rg_local} ]] && rg_local=" --rg_local "
[[ ${varf} ]] && varf=" --varf ${varf} "
[[ ${dfva} ]] && dfva=" --dfva ${dfva} "
[[ ${dfve} ]] && dfve=" --dfve ${dfve} "
[[ ${phyNum} ]] && phyNum=" --phyNum ${phyNum} "
[[ ${gen} ]] && gen=" --gen ${gen} "
[[ ${binf} ]] && binf=" --binf ${binf} "
[[ ${ped_var} ]] && ped_var=" --ped_var "
[[ ${pedf} ]] && pedf=" --pedf ${pedf} "
[[ ${tbv_col} ]] && tbv_col=" --tbv_col ${tbv_col} "
[[ ${nbin} ]] && nbin=" --nbin ${nbin} "
[[ ${rg} ]] && rg=" --rg ${rg} "
[[ ${dirPre} ]] && dirPre=" --prefix ${dirPre} "
[[ ${bfile} ]] && bfile_para=" --bfile ${bfile} "
[[ ${priorCor} ]] && priorCor=" --priorCor ${priorCor} "

if [[ ${VaSd} ]]; then
  if [[ ${VaSd} != "none" ]]; then
    VaSd=" --VaSd ${VaSd} "
  else
    unset VaSd
  fi
fi


if [[ ${R0Sd} ]]; then
  if [[ ${R0Sd} != "none" ]]; then
    R0Sd=" --R0Sd ${R0Sd} "
  else
    unset R0Sd
  fi
fi

## 由品种数确定的默认参数
nginds=${nginds:=$(printf "%${nbreed}s" | sed "s/ /600 /g" | sed 's/ *$//')}
last_litters=${last_litters:=$(printf "%${nbreed}s" | sed "s/ /10 /g" | sed 's/ *$//')}
last_females=${last_females:=$(printf "%${nbreed}s" | sed "s/ /200 /g" | sed 's/ *$//')}
founder_sel=${founder_sel:=$(printf "%${nbreed}s" | sed "s/ /rnd,/g" | sed 's/,$//')}
seg_sel=${seg_sel:=${founder_sel}}
last_sel=${last_sel:=${founder_sel}}
extentLDs=${extentLDs:=$(printf "%${nbreed}s" | sed "s/ /10 /g" | sed 's/ *$//')}
last_males=${last_males:=$(printf "%${nbreed}s" | sed "s/ /40 /g" | sed 's/ *$//')}
seg_gens=${seg_gens:=$(printf "%${nbreed}s" | sed "s/ /40 /g" | sed 's/ *$//')}

## 模拟表型时，性状数为1，无性状名
if [[ ${ntrait} == "0" ]]; then
  ntrait=1
  trait_array=('/')
  traits=/
fi

## 修改作业名
if [[ ${SLURM_CPUS_ON_NODE} ]]; then
  thread=${thread:=${SLURM_CPUS_ON_NODE}}

  ## 修改作业名
  scontrol update jobid=${SLURM_JOB_ID} name="${type}"
fi

## 并行作业数
[[ ! ${thread} && ${rep} && ${fold} ]] && thread=$((rep * fold))

## 修改工作文件夹
cd ${proj} || exit 5

## 从合并群中提取各个品种基因型信息(需要map/ped文件用于生成模拟表型)
if [[ ${bfile} ]]; then
  check_plink "${bfile}" ${nchr}
  unset bfiles
  for b in "${breeds_array[@]}"; do
    echo ${b} >fid.txt
    [[ ! -s ${b}m.map ]] &&
      plink --bfile ${bfile} --keep-fam fid.txt --chr-set ${nchr} --freq --recode --out ${b}m &>/dev/null
    bfiles="${bfiles} ${proj}/${b}m"
  done
  [[ -s fid.txt ]] && rm fid.txt
fi

## 文件夹
# mkdir -p ${proj}
# cd ${proj} || exit

## 日志文件
logf=${logp}/${type}_${SLURM_JOB_ID}.log

if [[ ${type} == "adj" ]]; then
  ## 计算校正表型（ebv+re）
  for ti in "${trait_array[@]}"; do # ti=MS;pi=1
    ## 获取要评估性状在所有性状中的位置索引
    pi=$(echo "${traits}" | tr ' ' '\n' | grep -n -w -m1 "${ti}" | cut -d':' -f1)
    phedir=${proj}/${ti}

    ## 修改作业名
    [[ ${SLURM_CPUS_ON_NODE} ]] && scontrol update jobid=${SLURM_JOB_ID} name="${type}_${ti}"

    ## 分析指定表型的文件夹
    mkdir -p ${phedir}
    cd "${phedir}" || exit

    [[ -s ${phedir}/phe_adj_BLUP.txt ]] && \
      echo "Warning: phe_adj_BLUP.txt existe! "

    for b in "${breeds_array[@]}"; do
      ## 分析指定表型中指定品种的文件夹
      mkdir -p ${phedir}/${b}
      cd ${phedir}/${b} || exit

      if [[ ${bfile} ]]; then
        if [[ ! -s ${b}_dmu.txt ]]; then
          ## 检查基因型文件是否存在及是否为二进制文件
          check_plink "${bfile}" ${nchr}

          ## 检查plink文件中有无该品种标识的家系id
          if [[ $(grep -c "${b}" ${bfile}.fam) -eq 0 ]]; then
            echo "no family id ${b} in ${bfile}.fam file! "
            exit 1
          fi

          ## 提取指定品种基因型(可能需要对基因型文件进行修改，所以也是复制基因型)
          echo ${b} >tmp_fid.txt
          plink --bfile ${bfile} --keep-fam tmp_fid.txt --chr-set ${nchr} --make-bed --out ${b} &>>${logf}
          rm tmp_fid.txt

          ## 品种相应表型文件
          awk 'FNR==NR{a[$2];next} $1 in a' ${b}.fam ${phef} >${b}_dmu.txt
          # grep ${b} ${phef} | awk '{$NF="";print}' > ${b}_dmu.txt # 表型文件中最后一列为品种
        fi
      else
        cp ${phef} ${b}_dmu.txt
      fi

      ## 计算校正表型
      $phe_adj \
        --phereal ${pi} \
        --DIR phe_adj_BLUP \
        --phef ${b}_dmu.txt \
        --all_eff "${all_eff}" \
        --ran_eff "${ran_eff}" \
        --out ${phedir}/phe_adj_BLUP.txt \
        ${bfile_para} \
        ${ped_var} \
        ${pedf} \
        ${dense} \
        ${varf} \
        ${dmu4} \
        ${debug} \
        ${append} &>>${logf}
    done
  done
elif [[ ${type} == "struct" ]]; then
  $genomic_structure \
    --proj "${proj}" \
    --bfile "${bfile}" \
    --breeds "${breeds}" \
    --maf ${maf}
elif [[ ${type} == "within" ]]; then
  ## 群体内评估
  for ti in "${trait_array[@]}"; do # ti=/;pi=1
    ## 获取要评估性状在所有性状中的位置索引
    pi=$(echo "${traits}" | tr ' ' '\n' | grep -n -w -m1 "${ti}" | cut -d':' -f1)
    phedir=${proj}/${ti}

    ## 真实育种值所在文件
    if [[ -s ${phedir}/phe_adj_BLUP.txt ]]; then
      tbvf=${phedir}/phe_adj_BLUP.txt
    else
      tbvf=${phef}
    fi

    for b in "${breeds_array[@]}"; do
      ## 修改作业名
      [[ ${SLURM_CPUS_ON_NODE} ]] && scontrol update jobid=${SLURM_JOB_ID} name="${type}_${ti}_${b}"

      ## 工作文件夹
      mkdir -p ${phedir}/${b}
      cd ${phedir}/${b} || exit 5

      ## 基因型文件
      if [[ -s ${phedir}/qtl_snpid.txt ]]; then
        ## 在基因型文件中去除qtl
        plink \
          --file ${phedir}/${b}m \
          --exclude ${phedir}/qtl_snpid.txt \
          --make-bed \
          --out ${phedir}/${b}mq &>>${logf}
        bfile=${phedir}/${b}mq
      else
        bfile=${proj}/${b}m
      fi

      ## 表型文件
      if [[ -s ${phedir}/pheno_sim.txt ]]; then
        grep "${b}" ../pheno_sim.txt | awk '{print $2="", $0}' | awk '$2="1"' >${b}_dmu_pheno.txt
        phef=${b}_dmu_pheno.txt
      fi

      ## 计算准确性
      $GP_single \
        --label ${b} \
        --phef ${phef} \
        --bfile ${bfile} \
        --method ${method} \
        --seed ${seed} \
        --all_eff "${all_eff}" \
        --ran_eff "${ran_eff}" \
        --rep ${rep} \
        --fold ${fold} \
        --phereal ${pi} \
        --thread ${thread} \
        --tbvf ${tbvf} \
        --gebvf ${gebvf} \
        ${pedf} \
        ${dense} \
        ${debug} \
        ${gen} \
        ${binf} \
        ${tbv_col} \
        --out ${out} &>>${logf}
    done
  done
elif [[ ${type} == "blend" || ${type} == "union" || ${type} == "multi" ]]; then
  for ti in "${trait_array[@]}"; do # ti=${trait_array[0]}
    ## 获取要评估性状在所有性状中的位置索引
    pi=$(echo "${traits}" | tr ' ' '\n' | grep -n -w -m1 "${ti}" | cut -d':' -f1)
    phedir=${proj}/${ti}

    ## 切换到工作文件夹
    cd ${phedir} || exit 5

    ## 修改作业名
    # mtype=$(basename "$(dirname "$phedir")")_$(basename "$phedir" | tr '/' '_')
    if [[ ${SLURM_CPUS_ON_NODE} ]]; then
      if [[ ${type} == "multi" ]]; then
        name="${type}_${ti}"
      else
        name="${type}_${ti}"
      fi
      scontrol update jobid=${SLURM_JOB_ID} name="${name}"
    fi

    ## 真实育种值所在文件
    if [[ -s ${phedir}/phe_adj_BLUP.txt ]]; then
      ## 校正表型
      tbvf=${phedir}/phe_adj_BLUP.txt
    else
      ## 表型文件中含"真实"育种值
      tbvf=${phef}
    fi

    ## 检查该品种组合中是否所有品种都完成品种内评估，否则跳过
    for breed in $breeds; do
      [[ ! -d ${phedir}/${breed} ]] && echo "${phedir}/${breed} not found." && not_run=true && break
    done
    [[ ${not_run} ]] && unset not_run && continue

    $GP_multi \
      --pops "${breeds}" \
      --type ${type} \
      --tbvf ${tbvf} \
      --phereal ${pi} \
      --nsnp_win ${nsnp_win} \
      --thread ${thread} \
      --iter ${iter} \
      --burnin ${burnin} \
      --ref ${ref} \
      --seed ${seed} \
      --gebvf ${gebvf} \
      --VaSori ${VaSori} \
      --R0Sori ${R0Sori} \
      --VaSSori ${VaSSori} \
      --R0SSori ${R0SSori} \
      ${phyNum} \
      ${priorCor} \
      ${dirPre} \
      ${rg_local} \
      ${binf} \
      ${rg} \
      ${suffix} \
      ${nbin} \
      ${tbv_col} \
      ${varSamp} \
      ${noCov} \
      ${dense} \
      ${dfva} \
      ${dfve} \
      ${VaSd} \
      ${R0Sd} \
      ${debug} &>>${logf}
  done
elif [[ ${type} == "accur" ]]; then
  ## 准确性统计
  $accuracy_summ \
    --proj "${proj}" \
    --breeds "${breeds}" \
    --traits "${traits}" \
    --code "${code}" \
    --cor "${rg_sim}" \
    --vg "${vg}" \
    --rgv "${rg_sim}" \
    --ref "${rg_sim}" \
    --rep "${rep}" \
    --dist "${rg_dist}" \
    ${out}
elif [[ ${type} == "var" ]]; then
  ## 遗传参数统计
  $varComp_summ \
    --proj "${proj}" \
    --breeds "${breeds}" \
    --traits "${traits}" \
    --code "${code}" \
    --vg "${vg}" \
    --rgv "${rg_sim}" \
    --ref "${rg_sim}" \
    --cor "${rg_sim}" \
    --rep "${rep}" \
    --dist "${rg_dist}" \
    ${out}
fi

# ## 删除没有信息的日志文件
# [[ $? -eq 0 ]] && [[ -s ${logf} ]] && rm ${logf}
