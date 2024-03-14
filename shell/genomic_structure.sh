#!/usr/bin/bash


########################################################################################################################
## 版本: 1.1.0
## 作者: 李伟宁 liwn@cau.edu.cn
## 日期: 2024-02-17
## 
## 根据基因型文件计算PCA、LD、基因频率和LD相关性
## 
## 使用: ./genomic_structure.sh.sh --help
## 
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


###################  参数处理  #####################
####################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## 参数名
TEMP=$(getopt -o h --long bfile:,proj:,breeds:,maf:,code:,win:,inter:,out:,wqc,help \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
## 解析参数
while true; do
  case "$1" in
    --bfile )  bfile="$2";  shift 2 ;; ## plink二进制基因型文件路径 [必要参数]
    --proj )   proj="$2";   shift 2 ;; ## 项目目录 [NULL]
    --breeds ) breeds="$2"; shift 2 ;; ## 品种标识符，如不提供则从基因型文件中获取 [NULL]
    --maf )    maf="$2";    shift 2 ;; ## 抽样SNP时允许的最小等位基因频率 [0.05]
    --win )    win="$2" ;   shift 2 ;; ## 计算LD时的SNP之间的最大距离/kb [10000]
    --inter )  inter="$2" ; shift 2 ;; ## 计算LD时的SNP之间间隔的最大SNP数 [99999]
    --r2 )     r2="$2" ;    shift 2 ;; ## 计算LD时需要过滤的最小r2值 [0]
    --nchr )   nchr="$2" ;  shift 2 ;; ## 染色体数目 [30]
    --code )   code="$2";   shift 2 ;; ## 脚本文件所在目录，如/BIGDATA2/cau_jfliu_2/liwn/code [NULL]
    --out )    out="$2";    shift 2 ;; ## 最终输出基因型文件前缀 [merge]
    --wqc )    wqc=true;    shift   ;; ## 品种内质控
  -h | --help)    grep ";; ##" $0 | grep -v help && exit 1 ;;
  -- ) shift; break ;;
  * ) shift; break ;;
  esac
done

## 检查必要参数是否提供
if [[ ! ${bfile} ]]; then
  echo "genotype file must be provied by --bfile! "
  exit 1
fi

## 避免执行R脚本时的警告("ignoring environment value of R_HOME")
unset R_HOME

## 工作路径
if [[ ${proj} ]]; then
  cd ${proj} || exit
fi

## 脚本所在文件夹
if [[ ${code} ]]; then
  [[ ! -d ${code} ]] && echo "${code} not exists! " && exit 5
else
  script_path=$(dirname "$(readlink -f "$0")")
  code=$(dirname "$script_path")
fi

## 脚本文件
func=${code}/shell/function.sh
PCA_cal=${code}/shell/pca_multi_pop.sh
PCA_plot=${code}/R/PCA_plot.R
geno_dist=${code}/shell/distance_multi_pop.sh
LD_cor=${code}/R/LD_decay_plot.R
corr_cal=${code}/R/columns_correlation.R

## 将程序路径加到环境变量中
export PATH=${code}/bin:$PATH

## 加载自定义函数
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${func}

## 检查需要的程序是否在环境变量中能检索到并且可执行
check_command plink

## 检查需要的脚本文件是否存在且具有执行权限
check_command $PCA_cal $geno_dist $PCA_plot $LD_cor $corr_cal

## 基因型文件格式检查
check_plink ${bfile} ${nchr}

## 默认参数
maf=${maf:=0.05}
win=${win:="10000"}
inter=${inter:="99999"}
r2=${r2:="0"}
nchr=${nchr:="30"}

## 品种信息
if [[ ! ${breeds} ]]; then
  mapfile -t breeds_array < <(awk '{print $1}' ${bfile}.fam | sort | uniq)
else
  read -ra breeds_array <<<"$breeds"
fi

## 品种数目
np=${#breeds_array[@]}

## 提取品种基因型文件并分品种质控
for b in "${breeds_array[@]}"; do
  ## 提取指定品种基因型信息
  echo "${b}" >keep_fid_tmp.txt
  plink \
    --bfile ${bfile} \
    --keep-fam keep_fid_tmp.txt \
    --chr-set ${nchr} \
    --make-bed \
    --out ${b}

  ## 质控
  [[ ${wqc} ]] && \
    plink \
      --bfile ${b} \
      --maf ${maf} \
      --chr-set ${nchr} \
      --make-bed \
      --out ${b}q
done

## 筛选出在所有品种中均通过质控的标记位点
if [[ ${wqc} ]]; then
  awk '{print $2}' "${breeds_array[@]/%/q.bim}" | sort | uniq -c | awk -v n=${np} '$1==n {print $2}' >common.snp
  for b in "${breeds_array[@]}"; do
    plink \
      --bfile ${b}q \
      --extract common.snp \
      -chr-set ${nchr} \
      --make-bed \
      --out ${b}
  done
fi

## pca计算
$PCA_cal --pre_list "${breeds_array[*]}" --fids "${breeds_array[*]}" --fid

## pca作图
$PCA_plot \
  --eigv "$(printf '%s_' "${breeds_array[@]}")pca.txt" \
  --out "$(printf '%s_' "${breeds_array[@]}")pca"

## LD计算
for b in "${breeds_array[@]}"; do
  plink \
    --bfile ${b} \
    --freq \
    --r2 \
    --ld-window-kb ${win} \
    --ld-window ${inter} \
    --ld-window-r2 ${r2} \
    --chr-set ${nchr} \
    --out ${b}
done

## LD结果统计、作图
$LD_cor \
  --files "${breeds_array[*]}" \
  --popN "${breeds_array[*]}" \
  --bin1 50 \
  --breaks 1000 \
  --bin2 100 \
  --max 5000 \
  --out g_"$(printf '%s_' "${breeds_array[@]}")"_5Mb

## 品种对之间的基因频率和LD相关性大小
for t in frq ld; do
  :>${t}_cor.txt
  for bi in $(seq 0 $((np - 1))); do
    for bj in $(seq ${bi} $((np - 1))); do
      [[ ${bi} == "${bj}" ]] && continue
      cor=$($corr_cal --file1 ${breeds_array[${bi}]}.${t} --file2 ${breeds_array[${bj}]}.${t})
      echo "${breeds_array[${bi}]} ${breeds_array[${bj}]} ${cor}" >>${t}_cor.txt
    done
  done
  echo "correlation of ${t}:"
  cat ${t}_cor.txt
done

## 群体间的遗传距离
$geno_dist --bfile ${bfile} --out ${proj}/dist.summ

## 合并所有品种的基因型
if [[ ${wqc} ]]; then
  : >plink_merge_list.txt
  for b in "${breeds_array[@]}"; do
    [[ ${b} == "${breeds_array[0]}" ]] && continue
    echo "${b}" >>plink_merge_list.txt
  done
  plink \
    --bfile ${breeds_array[0]} \
    --merge-list plink_merge_list.txt \
    --maf 0.05 \
    --make-bed \
    --chr-set ${nchr} \
    --out ${bfile}
fi
