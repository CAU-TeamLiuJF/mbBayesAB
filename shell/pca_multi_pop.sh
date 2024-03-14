#!/usr/bin/bash

########################################################################################################################
## 版本: 1.1.0
## 作者: 李伟宁 liwn@cau.edu.cn
## 日期: 2023-07-05
## 
## 对提供的任意多个群体进行PCA分析并作图
## 
## 使用: pca_multi_pop.sh --help
##
## 依赖软件/环境: 
##  1. R
##  2. plink/1.9
##  3. 其他R语言和Bash脚本
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


###################  参数处理  #####################
####################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## 参数名
TEMP=$(getopt -o h --long pre_list:,fids:,out:,nchr:,help,fid,plot \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
## 解析参数
while true; do
  case "$1" in
          --pre_list )   pre_list="$2";  shift 2 ;;  ## plink文件前缀，如"/public/home/popA /public/home/popB"
          --fids )       fids="$2"; 	   shift 2 ;;	 ## 品种（群体）标识符，如"popA popB"
          --nchr )       nchr="$2" ;     shift 2 ;;  ## 染色体数目 [30]
          --out )        out="$2" ;      shift 2 ;;  ## 输出文件名
          --fid )        fid=true;       shift   ;;  ## bfile中提供了品种（群体）标识符(即fid)
          --plot )       plot=true;      shift   ;;  ## 根据PCA结果作图
    -h | --help )        grep " shift " $0 && exit 1 ;;
    -- ) shift; break ;;
    * ) shift; break ;;
  esac
done

## 脚本所在文件夹
if [[ ${code} ]]; then
  [[ ! -d ${code} ]] && echo "${code} not exists! " && exit 5
else
  script_path=$(dirname "$(readlink -f "$0")")
  code=$(dirname "$script_path")
fi

## 脚本
pca_plot=${code}/R/PCA_plot.R
func=${code}/shell/function.sh

## 日志文件夹
logp=${code}/log

## 加载自定义函数
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${func}

## 检查需要的程序是否在环境变量中能检索到并且可执行
check_command plink

## 检查必要参数
[[ ! ${pre_list} ]] && echo "Open script $0 to view instructions" && exit 1

## 默认参数
nchr=${nchr:="30"}

## 只提供了一个文件，则根据fam文件中的fid进行分组
files_num=$(echo ${pre_list} | tr " " "\n" | wc -l)
if [[ ${files_num} -le 1 ]]; then
  unset fids
  if [[ ! ${fid} ]]; then
    echo "please provide the --fid option if the .ped(.fam) file contains fid"
    exit 1
  fi

  ## 检查文件是否存在
  check_plink "${pre_list}" ${nchr}

  ## plink基因型文件中的家系ID
  fid_uniq=$(awk '{print $1}' ${pre_list}.fam | sort | uniq)

  ## 提取不同品种的基因型信息
  for name in ${fid_uniq}; do
    ## 家系ID
    awk ${name} >${name}_fid.txt

    ## 提取指定家系ID的个体
    plink \
      --bfile ${pre_list} \
      --keep-fam ${name}_fid.txt \
      --chr-set ${nchr} \
      --make-bed --out ${name} >${logp}/plink_pca_keep_fam.log
    rm ${name}_fid.txt
  done

  pre_list=("${fid_uniq}")
  # Name_list="${pre_list[*]}"
fi

## 获取每个群体的标识符（family id）
# pre_list=("${pre_list[@]}")
IFS=' ' read -ra pre_list <<< "${pre_list[@]}"
index=0
for prefix in "${pre_list[@]}"; do # prefix=${pre_list[0]}
  ((index++))

  ## 检查文件是否存在
  check_plink "${prefix}" ${nchr}

  ## 群体id
  if [[ ${fids[*]} ]]; then
    fidi=$(echo "${fids[*]}" | cut -d ' ' -f ${index})
  else
    ## 提取第一行iid和fid，判断fid是否用iid表示或全为0，是的话则重新给个fid（如pop1）
    iid=$(head -n 1 ${prefix}.fam | awk '{print $2}')
    fidi=$(head -n 1 ${prefix}.fam | awk '{print $1}')
    if [[ ${fidi} == "${iid}" || ${fidi} == '0' ]]; then
      fidi=pop${index}
    fi
    # Name_list="${Name_list} ${fidi}"
  fi

  ## 输出iid和fid匹配表
  if [[ ${index} -eq 1 ]]; then
    awk '{print $2,"'${fidi}'"}' ${prefix}.fam >iid_fid_pca.txt
  else
    awk '{print $2,"'${fidi}'"}' ${prefix}.fam >>iid_fid_pca.txt

    ## 合并群体需要的文件名列表
    if [[ ${index} -ge 2 ]]; then
      echo ${prefix}.bed ${prefix}.bim ${prefix}.fam >>merge_list.txt
    else
      [[ -f merge_list.txt ]] && rm merge_list.txt
    fi
  fi
done

## 输出文件名
[[ ! ${out} ]] && out=$(echo "${pre_list[@]}" | tr " " "_")_pca.txt

## 合并所有群体plink文件
prefix1=$(echo "${pre_list[@]}" | cut -d ' ' -f 1)
num=$(echo "${pre_list[@]}" | wc -w)
plink \
  --bfile ${prefix1} \
  --merge-list merge_list.txt \
  --chr-set ${nchr} \
  --make-bed --out pop${num}_merge >${logp}/plink_pca_merge.log

## 存在多等位基因位点
if [[ $? -ne 0 ]]; then
  echo "please remove all offending variants (merged.missnp) in all populations"
  exit 2
fi

## pca计算
plink \
  --bfile pop${num}_merge \
  --chr-set ${nchr} \
  --pca 3 header \
  --out pop${num}_merge >${logp}/plink_pca.log

## 匹配群体标识符
if [[ ! ${fid} ]]; then
  sort -k 1n iid_fid_pca.txt -o iid_fid_pca.txt2
  sed '1d' pop${num}_merge.eigenvec | sort -k 2n > pop${num}_merge.eigenvec2
  join -1 2 -2 1 pop${num}_merge.eigenvec2 iid_fid_pca.txt2 >"${out}"
else
  cp pop${num}_merge.eigenvec "${out}" 
fi

## 作图
if [[ ${plot} ]]; then
  $pca_plot --eigv ${out}
fi

## 删除中间文件
rm pop${num}_merge.*
rm merge_list.*
rm iid_fid_pca.*
