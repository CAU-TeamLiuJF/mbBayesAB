#!/usr/bin/bash


########################################################################################################################
## 版本: 1.1.0
## 作者: 李伟宁 liwn@cau.edu.cn
## 日期: 2023-07-05
## 
## 对提供的任意多个群体进行平均遗传距离计算
## 
## 使用: distance_multi_pop.sh --help
## 
## 依赖软件/环境:
## R
## plink (1.90)
## 
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


###################  参数处理  #####################
####################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## 参数名
TEMP=$(getopt -o h --long bfile:,out:,help \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
## 解析参数
while true; do
  case "$1" in
    --bfile )   bfile="$2";  shift 2 ;;  ## plink文件前缀，如"/public/home/popA /public/home/popB"
    --nchr )    nchr="$2" ;  shift 2 ;;  ## 染色体数目 [30]
    --out )     out="$2" ;   shift 2 ;;  ## 输出文件名
    -h | --help )     grep " shift " $0 && exit 1 ;;
    -- ) shift; break ;;
    * ) shift; break ;;
  esac
done

## 默认参数
out=${out:="distance"}
nchr=${nchr:="30"}

## 脚本所在文件夹
if [[ ${code} ]]; then
  [[ ! -d ${code} ]] && echo "${code} not exists! " && exit 5
else
  script_path=$(dirname "$(readlink -f "$0")")
  code=$(dirname "$script_path")
fi

## 脚本
dist=${code}/R/mean_distance.R
func=${code}/shell/function.sh

## 载入自定义函数
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${func}

## 检查plink二进制文件是否存在
check_plink "${bfile}" ${nchr}

## 工作文件夹
dir=$(dirname ${bfile})
cd ${dir} || exit

## 个体按照FID排序
plink --bfile ${bfile} --chr-set ${nchr} --indiv-sort natural --make-bed --out distance_tmp

## 计算距离矩阵
plink --bfile distance_tmp --chr-set ${nchr} --distance square0 1-ibs --out distance_tmp

## 计算矩阵中品种配对对应的块的均值
$dist --prefix distance_tmp --out ${out}

## 删除中间文件
rm distance_tmp.*
