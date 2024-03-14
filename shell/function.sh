#!/bin/bash

########################################################################################################################
## 版本: 1.1.0
## 作者: 李伟宁 liwn@cau.edu.cn
## 日期: 2023-05-17
## 简介: 用于保存自定义的函数
## 
## 使用: ./function.sh
## 
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################

check_alphabet() {
  ## 检查文件中是否含有字符(DMU软件不接受字符类型的值)
  [[ ! -s ${1} ]] && echo "${1} not found! " && exit 1
  if [[ ${2} ]]; then
    NotNum=$(awk -vl=${2} '{print $l}' ${1} | grep -c "[a-zA-Z]")
  else
    NotNum=$(grep -v "[0-9]e" <${1} | grep -c "[a-zA-Z]")
  fi

  if [[ ${NotNum} -gt 0 ]]; then
    echo "Non numeric characters exist in ${1} file, please check! "
    exit 1
  fi
}

check_command() {
  # 检查脚本文件/软件是否存在并具有执行权限
  # 参数：
  # $@：需要检查的软件名称列表
  for cmd in "$@"; do
    ## 判断是否为文件且具有可执行权限
    if [[ -s ${cmd} ]]; then
      if [[ ! -x "$cmd" ]]; then
        echo "Error: \"$cmd\" does not have execute permission"
        return 1
      else
        return 0
      fi
    fi

    ## 判断是否为程序且具有可执行权限
    if ! command -v "$cmd" &> /dev/null; then
      echo "Error: \"$cmd\" not found (not installed or not in PATH)"
      return 1
    elif [[ ! -x $(command -v "$cmd") ]]; then
      echo "Error: \"$cmd\" does not have execute permission"
      return 1
    fi
  done
}

# check_command() {
#   # 检查脚本文件/软件是否存在并具有执行权限
#   # 参数：
#   # $@：需要检查的软件名称列表
#   for cmd in "$@"; do
#     ## 判断是否为文件且具有可执行权限
#     if [[ -x "${cmd}" ]]; then
#       return 0
#     fi

#     ## 判断是否为程序且具有可执行权限
#     if ! type -a "$cmd" &> /dev/null; then
#       echo "Error: \"$cmd\" not found (not installed or not in PATH)"
#       return 1
#     elif [[ ! -x $(type -p "$cmd") ]]; then
#       echo "Error: \"$cmd\" does not have execute permission"
#       return 1
#     fi
#   done
# }

check_plink() {
  ## 检查指定前缀的plink二进制文件是否存在
  ## 若存在map、ped则转化为二进制
  ## 参数1：plink文件名前缀，如/path/to/plink_prefix
  ## 参数2：染色体数目，默认为30

  ## 检查plink软件是否可用
  check_command plink

  ## 染色体数设置
  if [[ $2 ]]; then
    nchr=$2
  else
    nchr=30
  fi

  ## 文件格式
  for prefix in $1; do
    if [[ ! -f ${prefix}.fam ]]; then
      if [[ -f ${prefix}.map ]]; then
        plink \
          --file ${prefix} \
          --chr-set ${nchr} \
          --make-bed \
          --out ${prefix}
      else
        echo "Error: ${prefix}.bed not found"
        # exit 1
      fi
    fi
  done
}

merge_plink() {
  ## 合并plink文件
  ## 参数1：plink文件列表，以空格分隔，如"/path/to/plinkA_prefix /path/to/plinkB_prefix"，两侧需有双引号
  ## 参数2：输出文件名前缀，如/path/to/plink_merge_prefix
  ## 参数3：染色体数目，默认为30
  ## 依赖：plink

  ## 要合并的文件名列表
  local plink_list
  IF=" " read -r -a plink_list <<<"$1"

  ## 染色体数设置
  if [[ $3 ]]; then
    nchr=$3
  else
    nchr=30
  fi

  ## 检查参数个数
  nfiles=${#plink_list[@]}
  if [[ ${nfiles} -lt 2 ]]; then
    echo "Error: at least two plink files are needed to merge"
    return 1
  fi

  ## 临时文件，用于存放需要合并的群体的文件名前缀
  local seed_tmp=$RANDOM
  touch merge_${seed_tmp}.txt

  ## 检查各群体的plink文件是否存在
  check_plink "${plink_list[0]}" ${nchr}
  for i in $(seq 1 $((nfiles - 1))); do
    check_plink "${plink_list[${i}]}" ${nchr}
    echo ${plink_list[${i}]} >>merge_${seed_tmp}.txt
  done

  ## 检查输出文件前缀对应的plink文件是否已存在
  if [[ -f $2.bed ]]; then
    echo "warn: $2.bed already existd and will be overwrited! "
  fi

  ## 合并文件
  [[ -s merge_${seed_tmp}.txt ]] && \
    plink \
      --bfile ${plink_list[0]} \
      --chr-set ${nchr} \
      --merge-list merge_${seed_tmp}.txt \
      --make-bed \
      --out $2

  rm merge_${seed_tmp}.txt
}
