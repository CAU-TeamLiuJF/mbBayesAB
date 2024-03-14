#!/usr/bin/bash

########################################################################################################################
## 版本: 1.1.0
## 作者: 李伟宁 liwn@cau.edu.cn
## 日期: 2023-07-05
## 
## 用于初始化该项目脚本中的路径，主要包括
## 1.更换R语言脚本中的解释器路径
## 2.更换项目文件夹
## 
## 使用: ./initialize.sh
## 
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################

## 此脚本所在路径
if [[ ${code} ]]; then
  [[ ! -d ${code} ]] && echo "${code} not exists! " && exit 5
else
  main_path=$(dirname "$(readlink -f "$0")")
fi

## 切换到主脚本所在路径
cd ${main_path} || exit 5

## 将主脚本中的脚本路径进行替换
sed -i -E "s|(code=).*/(GP-cross-validation)|\1${main_path}|" main.sh
sed -i -E "s|(output=).*/(GP-cross-validation)|\1${main_path}|" ./shell/GP_cross_validation.sh

## 加载自定义函数
func=${main_path}/shell/function.sh
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${func}

## 检查R语言是否能在环境变量中找到
check_command Rscript
if [[ $? -ne 0 ]]; then
  echo "未在当前环境下发现R语言，请先安装R语言或加载环境变量（如module load R）后再运行本脚本"
  exit 1
fi

## R安装路径
R_PATH=$(which Rscript)

## 将R脚本中的第一行替换为脚本上述的脚本解释器
find ./R -name "*.R" | while read -r file; do
  # 检查脚本的第一行是否为脚本解释器路径
  if [[ $(head -n 1 "$file") =~ ^#!.*Rscript$ ]]; then
    # 替换第一行中的脚本解释器路径为 R_PATH 变量中保存的路径
    sed -i "1s|.*|#!$R_PATH|" "$file"
  else
    # 在脚本的第一行之前插入新的脚本解释器路径行
    sed -i "1i #!$R_PATH" "$file"
  fi
done

## 将所有的脚本赋予可执行权限
chmod u+x main.sh
chmod u+x ./bin/* ./R/* ./shell/*

## 检查需要的R包是否已安装
pkg_check=${main_path}/R/package_required.R
if [[ -s $pkg_check ]]; then
  ${pkg_check}
else
  echo "Error: file ${pkg_check} not found! " && exit 5
fi

## 创建需要的文件夹
mkdir -p ${main_path}/log

[[ $? -eq 0 ]] && echo "Initialization completed."
