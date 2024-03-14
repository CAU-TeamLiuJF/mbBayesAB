#!/bin/bash
#SBATCH --job-name=accuracy

########################################################################################################################
## 版本: 1.1.1
## 作者: 李伟宁 liwn@cau.edu.cn
## 日期: 2023-07-05
## 
## 统计各种情形下的交叉验证准确性结果
## 
## 使用: ./varcomp_summary.sh --help
## 
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################


###################  参数处理  #####################
####################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## 参数名
TEMP=$(getopt -o h --long code:,proj:,breeds:,rep:,dist:,vg:,rgv:,ref:,cor:,traits:,bin:,dirPre:,out:,help \
              -n 'javawrap' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
## 解析参数
while true; do
  case "$1" in
    --proj )     proj="$2";     shift 2 ;; ## 项目目录 [必要参数]
    --breeds )   breeds="$2";   shift 2 ;; ## 群体/品种标识符，如'YY DD' [必要参数]
    --traits )   traits="$2";   shift 2 ;; ## 性状名称，如"DF DPM" [""]
    --rep )      rep="$2";      shift 2 ;; ## 第几次重复 [""]
    --dist )     dist="$2";     shift 2 ;; ## 加性遗传相关服从的分布 [""]
    --vg )       vg="$2";       shift 2 ;; ## 品种间加性方差大小差异 [""]
    --rgv )      rgv="$2";      shift 2 ;; ## 加性方差大小与遗传相关大小差异 [""]
    --ref )      ref="$2";      shift 2 ;; ## 加性方差大小与遗传相关大小差异 [""]
    --cor )      cor="$2";      shift 2 ;; ## 加性遗传相关大小 [""]
    --dirPre )   dirPre="$2";   shift 2 ;; ## ebv文件夹额外前缀 [""]
    --bin )      bins="$2";     shift 2 ;; ## 多品种评估时区间划分方法，fix/frq/ld ["fix"]
    --code )     code="$2";     shift 2 ;; ## 脚本文件所在目录，如/BIGDATA2/cau_jfliu_2/liwn/code [NULL]
    --out )      out="$2";      shift 2 ;; ## 准确性输出文件名 [accuracy_$date.txt]
  -h | --help)    grep ";; ##" $0 | grep -v help && exit 1 ;;
  -- ) shift; break ;;
  * ) shift; break ;;
  esac
done

## 检查必要参数是否提供
if [[ ! -d ${proj} ]]; then
  echo "${proj} not found! "
  exit 1
elif [[ ! ${breeds} ]]; then
  echo "para --breeds is reduired! "
  exit 1
fi

## 日期
today=$(date +%Y%m%d)

## 避免执行R脚本时的警告("ignoring environment value of R_HOME")
unset R_HOME

## 脚本所在文件夹
if [[ ${code} ]]; then
  [[ ! -d ${code} ]] && echo "${code} not exists! " && exit 5
else
  script_path=$(dirname "$(readlink -f "$0")")
  code=$(dirname "$script_path")
fi

## 默认参数
out=${out:=${proj}/accuracy_${today}.txt}
bins=${bins:="fix"}
dirPre=${dirPre:=""}
traits=${traits:="/"}
rep=${rep:="/"}
dist=${dist:="/"}
cor=${cor:="/"}
vg=${vg:="/"}
rgv=${rgv:="/"}
ref=${ref:="/"}

## 解析参数
read -ra breeds_array <<<"$breeds"
read -ra bins_array <<<"$bins"
read -ra traits_array <<<"$traits"
read -ra reps_array <<<"$rep"
read -ra dists_array <<<"$dist"
read -ra cors_array <<<"$cor"
read -ra vgs_array <<<"$vg"
read -ra rgvs_array <<<"$rgv"
read -ra refs_array <<<"$ref"

##############  准确性结果统计(各种组合)  ##########
################################################
echo "模拟重复 其他处理 相关分布 相关大小 遗传方差 方差协方差 参考群比 模型 参考群 品种 性状 重复 交叉折数 准确性 无偏性 秩相关 验证群大小" >${out}
for t in "${traits_array[@]}"; do # t=${traits_array[0]};b=${breeds_array[0]}
  for r in "${reps_array[@]}"; do # r=${reps_array[0]};d=${dists_array[0]};c=${cors_array[0]}
    for d in "${dists_array[@]}"; do 
      for c in "${cors_array[@]}"; do
        for v in "${vgs_array[@]}"; do # v=${vgs_array[0]};g=${rgvs_array[0]};e=${refs_array[0]}
          for g in "${rgvs_array[@]}"; do
            for e in "${refs_array[@]}"; do
              for b in "${breeds_array[@]}"; do
                path=${proj}
                ## 模拟情形下的路径设置
                [[ ${r} != "/" ]] && path=${path}/rep${r}
                [[ ${t} != "/" ]] && path=${path}/${t}
                [[ ${d} != "/" ]] && path=${path}/${d}
                [[ ${c} != "/" ]] && path=${path}/cor${c}
                [[ ${v} != "/" ]] && path=${path}/v${v}
                [[ ${g} != "/" ]] && path=${path}/g${g}
                [[ ${e} != "/" ]] && path=${path}/e${e}

                ## 处理路径中存在多个斜杠的情况
                path=$(echo "$path" | sed 's#/\{2,\}#/#g; s#/$##')

                ## 判断文件夹是否存在
                [[ ! -d ${path} ]] && continue

                ## 交叉验证参数
                [[ ! -d ${path}/${b}/val1 ]] && continue
                rep=$(find ${path}/${b}/val1 -name "rep*" -type d | wc -l)
                fold=$(find ${path}/${b}/val* -name "rep1" -type d | wc -l)

                ## 群体内
                wf=${path}/${b}/accur_GBLUP.txt
                # bf=${path}/${b}/accur_BayesAS.txt
                # if [[ -s ${wf} ]]; then
                #   sed '1d' ${wf} | awk '{print "within","'${b}'","'${b}'","'${t}'", $0}' >>${out}
                # # else
                #   # echo "${accf} not found! "
                # fi

                ## 群体合并单性状
                accf=$(find ${path} -path "*/blen*/*" -name "accur_GBLUP*${b}*txt" 2>/dev/null)
                if [[ ${accf} ]]; then
                  # accf=${path}/blend/accur_GBLUP_${b}.txt
                  for f in ${accf}; do # f=${accf[0]}
                    if [[ -s ${f} ]]; then
                      ## 文件夹名
                      type=$(dirname ${f})
                      type=$(basename ${type})
                      type=${type/blend_/}
                      {
                        awk '{print "'${r}'","b-GBLUP","'${d}'","'${c}'","'${v}'","'${g}'","'${e}'","b-GBLUP","'${type}'","'${b}'","'${t}'",$0}' ${f}
                        # [[ -s ${wf} ]] && \
                        #   awk '{print "'${r}'","/","'${d}'","'${c}'","w-GBLUP","'${type}'","'${b}'","'${t}'", $0}' ${wf}
                        # [[ -s ${bf} ]] && \
                        #   awk '{print "'${r}'","/","'${d}'","'${c}'","w-BayesAS","'${type}'","'${b}'","'${t}'", $0}' ${bf}
                      } >>${out}
                    # else
                      # echo "${accf} not found! "
                    fi
                  done
                fi

                ## 群体合并多性状GBLUP
                accf=$(find ${path} -path "*/unio*/*" -name "accur_GBLUP*${b}*txt" 2>/dev/null)
                # accf=$(find ${path}/unio* -name "accur_*${b}*txt" 2>/dev/null)
                if [[ ${accf} ]]; then
                  for f in ${accf}; do # f=${accf[0]}
                    if [[ -s ${f} ]]; then
                      ## 文件夹名
                      type=$(dirname ${f})
                      type=$(basename ${type})
                      type=${type/union_/}
                      {
                        awk '{print "'${r}'","u-GBLUP","'${d}'","'${c}'","'${v}'","'${g}'","'${e}'","u-GBLUP","'${type}'","'${b}'","'${t}'",$0}' ${f}
                        # [[ -s ${wf} ]] && \
                        #   awk '{print "'${r}'","/","'${d}'","'${c}'","w-GBLUP","'${type}'","'${b}'","'${t}'", $0}' ${wf}
                        # [[ -s ${bf} ]] && \
                        #   awk '{print "'${r}'","/","'${d}'","'${c}'","w-BayesAS","'${type}'","'${b}'","'${t}'", $0}' ${bf}
                      } >>${out}
                    # else
                    #   echo "${accf} not found! "
                    fi
                  done
                fi

                ## MT-bayesAS
                for bin in "${bins_array[@]}"; do # bin=${bins_array[0]}
                  accf=$(find ${path}/*mult* -name "accur*${bin}*${b}.txt" 2>/dev/null)
                  [[ ! ${accf} ]] && continue

                  for f in ${accf}; do # f=$(echo ${accf} | awk '{print $1}')
                    ## 获取文件夹名称
                    type=$(dirname ${f})
                    type=$(basename ${type})

                    ## 分割文件夹名称
                    read -ra types <<<"${type//multi_/ }"

                    ## 路径前缀处理
                    if [[ ${#types[@]} -gt 1 ]]; then
                      dirPre=${types[0]}
                      dirPre=${dirPre/%_/}
                      breeds=${types[1]}
                    else
                      dirPre=/
                      breeds=${types[0]}
                    fi

                    bf=${path}/${b}/accur_BayesAS_${bin}_${type}.txt
                    # echo "bf=$bf" && exit
                    if [[ -s ${f} ]]; then
                      {
                        awk '{print "'${r}'","'${dirPre}'","'${d}'","'${c}'","'${v}'","'${g}'","'${e}'","'mbBayesAB-${bin}'","'${breeds}'","'${b}'","'${t}'",$0}' ${f}
                        [[ -s ${wf} ]] && \
                          awk '{print "'${r}'","w-GBLUP","'${d}'","'${c}'","'${v}'","'${g}'","'${e}'","w-GBLUP","'${breeds}'","'${b}'","'${t}'", $0}' ${wf}
                        [[ -s ${bf} ]] && \
                          awk '{print "'${r}'","w-BayesAS","'${d}'","'${c}'","'${v}'","'${g}'","'${e}'","w-BayesAS","'${breeds}'","'${b}'","'${t}'", $0}' ${bf}
                      } >>${out}
                    # else
                    #   echo "${f} not found! "
                    fi
                  done
                done
              done
            done
          done
        done
      done
    done
  done
done

## 去掉首位行
sed -i '/rep/d' ${out}
sed -i '/mean/d' ${out}
## 替换空格
sed -i 's/ /\t/g' ${out}
