#!/bin/bash
#SBATCH --job-name=accur_GBLUP

########################################################################################################################
## 版本: 1.1.1
## 作者: 李伟宁 liwn@cau.edu.cn
## 日期: 2023-07-05
## 
## 统计各种情形下的方差组分结果
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
    --traits )   traits="$2";   shift 2 ;; ## 性状名称，如"DF DPM" ["/"]
    --rep )      rep="$2";      shift 2 ;; ## 第几次重复 ["/"]
    --dist )     dist="$2";     shift 2 ;; ## 加性遗传相关服从的分布 ["/"]
    --vg )       vg="$2";       shift 2 ;; ## 品种间加性方差大小差异 [""]
    --rgv )      rgv="$2";      shift 2 ;; ## 加性方差大小与遗传相关大小差异 [""]
    --ref )      ref="$2";      shift 2 ;; ## 加性方差大小与遗传相关大小差异 [""]
    --cor )      cor="$2";      shift 2 ;; ## 加性遗传相关大小 ["/"]
    --dirPre )   dirPre="$2";   shift 2 ;; ## ebv输出文件夹增加的前缀 [""]
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

## 默认参数
out=${out:=${proj}/varcomp_${today}.txt}
bins=${bins:="fix"}
dirPre=${dirPre:=""}
traits=${traits:="/"}
rep=${rep:="/"}
dist=${dist:="/"}
cor=${cor:="/"}
vg=${vg:="/"}
rgv=${rgv:="/"}
ref=${ref:="/"}

## 避免执行R脚本时的警告("ignoring environment value of R_HOME")
unset R_HOME

## DMU参数卡名称
DIR_full=phe_adj_PBLUP
DIR_within=within
DIR_blend=blend
DIR_union=union

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

##############  方差组分整理  ##############
###########################################
echo "模拟重复 其他处理 相关分布 相关大小 遗传方差 方差协方差 参考群比 模型 参考群 品种 性状 交叉重复 交叉折数 类型 值" >${out}
for r in "${reps_array[@]}"; do # b=${breeds_array[0]};rep=1;f=1;bin=fix
  for t in "${traits_array[@]}"; do # t=${traits_array[0]};r=${reps_array[0]};d=${dists_array[0]};c=${cors_array[0]}
    for d in "${dists_array[@]}"; do
      for c in "${cors_array[@]}"; do
        for v in "${vgs_array[@]}"; do # v=${vgs_array[0]};g=${rgvs_array[0]};e=${refs_array[0]}
          for g in "${rgvs_array[@]}"; do
            for e in "${refs_array[@]}"; do
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

              ## 品种内遗传力
              for b in "${breeds_array[@]}"; do
                ## 完整数据集遗传力
                if [[ -s ${path}/${b}/${DIR_full}.DIR ]]; then
                  h2_full=$(grep -A2 'Trait  correlation' ${path}/${b}/${DIR_full}.lst | tail -n 1 | awk '{print $2}')
                  echo "${r} / ${d} ${c} ${v} ${g} ${e} GBLUP ${b} ${b} ${t} / / h2 ${h2_full}" >>${out}
                fi

                ## 交叉验证参数
                rep=$(find ${path}/${b}/val1 -name "rep*" -type d | wc -l)
                fold=$(find ${path}/${b}/val* -name "rep1" -type d | wc -l)

                ## 每个子集方差组分
                for rep in $(seq 1 ${rep}); do # rep=1;f=1
                  for f in $(seq 1 ${fold}); do
                      ## GBLUP
                      lst=${path}/${b}/val${f}/rep${rep}/${DIR_within}.lst
                      if [[ -s ${lst} ]]; then
                        h2_within=$(grep -A2 'Trait  correlation' ${lst} | tail -n 1 | awk '{print $2}')
                        echo "${r} / ${d} ${c} ${v} ${g} ${e} w-GBLUP ${b} ${b} ${t} ${rep} ${f} h2 ${h2_within}" >>${out}
                      fi

                      ## Bayes
                      ebvf=${path}/${b}/val${f}/rep${rep}/EBV_fix_y1.txt
                      varf=${path}/${b}/val${f}/rep${rep}/var_fix.txt
                      if [[ -s ${ebvf} && -s ${varf} ]]; then
                        varg=$(sed '1d' ${ebvf} | awk '{sum+=$2; sumsq+=($2)^2} END {print (sumsq/NR-(sum/NR)^2)}')
                        vare=$(tail -n 1 ${varf})
                        h2_Bayes=$(echo "scale=4; $varg / ($varg + $vare)" | bc)
                        echo "${r} / ${d} ${c} ${v} ${g} ${e} w-fix ${b} ${b} ${t} ${rep} ${f} h2 ${h2_Bayes}" >>${out}
                      fi
                  done
                done
              done

              ## 单(多)性状GBLUP
              mapfile -t blends < <(find ${path} -name "blend*" -type d 2>/dev/null)
              for type in "${blends[@]}"; do # type=${blends[0]}
                paths=$(basename ${type})
                types=${paths/blend_/}
                IFS='_' read -r -a breeds_sub <<<"$types"
                n=${#breeds_sub[@]}

                ## types
                if [[ ${types} == "blend" ]]; then
                  types=""
                  comb=${breeds// /_}
                else
                  types="_${types}"
                  comb=$(IFS="_"; echo "${breeds_sub[*]}")
                fi

                ## 组合中每个品种
                for i in "${!breeds_sub[@]}"; do # i=0
                  ## 子集
                  for rep in $(seq 1 ${rep}); do # r=1;f=1
                    for f in $(seq 1 ${fold}); do
                      ## 单性状GBLUP联合评估
                      lst=${path}/blend${types}/val${f}/rep${rep}/${DIR_blend}.lst
                      if [[ -s ${lst} ]]; then
                        h2_blend=$(grep -A2 'Trait  correlation' ${lst} | tail -n 1 | awk '{print $2}')
                      # else
                      #   h2_blend=""
                      #   echo "${lst} not found! "
                      fi

                      ## 双性状联合评估
                      lst=${path}/union${types}/val${f}/rep${rep}/${DIR_union}.lst
                      if [[ -s ${lst} ]]; then
                        message="Correlation matrix for random"
                        h2_union=$(grep -A $((i + 2)) 'Trait  correlation' ${lst} | tail -n 1 | awk '{print $2}')

                        for j in $(seq $((i + 1)) $((n - 1))); do
                          rg_union=$(grep -A $((j + 2)) "${message}" ${lst} | tail -n 1 | awk -v col=$((i + 2)) '{print $col}')
                          echo "${r} / ${d} ${c} ${v} ${g} ${e} MT-GBUP ${comb} ${breeds_sub[i]}_${breeds_sub[j]} ${t} ${rep} ${f} rg ${rg_union}" >>${out}
                        done
                      # else
                      #   h2_union=""
                      #   echo "${lst} not found! "
                      fi

                      ## 写出到文件
                      {
                        echo "${r} / ${d} ${c} ${v} ${g} ${e} ST-GBLUP ${comb} ${breeds_sub[i]} ${t} ${rep} ${f} h2 ${h2_blend}"
                        echo "${r} / ${d} ${c} ${v} ${g} ${e} MT-GBLUP ${comb} ${breeds_sub[i]} ${t} ${rep} ${f} h2 ${h2_union}"
                      } >>${out}
                    done
                  done
                done
              done

              for bin in "${bins_array[@]}"; do
                accf=$(find ${path}/*ulti* -name "accur*${bin}*${b}.txt" 2>/dev/null)
                [[ ! ${accf} ]] && continue

                for fi in ${accf}; do # fi=$(echo ${accf} | awk '{print $1}')
                  ## 获取文件夹名称
                  dir=$(dirname ${fi})
                  dir=$(basename ${dir})

                  ## 分割文件夹名称
                  read -ra types <<<"${dir//multi_/ }"

                  ## 路径前缀处理
                  if [[ ${#types[@]} -gt 1 ]]; then
                    dirPre=${types[0]}
                    dirPre=${dirPre/%_/}
                    breeds=${types[1]}
                  else
                    dirPre=/
                    breeds=${types[0]}
                  fi

                  ## 品种
                  IFS='_' read -r -a breeds_sub <<<"$breeds"

                  for i in "${!breeds_sub[@]}"; do # i=0
                    ## 子集
                    for rep in $(seq 1 ${rep}); do # r=1;f=1
                      for f in $(seq 1 ${fold}); do
                        # for dirPre in /; do
                        multi_path=${path}/${dir}/val${f}/rep${rep}

                        ebvfi=$(find ${multi_path} -name "EBV_${bin}_y$((i + 1)).txt" 2>/dev/null)
                        varfi=$(find ${multi_path} -name "var_${bin}*.txt" 2>/dev/null)
                        [[ ! -s ${ebvfi} || ! -s ${varfi} ]] && continue

                        ## 遗传力
                        varg=$(sed '1d' ${ebvfi} | awk '{sum+=$2; sumsq+=($2)^2} END {print (sumsq/NR-(sum/NR)^2)}')
                        vare=$(tail -n 1 ${varfi} | awk -v col=$((i * n + i + 1)) '{print $col}')
                        h2_multi=$(echo "scale=4; $varg / ($varg + $vare)" | bc)
                        echo "${r} ${dirPre} ${d} ${c} ${v} ${g} ${e} ${bin} ${comb} ${breeds_sub[i]} ${t} ${rep} ${f} h2 ${h2_multi}" >>${out}

                        ## 遗传相关
                        for j in $(seq $((i + 1)) $((n - 1))); do
                          ebvfj=$(find ${multi_path} -name "EBV_${bin}_y$((j + 1)).txt" 2>/dev/null)
                          [[ ! -s ${ebvfj} ]] && continue
  
                          # 计算A文件中第2列减去均值的值
                          meanA=$(awk 'NR>1{sum+=$2}END{print sum/(NR-1)}' ${ebvfi})
                          awk -v meanA="$meanA" 'NR>1{print $2-meanA}' ${ebvfi} >A.tmp
                          # 计算B文件中第2列减去均值的值
                          meanB=$(awk 'NR>1{sum+=$2}END{print sum/(NR-1)}' ${ebvfj})
                          awk -v meanB="$meanB" 'NR>1{print $2-meanB}' ${ebvfj} >B.tmp
                          # 计算协方差
                          cov=$(paste A.tmp B.tmp | awk '{sum+=($1*$2)}END{print sum/(NR - 1)}')
                          var2=$(sed '1d' ${ebvfj} | awk '{sum+=$2; sumsq+=($2)^2} END {print (sumsq/NR-(sum/NR)^2)}')
  
                          ## 计算遗传相关
                          if [[ $(echo "$varg <= 0" | bc -l) -eq 1 || $(echo "$var2 <= 0" | bc -l) -eq 1 ]]; then
                            echo "varg=$varg var2=$var2"
                            rg_multi=0
                          else
                            rg_multi=$(echo "scale=4; $cov / sqrt($varg * $var2)" | bc | xargs printf "%.4f")
                          fi
  
                          echo "${r} ${dirPre} ${d} ${c} ${v} ${g} ${e} ${bin} ${comb} ${breeds_sub[i]}_${breeds_sub[j]} ${t} ${rep} ${f} rg ${rg_multi}" >>${out}
                        done
                      # done
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
  done
done

## 去除na值行
sed -i '/ $/d' ${out}
sed -i 's/ /\t/g' ${out}

[[ -s A.tmp ]] && rm A.tmp B.tmp
