#!/bin/usr/bash

########################################################################################################################
## 版  本: 1.0.0
## 作  者: 李伟宁 liwn@cau.edu.cn
## Orcid: 0000-0002-0578-3812
## 单  位: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## 日  期: 2024-03-14
##
## 功能：
##  用于复现文章"Multi-trait Bayesian models enhance the accuracy of genomic prediction in 
##  multi-breed reference populations"文章中的研究结果
##
## 数据来源：
##  Xie, L., J. Qin, L. Rao, X. Tang and D. Cui et al., 2021 Accurate prediction and genome-wide 
##  association analysis of digital intramuscular fat content in longissimus muscle of pigs. Animal Genetics 
##  52: 633-644. https://doi.org/10.1111/age.13121
##
##  [Lee, J.](https://doi.org/10.1111/age.12846), Kim, J. M., & Garrick, D. J. (2019). Increasing the accuracy
##  of genomic prediction in pure-bred Limousin beef cattle by including cross-bred Limousin data and accounting 
##  for an F94L variant in MSTN. Animal genetics, 50(6), 621–633. 
##
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################

## ***注意***：使用前需要以绝对路径形式运行根目录下的initialize.sh脚本，以修改脚本中相应的路径，如：
## $ /public/home/liwn/GitHub/mbBayesAB/initialize.sh

## 主脚本路径
code=/public/home/liujf/liwn/code/GitHub/mbBayesAB
GP_cross=${code}/shell/GP_cross_validation.sh

## 将程序路径加到环境变量中
export PATH=${code}/bin:$PATH

######################## 数据分析 ########################
for source in Xie2021 Lee2019; do
  pro=${code}/data/${source}
  cd ${pro} || exit

  ## 表型和基因型文件
  phef=${pro}/phenotype.txt
  bfile=${pro}/genotype

  ## 品种和性状名称
  if [[ ${source} == "Xie2021" ]]; then
    all_eff="2 3 1"  ## intercept sex additive
    breeds=(YY LL)
    traits_all=(PFAI MS)
  else
    all_eff="2 1"    ## intercept additive
    breeds=(AAN LIM)
    traits_all=(WWT YWT)
  fi

  ## 群体结构
  $GP_cross \
    --proj ${pro} \
    --breeds "${breeds[*]}" \
    --bfile ${bfile} \
    --code ${code} \
    --type struct

  ## 计算校正表型
  for trait in "${traits_all[@]}"; do
    mkdir -p ${pro}/${trait}

    ## 需根据Linux服务器特性和需求调整参数
    ##  注意：如果在个人电脑而不是装有slurm作业管理系统的服务器上运行，请注释掉sbatch所在行，下同
    sbatch -c2 --mem=4G \
    $GP_cross \
      --proj ${pro} \
      --breeds "${breeds[*]}" \
      --traits "${traits_all[*]}" \
      --trait "${trait}" \
      --bfile ${bfile} \
      --phef ${phef} \
      --all_eff "${all_eff}" \
      --code ${code} \
      --type adj \
      --append
    sleep 5
  done

  ## 等待校正表型计算结束
  while [[ $(wc -l 2>/dev/null <${pro}/${traits_all[-1]}/${breeds[-1]}/phe_adj_BLUP.SOL) -lt 10 ]]; do
    sleep 3
  done

  ## 用于划分参考/验证集的随机数种子
  if [[ ! -s "${pro}/random.seed" ]]; then
    seed=$RANDOM
    echo ${seed} >"${pro}/random.seed"
  else
    seed=$(cat "${pro}/random.seed")
  fi

  ## 计算品种内评估准确性GBLUP
  for trait in "${traits_all[@]}"; do
    [[ -s "${pro}/${trait}/${breeds[-1]}/accur_GBLUP.txt" ]] && continue
    sbatch -c50 --time=2:30:00 --mem=150G \
    $GP_cross \
      --proj ${pro} \
      --breeds "${breeds[*]}" \
      --traits "${traits_all[*]}" \
      --trait "${trait}" \
      --bfile ${bfile} \
      --phef ${phef} \
      --all_eff "${all_eff}" \
      --seed ${seed} \
      --code ${code} \
      --rep 10 \
      --fold 5 \
      --thread 50 \
      --type within
    sleep 5
  done

  ## 等待
  while [[ ! -s ${pro}/${traits_all[-1]}/${breeds[-1]}/val1/rep1/pheno.txt ]]; do
    sleep 2
  done

  ## 计算多品种GBLUP评估准确性
  for trait in "${traits_all[@]}"; do
    for type in blend union; do
      sbatch -c50 --time=4:00:00 --mem=150G \
      $GP_cross \
        --type ${type} \
        --proj ${pro} \
        --breeds "${breeds[*]}" \
        --traits "${traits_all[*]}" \
        --trait "${trait}" \
        --code ${code} \
        --dense \
        --phef ${phef} \
        --seed ${seed} \
        --thread 50 \
        --suffix
    done
  done

  ## 计算多品种多性状Bayes评估准确性 —————————— 不同的尺度参数
  for trait in "${traits_all[@]}"; do
    for So in iden3 iden iden0.01 phe; do
      sbatch -c50 --exclude=cnode1002 --time=2:30:00 --mem=150G \
      $GP_cross \
        --type multi \
        --proj ${pro} \
        --breeds "${breeds[*]}" \
        --traits "${traits_all[*]}" \
        --trait "${trait}" \
        --code ${code} \
        --seed ${seed} \
        --phef ${phef} \
        --thread 50 \
        --VaSori ${So} \
        --dirPre ${So}_ \
        --suffix
        sleep 10
    done
  done

  ## 计算多品种多性状Bayes评估准确性 —————————— 不同的自由度参数
  for trait in "${traits_all[@]}"; do
    for df in 2 3 4 5 6 100; do
      sbatch -c50 --exclude=cnode1002 --time=2:30:00 --mem=150G \
      $GP_cross \
        --type multi \
        --proj ${pro} \
        --breeds "${breeds[*]}" \
        --traits "${traits_all[*]}" \
        --trait "${trait}" \
        --code ${code} \
        --seed ${seed} \
        --phef ${phef} \
        --thread 50 \
        --dfva ${df} \
        --dirPre df${df}_ \
        --suffix
        sleep 10
    done
  done

  ## 计算多品种多性状Bayes评估准确性 —————————— HIW-WI先验
  for trait in "${traits_all[@]}"; do
    sbatch -c50 --exclude=cnode1002 --time=2:30:00 --mem=150G \
    $GP_cross \
      --type multi \
      --proj ${pro} \
      --breeds "${breeds[*]}" \
      --traits "${traits_all[*]}" \
      --trait "${trait}" \
      --code ${code} \
      --seed ${seed} \
      --phef ${phef} \
      --thread 50 \
      --VaSd wishart \
      --dfva 3 \
      --dirPre wishart_df3_ \
      --suffix
  done

  ## 计算多品种多性状Bayes评估准确性 —————————— HIW-IG先验
  for trait in "${traits_all[@]}"; do
    sbatch -c50 --exclude=cnode1002 --time=2:30:00 --mem=150G \
    $GP_cross \
      --type multi \
      --proj ${pro} \
      --breeds "${breeds[*]}" \
      --traits "${traits_all[*]}" \
      --trait "${trait}" \
      --code ${code} \
      --seed ${seed} \
      --phef ${phef} \
      --thread 50 \
      --VaSd invgamma \
      --dfva 2 \
      --dirPre invgamma_df2_ \
      --suffix
  done

  ## 计算多品种多性状Bayes评估准确性 —————————— 基因组区块大小
  for trait in "${traits_all[@]}"; do
    for win in 1 25 50 100 200 whole; do
      sbatch -c50 --exclude=cnode1002 --time=2:30:00 --mem=150G \
      $GP_cross \
        --type multi \
        --proj ${pro} \
        --breeds "${breeds[*]}" \
        --traits "${traits_all[*]}" \
        --trait "${trait}" \
        --code ${code} \
        --seed ${seed} \
        --phef ${phef} \
        --thread 50 \
        --nsnp_win ${win} \
        --dirPre win${win}_ \
        --suffix
    done
  done

  ## 统计准确性和方差组分结果
  for type in accur var; do
    $GP_cross \
      --type ${type} \
      --proj ${pro} \
      --breeds "${breeds[*]}" \
      --bin "fix" \
      --traits "${traits_all[*]}" \
      --code ${code}
  done
done
