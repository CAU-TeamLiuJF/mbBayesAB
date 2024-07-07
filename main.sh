#!/bin/usr/bash

########################################################################################################################
## Version:   1.2.0
## Author:    Liweining liwn@cau.edu.cn
## Orcid:     0000-0002-0578-3812
## Institute: College of Animal Science and Technology, China Agricul-tural University, Haidian, 100193, Beijing, China
## Date:      2024-07-07
##
## Function：
##  Used to reproduce the research results in the article "Multi-trait Bayesian models enhance the accuracy of genomic 
## prediction in multi-breed reference populations"
##
##
## Data sources ：
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

## Note：Initialize the script, modify some paths in the script
/public/home/liujf/liwn/code/GitHub/debug/mbBayesAB/initialize.sh

## Path of main script
code=/public/home/liujf/liwn/code/GitHub/debug/mbBayesAB
GP_cross=${code}/shell/GP_cross_validation.sh

## Add program path to environment variable
export PATH=${code}/bin:$PATH

######################## 数据分析 ########################
for source in Xie2021 Lee2019; do
  pro=${code}/data/${source}
  cd ${pro} || exit

  ## Phenotype and genotype files
  phef=${pro}/phenotype.txt
  bfile=${pro}/genotype

  ## breeds (populations) and trait names
  if [[ ${source} == "Xie2021" ]]; then
    all_eff="2 3 1"  ## intercept sex additive
    breeds=(YY LL)
    traits_all=(PFAI MS)
  else
    all_eff="2 1"    ## intercept additive
    breeds=(AAN LIM)
    traits_all=(WWT YWT)
  fi

  ## Analysis of population genetic structure
  $GP_cross \
    --proj ${pro} \
    --breeds "${breeds[*]}" \
    --bfile ${bfile} \
    --code ${code} \
    --type struct

  ## Calculate corrected phenotype
  for trait in "${traits_all[@]}"; do
    mkdir -p ${pro}/${trait}

    ## Parameters need to be adjusted according to the characteristics and requirements of Linux servers
    ## Note: If the Slurm workload manager is not installed in the system, please comment out the line where 'sbatch' is located in whole script
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

  ## Wait for the correction phenotype calculation to be completed
  while [[ $(wc -l 2>/dev/null <${pro}/${traits_all[-1]}/${breeds[-1]}/phe_adj_BLUP.SOL) -lt 10 ]]; do
    sleep 3
  done

  ## Random number seed for splitting reference/validation sets
  if [[ ! -s "${pro}/random.seed" ]]; then
    seed=$RANDOM
    echo ${seed} >"${pro}/random.seed"
  else
    seed=$(cat "${pro}/random.seed")
  fi

  ## Calculate the accuracy of GBLUP model in within-breed prediction
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

  ## wait for the GBLUP accuracy calculation to be completed
  while [[ ! -s ${pro}/${traits_all[-1]}/${breeds[-1]}/val1/rep1/pheno.txt ]]; do
    sleep 2
  done

  ## Calculate the accuracy of GBLUP model in multi-breed prediction
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
      sleep 5
    done
  done

  ## Calculate the accuracy of Bayes multi-breed prediction —————————— differen scale matrix
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
        sleep 5
    done
  done

  ## Calculate the accuracy of Bayes multi-breed prediction —————————— differen degree of freedom
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
        sleep 5
    done
  done

  ## Calculate the accuracy of Bayes multi-breed prediction —————————— HIW-WI prior
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
    sleep 5
  done

  ## Calculate the accuracy of Bayes multi-breed prediction —————————— HIW-IG prior
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
    sleep 5
  done

  ## Calculate the accuracy of Bayes multi-breed prediction —————————— Different genome block sizes
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
      sleep 5
    done
  done

  ## Statistical Accuracy and Variance Component Results
  for type in accur var; do
    $GP_cross \
      --type ${type} \
      --proj ${pro} \
      --breeds "${breeds[*]}" \
      --bin "fix" \
      --traits "${traits_all[*]}" \
      --code ${code}
    sleep 5
  done
done
