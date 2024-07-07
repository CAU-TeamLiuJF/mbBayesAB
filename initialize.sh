#!/usr/bin/bash

########################################################################################################################
## Version: 1.2.0
## Author:  Liweining liwn@cau.edu.cn
## Date:    2023-07-07
## 
## Used to initialize the path in the project script, including
## 1.Replace the interpreter path in the R language script
## 2.Change project folder
## 
## Usage: ./initialize.sh
## 
## License:
##  This script is licensed under the GPL-3.0 License.
##  See https://www.gnu.org/licenses/gpl-3.0.en.html for details.
########################################################################################################################

## The path where this script is located
if [[ ${code} ]]; then
  [[ ! -d ${code} ]] && echo "${code} not exists! " && exit 5
else
  main_path=$(dirname "$(readlink -f "$0")")
fi

## Change to the path where the main script is located
cd ${main_path} || exit 5

## Replace the script path in the main script
sed -i -E "s|.*/(mbBayesAB)|\1${main_path}|" main.sh

## Load custom functions
func=${main_path}/shell/function.sh
[[ ! -s ${func} ]] && echo "Error: ${func} not found! " && exit 5
source ${func}

## Check if R language can be found in environment variables
check_command Rscript plink
if [[ $? -ne 0 ]]; then
  echo "Error: R language was not found in the current environment! "
  echo "Please install R language or load environment variables (such as module load R) before running this script"
  exit 1
fi

## R installation path
R_PATH=$(which Rscript)

## Replace the first line in the R script with the script interpreter mentioned above
find ./R -name "*.R" | while read -r file; do
  ## Check if the first line of the script is the script interpreter path
  if [[ $(head -n 1 "$file") =~ ^#!.*Rscript$ ]]; then
    ## Replace the script interpreter path in the first line with the path saved in the R_PATH variable
    sed -i "1s|.*|#!$R_PATH|" "$file"
  else
    # Insert a new script interpreter path line before the first line of the script
    sed -i "1i #!$R_PATH" "$file"
  fi
done

## Assign executable permissions to all scripts
chmod u+x main.sh
chmod u+x ./bin/* ./R/* ./shell/*

## Check if the required R package has been installed
pkg_check=${main_path}/R/package_required.R
if [[ -s $pkg_check ]]; then
  ${pkg_check}
else
  echo "Error: file ${pkg_check} not found! " && exit 5
fi

## Create the required folder
mkdir -p ${main_path}/log

[[ $? -eq 0 ]] && echo "Initialization completed."
