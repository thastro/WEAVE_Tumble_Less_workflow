#!/bin/bash

#
# Copyright (C) 2020 Cambridge Astronomical Survey Unit
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <https://www.gnu.org/licenses/>.
#

echo "======================================================"
echo " Simple FITS catalogue generation tool for WEAVE LIFU"
echo
echo " Please run it from a directory containing:"
echo
echo " <base_name>_params.csv"
echo " <base_name>_target_list.csv"
echo " <*>_CatalogueTemplate.fits"
echo "======================================================"

#-------------------------------------------------------------------------------
# Set directories
#-------------------------------------------------------------------------------

working_dir=$(pwd)
echo "- Current (working) directory: $working_dir"
#--------------------------------------------------

script_name="${BASH_SOURCE[0]}"
get_out()
{
  echo
  echo " Something went wrong :^( "
  echo
  if [ "$(basename "$script_name")" = "$(basename "$0")" ]; then # script was called directly
    exit
  else # script was called with 'source', and 'exit' would kill the shell
    return
  fi
}
script_fullpath="$(readlink -f "$script_name")"
script_dir=$(dirname $script_fullpath)
script_name=$(basename $script_fullpath)
echo "- Actual script location (resolving symlinks): $script_fullpath"
#--------------------------------------------------

ifu_dir=$(dirname $script_dir)
export PYTHONPATH=$ifu_dir
echo "- 'ifu' directory, assumed to contain the IFU workflow: $ifu_dir"
#--------------------------------------------------

workflow_dir="$ifu_dir/workflow"
echo "  IFU workflow: $workflow_dir"

driver_template="$workflow_dir/aux/ifu_driver_template.fits";
echo "  Driver template: $driver_template" 
if [ ! -f "$driver_template" ]; then get_out; return; fi
progtemp_file="$workflow_dir/aux/progtemp.dat"
echo "  PROGTEMP file: $progtemp_file"
if [ ! -f "$progtemp_file" ]; then
  echo "                 WARNING: You will need an Internet conection (and 'whitelisted' IP) to download the PROGTEMP file"
fi
obstemp_file="$workflow_dir/aux/obstemp.dat"
echo "  OBSTEMP file: $obstemp_file"
if [ ! -f "$obstemp_file" ]; then
  echo "                 WARNING: You will need an Internet conection (and 'whitelisted' IP) to download the OBSTEMP file"
fi
xml_template="$workflow_dir/aux/BlankXMLTemplate.xml"
echo "  XML template: $obstemp_file"
if [ ! -f "$xml_template" ]; then
  echo "                 WARNING: You will need an Internet conection (and 'whitelisted' IP) to download the XML template"
fi
#--------------------------------------------------

echo "- Looking for a directory named 'configure', containing/linking to the last version:"

# Working directory:
configure_dir="$working_dir/configure"
echo "  Trying $configure_dir..."

# Same as 'ifu' directory:
if [ ! -e "$configure_dir" ]; then
  configure_dir="$(dirname $ifu_dir)/configure"
  echo "  Trying $configure_dir..."
fi

if [ ! -e "$configure_dir" ]; then
  echo
  echo "ERROR: configure directory not found. Please create a link in the working directory:"
  echo "$working_dir"
  get_out; return
fi

echo "  Found!"
$configure_dir/configure --version
if [ $? -gt 0 ]; then
  get_out; return
fi

echo "======================================================"


#-------------------------------------------------------------------------------
# Set base_name (root for all input and output files).
# WEAVE standard:
# - '<TARGSRVY>_<trimester>' for Science Team catalogues
# - 'W<S/V><trimester>-<xxx>' for SV and Open Time
#-------------------------------------------------------------------------------

if test "$#" -eq 1; then
  base_name=$1
else
  files=($(ls *_params.csv))
  if [[ ${#files[@]} -eq 1 ]]; then
    base_name=$(basename "${files[0]}" _params.csv)
  elif [[ ${#files[@]} -eq 0 ]]; then
    echo "Sorry, no <base_name>_params.csv found in $working_dir"
    echo "======================================================"
    get_out; return
  else
    echo "ERROR: More than one <base_name>_params.csv found in"
    echo "       $working_dir"
    echo
    ls -1 *_params.csv
    echo
    echo "Please select one of them by passing an argument,"
    echo "and/or run all the following commands as a script:"
    echo
    for file in ${files[@]}; do
      filename=$(basename "$file" _params.csv)
      echo "$script_name \"$filename\""
    done
    echo
    echo "======================================================"
    get_out; return
  fi
fi

#-------------------------------------------------------------------------------
# Process CSV files
#-------------------------------------------------------------------------------

echo "- parameter file: ${base_name}_params.csv"
get_parameter()
{
  line=($(grep $1 ${base_name}_params.csv))
  param="${line[1]}"
}
get_parameter "PI_first_name"; PI_first_name="$param"
get_parameter "PI_surname"; PI_surname="$param"
get_parameter "catalogue_author"; catalogue_author="$param"
get_parameter "cc_report"; cc_report="$param"
get_parameter "GAIA_DR"; gaia_dr="$param"
get_parameter "TACALLOC"; tacalloc="$param"
#line=($(grep PI_first_name ${base_name}_params.csv)); PI_first_name=${line[1]}
#line=($(grep PI_surname ${base_name}_params.csv)); PI_surname=${line[1]}
#line=($(grep catalogue_author ${base_name}_params.csv)); catalogue_author=${line[1]}
#line=($(grep cc_report ${base_name}_params.csv)); cc_report=${line[1]}
#line=($(grep GAIA_DR ${base_name}_params.csv)); gaia_dr=${line[1]}
echo "                  survey PI: $PI_surname, $PI_first_name"
echo "                  catalogue author: $catalogue_author"
echo "                  Cc report: $cc_report"
echo "                  GAIA DR: $gaia_dr"
echo "                  TACALLOC: $tacalloc"

#--------------------------------------------------

if test -f "${base_name}_target_list.csv"; then
  echo "- target list:    ${base_name}_target_list.csv"
else
  echo "ERROR: '${base_name}_target_list.csv' not found in $working_dir"
  echo "======================================================"
  get_out; return
fi

#-------------------------------------------------------------------------------
# Set FITS catalog template
#-------------------------------------------------------------------------------

catalogue_template="${base_name}_CatalogueTemplate.fits"
if test -f "$catalogue_template"; then
  echo "- FITS catalogue template: $catalogue_template"
else
  files=($(ls *_CatalogueTemplate.fits))
  if test "${#files[@]}" -eq 1; then
    catalogue_template=${files[0]}
    echo "- FITS catalogue template: $catalogue_template"
  elif test "${#files[@]}" -eq 0; then
    echo "Sorry, no FITS catalogue templates found in $working_dir"
    echo "======================================================"
    get_out; return
  elif test "${#files[@]}" -gt 1; then
    echo "ERROR: More than one FITS catalogue template found in"
    echo "       $working_dir"
    echo
    ls -1 *_CatalogueTemplate.fits
    echo
    echo "Please rename (or link) one of them as $catalogue_template"
    echo "======================================================"
    get_out; return
  fi
fi
echo

# Set output directory

output_dir=$working_dir/$base_name
output_base=$output_dir/$base_name
echo "- output directory: $output_dir"
echo "======================================================"

# Run workflow

<<'# BLOCK-COMMENT'
# BLOCK-COMMENT
echo
echo "---------------------- Stage 1 ----------------------"
rm -rf $base_name
python "$workflow_dir/stage1/create_ifu_driver_cat_from_csv.py" $base_name
driver_catalogue="$output_base-ifu_driver_cat.fits"
if [ ! -f "$driver_catalogue" ]; then get_out; return; fi
python $workflow_dir/stage1/check_ifu_driver_cat.py $driver_catalogue --template $driver_template --progtemp_file $progtemp_file --obstemp_file $obstemp_file
#python stage1/plot_ifu_driver_cat.py $driver_catalogue

echo
echo "---------------------- Stage 2 ----------------------"
driver_catalogue="$output_base-ifu_driver_cat.fits"
if [ ! -f "$driver_catalogue" ]; then get_out; return; fi
python $workflow_dir/stage2/create_xml_files.py $driver_catalogue --outdir $output_dir --progtemp_file $progtemp_file --obstemp_file $obstemp_file --xml_template $xml_template

echo
echo "---------------------- Stage 3 ----------------------"
nfiles="$(ls $output_base*-t.xml | wc -l)"
if [[ $nfiles -eq 0 ]]; then get_out; return; fi
python $workflow_dir/stage3/add_guide_and_calib_stars.py $output_base*-t.xml  --outdir $output_dir

echo
echo "---------------------- Stage 4 ----------------------"

nfiles="$(ls $output_base*-tgc.xml | wc -l)"
if [[ $nfiles -eq 0 ]]; then get_out; return; fi
for (( counter=1; counter<=$nfiles; counter++ )) do
  num="$(printf "%02g" $counter)"
  file_base=$output_base-lifu_$num
  echo "Input XML: $output_base-tgc.xml"
  $configure_dir/lifu_configure --epoch 2020 -i $file_base-tgc.xml -o $file_base-tgcs-tmp.xml
  $configure_dir/dither -i $file_base-tgcs-tmp.xml -o $file_base-tgcs.xml
done

echo
echo "---------------------- Stage 5 ----------------------"
nfiles="$(ls $output_base*-tgcs.xml | wc -l)"
if [[ $nfiles -eq 0 ]]; then get_out; return; fi
python $workflow_dir/stage5/create_ifu_fits_cat.py $catalogue_template $output_base-*-tgcs.xml --out $output_base-ifu_from_xmls.fits --cat_nme1=$PI_first_name --cat_nme2=$PI_surname --gaia_dr $gaia_dr

echo
echo "---------------------- Stage 6 ----------------------"
if [ ! -f "$output_base-ifu_from_xmls.fits" ]; then get_out; return; fi
cp $output_base-ifu_from_xmls.fits $output_base.fits
python $workflow_dir/stage6/add_photometry_from_ps.py $output_base.fits
if [ ! -f "$output_base.fits" ]; then get_out; return; fi
python $workflow_dir/stage6/add_talloc_info.py $output_base.fits ${tacalloc}
echo "- Added TACALLOC code '$tacalloc' to catalogue"
