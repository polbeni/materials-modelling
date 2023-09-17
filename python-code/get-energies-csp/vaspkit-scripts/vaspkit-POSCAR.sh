#!/bin/bash

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

file_path="$script_dir/input_vaspkit_POSCAR"
log_path="$script_dir/vaspkit_POSCAR.log"

cd $script_dir

vaspkit < $file_path > $log_path
