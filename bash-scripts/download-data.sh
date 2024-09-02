#!/bin/bash

# remote computer
remote_name="office_upc"
remote_base_dir="/home/pol/work/materials-modelling/materials-simulations/chalcohalide-antiperovskites/silver/solid-solutions/VCA-grid/VCA_structures"

# local computer
local_dest_dir="/home/pol/work/upc/solid-solutions/data-results-relaxation/"

# Ask for the password
read -sp "Enter password for ${remote_user}@${remote_host}: " password

echo

# Loop through each directory from vca-001 to vca-121
for i in $(seq -f "vca-%03g" 1 121)
do
    echo "Directory ${i}"

    # Construct the remote directory path
    remote_dir="${remote_base_dir}/${i}/relaxation/*"

    
    # Use scp to copy files from the remote machine
    sshpass -p "$password" scp -r "${remote_name}:${remote_dir}" "${local_dest_dir}/${i}"

    # Optionally, add a delay between iterations to avoid overloading the remote server
    sleep 1
done
