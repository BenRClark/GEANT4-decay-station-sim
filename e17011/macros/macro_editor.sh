#!/bin/bash


energies=("43" "87" "105" "123" "176" "247" "428" "444" "463" "478" "582" "591" "601" "612" "625" "636" "677" "692" "723" "757" "873" "893" "904" "996" "1004" "1119" "1246" "1274" "1596")


# Assign the filename


# Take the search string
read -p "Enter the search string: " search

# Take the replace string
read -p "Enter the replace string: " replace

for energy in ${!energies[@]}; do
filename="e17011_Eff_${energies[energy]}keV_no_Vis_face_depth.mac"
echo "$filename"
if [[ $search != "" && $replace != "" ]]; then
sed -i "s|$search|$replace|" $filename
fi
done
