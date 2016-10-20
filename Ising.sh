#!/bin/bash

array=()

dimensions=("1D" "2D" "3D")
echo "Welcome. Please select the dimensions of the isotropic lattice."
echo -n "Type 1D, 2D or 3D: "
read var_one



if [[ " ${dimensions[*]} " != *" $var_one "* ]];
then
    echo "Option not allowed. Please try again"
    read var_right
else
    var_right=$var_one
fi
array+=($var_right)

echo "Please enter minimal temperature [J/k]"
echo -n "Tmin = " 
read t
array+=($t)

echo "Please enter maximal temperature [J/k]"
echo -n "Tmax = " 
read T
array+=($T)

echo "Please enter the values of lattice side size N, J and h."
echo -n "N = " 
read N
array+=($N)
echo -n "J = "
read J
array+=($J)
echo -n "h = "
read h
array+=($h)

echo "Please wait. Plots of observables are being produced."
./Ising.py  `echo ${array[@]}`

