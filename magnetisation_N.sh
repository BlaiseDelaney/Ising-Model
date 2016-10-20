#!/bin/bash

echo "Please enter values of N separated by a space."

read -a values
echo "Thank you. Please wait for the graph of magnetization vs temperature to be produced."

./magnetisation_N.py echo ${values[@]}

