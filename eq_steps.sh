#!/bin/bash

echo "Please enter a value for lattice side size at T = 0.3 J/K."
echo -n "N = "
read -a value
echo "Thank you. Please wait for the relaxation time to be produced."

./eq_steps.py echo $value

