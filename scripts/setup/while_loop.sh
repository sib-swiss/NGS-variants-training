#!/usr/bin/env bash

cat animals.txt | while read animal behaviour leg_number
do
    echo "The $animal $behaviour, and has $leg_number legs" 
done
