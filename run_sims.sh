#!/bin/bash

for value in 400 600 800 1000 1200; do
    Rscript simulation1.r "$value" 30 &
    Rscript simulation2.r "$value" 30 &
    Rscript simulation1.r "$value" 80 &
    Rscript simulation2.r "$value" 80 &
done

