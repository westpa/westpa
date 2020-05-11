#!/bin/bash

gmx distance \
  -f minimization.gro\
  -s minimization.tpr\
  -pbc \
  -select 1\
  -oall dist.xvg

cat dist.xvg | tail -n 1 | awk '{print $2}'


