#!/usr/bin/env bash


qsub  -l m_mem_free=50G -pe threads 4 -cwd -V homer_tag.sh $1 $2

