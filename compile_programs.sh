#!/bin/sh
cd linalg
make 
cd ../subs
make
cd ../programs
sh ./compile_all.sh
