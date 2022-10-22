#!/bin/sh
for file in *OUTPUT* iic2cicsc
do
    if test -r $file
    then
        rm $file
    fi
done
