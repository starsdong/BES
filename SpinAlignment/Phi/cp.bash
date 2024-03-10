#!/bin/bash

SOURCE=/Users/starsdong/work/work/BES/SpinAlignment/Phi


for FILE in $(cat ${SOURCE}/test.list); do
  echo $FILE
  cp $SOURCE/$FILE ./$FILE
done
