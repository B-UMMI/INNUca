#!/usr/bin/env bash

#script to update the mlst database

currentDir=$(pwd)

cd /NGStools/mlst/scripts/
./mlst-download_pub_mlst | bash
mv ../db/pubmlst ../db/pubmlst.old
mv ./pubmlst ../db/
./mlst-make_blast_db

cd $currentDir
