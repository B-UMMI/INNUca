#!/usr/bin/env bash

#script to update the mlst database

/NGStools/mlst/scripts/mlst-download_pub_mlst | bash
rm -r /NGStools/mlst/db/pubmlst/
mv /NGStools/mlst/scripts/pubmlst /NGStools/mlst/db/
/NGStools/mlst/scripts/mlst-make_blast_db

