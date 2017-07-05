#!/usr/bin/env bash

#script to update the mlst database

/tmp/mlst/scripts/mlst-download_pub_mlst | bash
rm -r /tmp/mlst/db/pubmlst/
mv /tmp/mlst/scripts/pubmlst /tmp/mlst/db/
/tmp/mlst/scripts/mlst-make_blast_db

