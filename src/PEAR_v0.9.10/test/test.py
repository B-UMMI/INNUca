#! /usr/bin/env python
import sys
import os
import filecmp
import subprocess
import shutil

def unit_test(f, r, o, fcompare):
    subprocess.call(["../src/pear", "-f", f, "-r", r, "-o", o, "-c", "0"], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
    correct = filecmp.cmp(o+".assembled.fastq", fcompare)
    os.remove(o+".assembled.fastq")
    os.remove(o+".discarded.fastq")
    os.remove(o+".unassembled.forward.fastq")
    os.remove(o+".unassembled.reverse.fastq")
    if correct:
        return "passed"
    else:
        return "failed"

def gen(f, r, o, fcompare):
    subprocess.call(["../src/pear", "-f", f, "-r", r, "-o", o, "-c", "0"], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
    shutil.copyfile(o+".assembled.fastq", fcompare)
    os.remove(o+".discarded.fastq")
    os.remove(o+".unassembled.forward.fastq")
    os.remove(o+".unassembled.reverse.fastq")
    

if __name__ == "__main__":
    test_files = [["raw/16S_101bp1.fq", "raw/16S_101bp2.fq", "raw/compare/101bp.assembled.fastq"], 
                  ["raw/16S_101bp1.fq.through.fastq", "raw/16S_101bp2.fq.through.fastq", "raw/compare/101bp.through.assembled.fastq"], 
                  ["raw/16S_150bp1.fq", "raw/16S_150bp2.fq", "raw/compare/150bp.assembled.fastq"], 
                  ["raw/16S_165bp1.fq", "raw/16S_165bp2.fq", "raw/compare/165bp.assembled.fastq"], 
                  ["raw/16S_180bp1.fq", "raw/16S_180bp2.fq", "raw/compare/180bp.assembled.fastq"], 
                  ["raw/16S_190bp1.fq", "raw/16S_190bp2.fq", "raw/compare/190bp.assembled.fastq"], 
                  ["raw/16S_250bp1.fq", "raw/16S_250bp2.fq", "raw/compare/250bp.assembled.fastq"], 
                  ["nooverlap/equall_1.fq", "nooverlap/equall_2.fq", "nooverlap/compare/equall.assembled.fastq"], 
                  ["nooverlap/f_long_1.fq", "nooverlap/f_long_2.fq", "nooverlap/compare/f_long.assembled.fastq"], 
                  ["nooverlap/r_long_1.fq", "nooverlap/r_long_2.fq", "nooverlap/compare/r_long.assembled.fastq"], 
                  ["normal/equall_1.fq", "normal/equall_2.fq", "normal/compare/equall.assembled.fastq"], 
                  ["normal/f_long_1.fq", "normal/f_long_2.fq", "normal/compare/f_long.assembled.fastq"], 
                  ["normal/r_long_1.fq", "normal/r_long_2.fq", "normal/compare/r_long.assembled.fastq"], 
                  ["runthrough/equall_1.fq", "runthrough/equall_2.fq", "runthrough/compare/equall.assembled.fastq"], 
                  ["runthrough/f_long_c2_1.fq", "runthrough/f_long_c2_2.fq", "runthrough/compare/f_long_c2.assembled.fastq"], 
                  ["runthrough/f_long_c3_1.fq", "runthrough/f_long_c3_2.fq", "runthrough/compare/f_long_c3.assembled.fastq"], 
                  ["runthrough/r_long_c2_1.fq", "runthrough/r_long_c2_2.fq", "runthrough/compare/r_long_c2.assembled.fastq"], 
                  ["runthrough/r_long_c3_1.fq", "runthrough/r_long_c3_2.fq", "runthrough/compare/r_long_c3.assembled.fastq"], 
                  ]
    i = 1
    for test in test_files:
        #gen(test[0], test[1], "temp", test[2])
        flag = unit_test(test[0], test[1], "temp", test[2])
        print("Test " + repr(i) + " " + flag)
        i = i + 1
