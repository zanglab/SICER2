#!/usr/bin/env python
# Copyright (c) 2010 The George Washington University
# Author: Weiqun Peng
#
# This software is distributable under the terms of the GNU General
# Public License (GPL) v2, the text of which can be found at
# http://www.gnu.org/copyleft/gpl.html. Installing, importing or
# otherwise using this module constitutes acceptance of the terms of
# this License.
#
# Disclaimer
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).
#
# Version 1.1  6/9/2010


import re, os, sys, shutil
from math import *
from string import *

def get_total_tag_counts(tag_bed_file):
    """
    Get total tag counts given the current experimental run
    file should be a bed file.
    """
    counts =0;
    infile = open(tag_bed_file,'r');
    for line in infile:
            """ check to make sure not a header line """
            if not re.match("track", line):
                counts += 1;
    infile.close();
    return counts;

def get_total_tag_counts_bed_graph(summary_graph_file, bed_val={}, threshold = 0):
    """
    Get total tag counts given the current experimental run
    file should be a summary bed graph file.
    can be used for bed summary or islands file

    """
    count = 0;
    if summary_graph_file != "": #if summary graph exists
        infile = open(summary_graph_file,'r');
        for line in infile:
            """ check to make sure not a header line """
            if not re.match("track", line):
                line = line.strip();
                sline = line.split();
                assert ( len(sline) == 4 );
                value = float(sline[3]);
                if (value >= threshold):
                    count += value;
        infile.close();
    elif len(bed_val)>0:
        for chrom in bed_val.keys():
            for item in bed_val[chrom]:
                if (item.value >= threshold):
                    count +=item.value;
    else:
        print ("wrong input!!");
    return count;

def fileExists(f):
    try:
        file = open(f)
    except IOError:
        exists = 0
    else:
        exists = 1
    return exists

def is_bed_sorted(list):
    """
    Check if sorted in ascending order.
    input is a list of BED with chrom, start, end and value.
    output: sorted =1 or 0
    """
    sorted = 1;
    for index in range(0, len(list)-1):
        if list[index].start > list[index+1].start:
            sorted = 0;
    return sorted;

def is_list_sorted(list):
        """
        Check if sorted in ascending order.
        input is a list of values.
        output: sorted =1 or 0
        """
        sorted = 1;
        for index in range(0, len(list)-1):
                if list[index] > list[index+1]:
                        sorted = 0;
        return sorted;


def rescale_a_column(input_file, c, rescale_factor, output_file):
    """
    c is the 0-based column number
    Return a list of names

    """
    infile = open(input_file,'r')
    outfile = open(output_file, 'w')
    for line in infile:
        if not re.match("#", line):
            line = line.strip()
            sline = line.split()
            if (len(sline)>0):
                new_value = atof(sline[c]) * rescale_factor;
                sline[c] = str(new_value);
                outfile.write('\t'.join(sline)+'\n')
    infile.close()
    outfile.close()

def normalize_a_column(input_file, c, output_file):
    """
    c is the 0-based column number
    Return a list of names

    """
    line_number = 0;
    infile = open(input_file,'r')
    outfile = open(output_file, 'w')
    for line in infile:
        if not re.match("#", line):
            line = line.strip()
            sline = line.split()
            if (len(sline)>0):
                if (line_number == 0):
                    rescale_factor = atof(sline[c])
                new_value = atof(sline[c]) / rescale_factor;
                sline[c] = str(new_value);
                outfile.write('\t'.join(sline)+'\n')
                line_number += 1;
    infile.close()
    outfile.close()

def add_column (infile, c, outfile, const=-1):
    """
    c is the 0-based column number
    add a column to the original file
    default value would be the line number
    """
    file = open(infile,'r')
    ofile = open(outfile, 'w')
    counter = 0;
    for line in file:
        if not re.match("#", line):
            counter += 1
            line = line.strip()
            sline = line.split()
            if const == -1:
                sline.insert(c, "L" + str(counter));
            else:
                sline.insert(c, str(const));
            line =  '\t '.join(sline) + '\n';
            ofile.write(line);
    file.close()
    ofile.close()

def extract_two_columns(gene_file, c1, c2, outfile):
    """
    c is the 0-based column number
    """
    maxi = max (c1, c2);
    mini = min (c1, c2);
    file = open(gene_file,'r')
    ofile = open(outfile, 'w')
    for line in file:
        line = line.strip()
        sline = line.split()
        if len(sline)> maxi:
            outline = sline[c1] + '\t' + sline[c2] + '\n';
        elif len(sline)>mini:
            outline = sline[mini]+ '\n';
        ofile.write(outline);
    ofile.close();
    file.close();
