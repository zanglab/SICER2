#!/usr/bin/env python
# Authors: Weiqun Peng, Chongzhi Zang, Dustin E Schones and Keji Zhao
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
# Version 1.1  11/9/2010


import re
from math import *
from optparse import OptionParser

"""
* identical tags are counted in
* tags on different strands are counted differently
    -> tag start + fragment size/2
    <- tag start - fragment size/2

* There are two sizes that matter.
1) a bin size for binning the tags: window_size
2) the fragment size after CHIP, which determine the physical
resolution of the CHIP experiment: fragment_size

================================================================

Read in the tag coords from the Solexa results, count the number of
tags in windows and make a bed graph plot showing the number of tags
in a window.

IMPORTANT: it is assumed that chromosomes are already separated

11/2/2010: the issue of the tag position potentially outside of chromsome is fixed.
Weiqun Peng

"""

plus = re.compile("\+");
minus = re.compile("\-");


def get_bed_coords(file, chrom_length, fragment_size):
    """
    *This takes into account the identical tags
    *Tags on different strands are positioned differently
        -> tag start (atoi(sline[1])) + fragment_size/2
        <- tag start (atoi(sline[2])) - 1 - fragment_size/2, the extra -1 is because that bed format has open-half, the sline[2] is not included.
    The stored positions are not the midpoint rathan than the start
    The interface is no longer the same as that for getBedCoords(file)
    input:
        file:  the file that has the raw tag data from one chromosome
        fragment_size: the fragment size after CHIP experiment.
    output:
        return: a sorted list of positions which might have redundent entries
    """
    infile = open(file);
    postive_tag_counts = 0.0;
    negative_tag_counts = 0.0;
    shift = int(round(fragment_size / 2));
    taglist = [];
    for line in infile:
        """ check to make sure not a header line """
        if not re.match("track", line):
            line = line.strip();
            sline = line.split();
            if int(sline[1]) < 0:
                print
                "Ilegitimate read with start less than zero is ignored";
                print
                line;
            elif int(sline[2]) >= chrom_length:
                print
                "Ilegitimate read with end beyond chromosome length ", chrom_length, " is ignored";
                print
                line;
            else:
                if plus.match(sline[5]):
                    position = int(sline[1]) + shift;
                    # If the position is beyond limit then don't shift.
                    if position >= chrom_length:
                        position = chrom_length - 1;
                    taglist.append(position);
                    postive_tag_counts += 1.0
                    elif minus.match(sline[5]):
                    position = int(sline[2]) - 1 - shift;
                    # in case the shift move the positions
                    # beyond zero, use zero
                    if position < 0: position = 0;  # UCSC genome coordinate is 0-based
                    taglist.append(position);
                    negative_tag_counts += 1.0;
    taglist.sort();
    """
    """
    total_tag_counts = postive_tag_counts + negative_tag_counts;
    print
    'total tag count in ' + file + ' is: ' + str(total_tag_counts) + ' = ' + str(postive_tag_counts) + '+' + str(
        negative_tag_counts);
    infile.close();
    return taglist;


def Generate_windows_and_count_tags(taglist, chrom, chrom_length, window_size, file):
    """
    taglist: sorted list of positions that includes every tag on a chromosome
    window_size: the artificial bin size for binning the tags
    bed_vals: a dictionary keyed by the start of tag_containing
        windows, with value being the tag count in the window.

    In this function, the bins are set up using an absolute coordinate
    system.  Namely [0, window_size-1),[window_size,
    2*window_size-1). If the last window goes beyond the limit of the chromosome,
    that window is ignored.

    The result writen into the file is guaranteed to be already sorted
    within a chromosome.
    """

    bed_vals = {};
    outfile = open(file, 'w');

    if (len(taglist) > 0):
        current_window_start = (taglist[0] / window_size) * window_size;
        tag_count_in_current_window = 1;
        for i in range(1, len(taglist)):
            start = (taglist[i] / window_size) * window_size;
            if start == current_window_start:
                tag_count_in_current_window += 1;
            elif start > current_window_start:
                # All the tags in the previous window have been counted
                current_window_end = current_window_start + window_size - 1;
                # if the window goes beyond the chromsome limit, it is discarded.
                if current_window_end < chrom_length:
                    bed_vals[current_window_start] = tag_count_in_current_window;
                    # write the window to file
                    outline = chrom + "\t" + str(current_window_start) + \
                              "\t" + str(current_window_end) + "\t" + \
                              str(tag_count_in_current_window) + "\n";
                    outfile.write(outline);
                current_window_start = start;
                tag_count_in_current_window = 1;
            else:
                print
                'Something is wrong!!!!!!!';

        current_window_end = current_window_start + window_size - 1;
        # if the window goes beyond the chromsome limit, it is discarded.
        if current_window_end < chrom_length:
            bed_vals[current_window_start] = tag_count_in_current_window;
            outline = chrom + "\t" + str(current_window_start) + \
                      "\t" + str(current_window_end) + "\t" + \
                      str(tag_count_in_current_window) + "\n";
            outfile.write(outline);
    outfile.close();
    return bed_vals;


def Total_number_of_windows(bed_vals):
    return len(bed_vals.keys());


def make_graph_file(tagfile, chrom, chrom_length, window_size, fragment_size, outfile):
    tag_list = get_bed_coords(tagfile, chrom_length, fragment_size);
    bed_vals = Generate_windows_and_count_tags(tag_list, chrom, chrom_length, window_size, outfile);


def main(argv):
    parser = OptionParser()
    parser.add_option("-f", "--tagfile", action="store", type="string",
                      dest="tagfile", help="file with tag coords in bed format",
                      metavar="<file>")
    parser.add_option("-c", "--chrom", action="store", type="string",
                      dest="chrom", help="chromosome name for graph",
                      metavar="<string>")
    parser.add_option("-l", "--chrom_length", action="store", type="int",
                      dest="chrom_length", help="length of chromosome",
                      metavar="<string>")
    parser.add_option("-w", "--window_size", action="store", type="int",
                      dest="window_size", help="window size to make summary",
                      metavar="<int>")
    parser.add_option("-i", "--fragment_size", action="store", type="int",
                      dest="fragment_size", metavar="<int>",
                      help="average size of a fragment after CHIP experiment")
    parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="output file name",
                      metavar="<file>")

    (opt, args) = parser.parse_args(argv)
    if len(argv) < 12:
        parser.print_help()
        sys.exit(1)

    tag_list = get_bed_coords(opt.tagfile, opt.chrom_length, opt.fragment_size);
    bed_vals = Generate_windows_and_count_tags(tag_list, opt.chrom, opt.chrom_length,
                                               opt.window_size, opt.outfile);


if __name__ == "__main__":
    main(sys.argv)
