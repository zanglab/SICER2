#!/usr/bin/env python
# Authors: Chongzhi Zang, Weiqun Peng
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


import bisect


def tag_position(line, fragment_size):
    shift = int(round(fragment_size / 2))
    strand = line[5]
    if strand == '+':
        return int(line[1]) + shift
    elif strand == '-':
        return int(line[2]) - 1 - shift


def find_readcount_on_islands(island_start_list, island_end_list, tag_position):
    """
    Make sure the islands are sorted.
    Islands are non-overlapping!
    Returns the index of the island on which the tag lands, or -1.
    """

    index = bisect.bisect_right(island_start_list, tag_position);
    if index - bisect.bisect_left(island_end_list, tag_position) == 1:
        return index - 1;
    else:
        return -1;
