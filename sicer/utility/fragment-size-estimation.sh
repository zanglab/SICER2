#!/bin/bash
# 
# Authors: Chongzhi Zang, Weiqun Peng
#
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).
#
# Version 1.1  6/9/2010

##############################################################
# ##### Please replace PATHTO with your own directory ###### #
##############################################################

SICER=PATHTO/SICER
PYTHONPATH=$SICER/lib
export PYTHONPATH

if [ $# -lt 1 ]; then
    echo ""
    echo 1>&2 Usage: $0 ["bed file"]
    echo ""
    exit 1
fi

# INPUT is the sample bed file
INPUT=$1
SAMPLE=${INPUT%.*}

# Directory where the INPUT bed file is. local as default
DIR=.

# RESOLUTION is the data resolution of the estimated fragment size (in bp). It is the window size of the summary graph file generated and the step size of the cross correlation calculation. Decreasing this number would significantly increase the computing time.
RESOLUTION=20
echo "########################################################################################"
echo "##  Estimating ChIP fragment size by cross correlation between watson and crick tags  ##"
echo "##                          An affiliated tool of SICER v1.1                          ##"
echo "########################################################################################"
echo "Data resolution: $RESOLUTION"

# Calculating tag correlation in chromosome 1
CHROM=chr1

# Output file name
DATAFILE=$SAMPLE-$CHROM-tag-correlation

# Grep chr1 tags
grep ${CHROM}[[:space:]] $DIR/$INPUT > temp1.bed
echo "Total tag count in $CHROM is"
wc -l temp1.bed

# Separate watson and crick reads
grep + temp1.bed > temp_plus.bed
grep - temp1.bed > temp_minus.bed

# Make summary graph files
echo "Positive"
python $SICER/lib/make_graph_file.py -f temp_plus.bed -c $CHROM -w $RESOLUTION -i 0 -o temp_plus.graph

echo "Negative"
python $SICER/lib/make_graph_file.py -f temp_minus.bed -c $CHROM -w $RESOLUTION -i 0 -o temp_minus.graph

# Calculate cross correlation between watson summary graph and crick summary graph
echo "Calculating correlation"
python $SICER/utility/calculate_cross_correlation_long_range.py -s hg18 -a temp_plus.graph -b temp_minus.graph -i $RESOLUTION -d $RESOLUTION -o $DATAFILE.txt

rm temp*

# Plotting
gnuplot << EOF
set terminal postscript eps color enhanced
set style fill solid 1.000000 border -1
set output "$DATAFILE.eps"
set size 0.6, 0.6
set title "$DATAFILE" font "Times-Roman, 20"
set xrange[:400]
plot "$DATAFILE.txt" using 1:2 notitle w l 1

EOF

ps2pdf $DATAFILE.eps
echo "Done! Please find the estimated ChIP fragment size by reading the peak location in the plot $DATAFILE.eps"
