#!/usr/bin/python
# Authors: Weiqun Peng
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


from math import *


class Background_island_probscore_statistics:
    #  External genomeLength and gapSize are in units of bps
    #  Internal genomeLength and gapSize are in units of windows
    #  only look at enrichment!
    def __init__(self, total_tags, windowSize, gapSize, window_pvalue, genomeLength, bin_size):
        self.tag_density = total_tags * 1.0 / genomeLength;
        self.window_size = windowSize;  # In bps.
        self.gap_size = gapSize // windowSize;  # maximum number of windows allowed in a gap
        self.genome_length = int(ceil(float(genomeLength) / windowSize));
        self.average = self.tag_density * windowSize;  # lambda parameter for poisson distribution
        self.bin_size = bin_size;

        # Precalculate the poisson, cumulative poisson values up to max (500, 2*self.average) .
        self.max_index = max(500, int(2 * self.average));
        # print self.average, self.max_index;
        # self.fact=[];
        self.poisson_value = [];
        self.window_score = [];
        self.window_scaled_score = [];
        for index in range(self.max_index):
            # self.fact.append(self.factorial(index));
            prob = self.poisson(index, self.average);
            self.poisson_value.append(prob);
            if (index < self.average):  # only want to look at enrichment
                self.window_score.append(0);
                self.window_scaled_score.append(0);
            else:
                if prob > 0:
                    self.window_score.append(-log(prob));
                    # scaled_score =int(-log(prob)/self.bin_size);
                    scaled_score = int(round(-log(prob) / self.bin_size))
                    self.window_scaled_score.append(scaled_score);
                else:  # prob is too small and deemed 0 by the system
                    self.window_score.append(1000);
                    scaled_score = int(round(1000 / self.bin_size))
                    self.window_scaled_score.append(scaled_score);
            # print index, self.poisson_value[index], self.window_score[index];
        self.max_index = len(self.poisson_value);
        # print "max_index ", self.max_index;
        # gap_contribution needs min_tags_in_window
        # So the position of this line is critical.
        self.min_tags_in_window = 0;
        sf = 1;
        # print "self.poisson_value[0]=", self.poisson_value[0];
        while (sf > window_pvalue):
            # print self.min_tags_in_window, sf;
            sf -= self.poisson_value[self.min_tags_in_window]
            self.min_tags_in_window += 1;
        # An alternative approach that uses the scipy package,
        # poisson.sf (n, lambda) = \sum_{i= n+1}^{\infty} p(i, lambda)
        # self.min_tags_in_window = int(self.average);
        # while (scipy.stats.poisson.sf(self.min_tags_in_window-1) > window_pvalue):
        #   self.min_tags_in_window += 1;

        # print "Window read count threshold: ", self.min_tags_in_window;

        self.gap_contribution = self.gap_factor();
        self.boundary_contribution = self.boundary();

        self.cumulative = [];
        # new method, first fill the lowest score.
        prob = self.boundary_contribution * self.poisson_value[self.min_tags_in_window];
        score = -log(self.poisson_value[self.min_tags_in_window]);
        # scaled_score = int(score/self.bin_size);
        scaled_score = int(round(score / self.bin_size));
        self.island_expectation = [0] * (scaled_score + 1);
        self.island_expectation[scaled_score] = prob * self.genome_length;

        #       if len(self.island_expectation) < scaled_score:
        #           self.island_expectation += [0] *(scaled_score-len(self.island_expectation)+1);
        #           self.island_expectation[scaled_score] = prob*self.genome_length;
        # initial condition
        self.island_expectation[0] = self.boundary_contribution * self.genome_length / self.gap_contribution;

        self.root = self.find_asymptotics_exponent();
        # print "Exponent for Asymptotics: ", self.root;

    def factorial(self, m):
        value = 1.0
        if m != 0:
            while m != 1:
                value = value * m
                m = m - 1
        return value

    # Return the log of a factorial, using Srinivasa Ramanujan's approximation
    def factln(self, m):
        if m < 20:
            value = 1.0
            if m != 0:
                while m != 1:
                    value = value * m
                    m = m - 1
            return log(value)
        else:
            return m * log(m) - m + log(m * (1 + 4 * m * (1 + 2 * m))) / 6.0 + log(pi) / 2

    def poisson(self, i, average):
        if i < 20:
            return exp(-average) * average ** i / self.factorial(i)
        else:
            exponent = -average + i * log(average) - self.factln(i)
            return exp(exponent)

    """
        gap is in the unit of windows. In each window in the gap, the
        window could have 0, 1, min_tags_in_windows-1 tags.

        say gap = 1, min_tags_in_window= 2, gap_factor = 1 +
        poission(0,a) + poisson(1, a), where 1 represents no gap,
        poisson(0,a) represents a window with 0 tag,
        poisson(1,a) represents a window with 1 tag,

        The gap contribution from each window is not independent
    """

    def single_gap_factor(self):
        my_gap_factor = 0
        for i in range(self.min_tags_in_window):
            my_gap_factor += self.poisson_value[i];
        return my_gap_factor

    # gap contribution is bigger than 1
    def gap_factor(self):
        if self.gap_size == 0:
            return 1
        else:
            i = 1
            gap_contribution = 1  # contribution from no gap
            my_gap_factor = self.single_gap_factor()
            for i in range(1, self.gap_size + 1):
                gap_contribution += pow(my_gap_factor, i)
            return gap_contribution

    def boundary(self):
        """
        The condition for boundary is a continuous region of
        unqualified windows longer than gap
        """
        temp = self.single_gap_factor()
        temp = pow(temp, self.gap_size + 1)
        return temp * temp  # start & end

    # forward method that memorize the calculated results.
    def background_island_expectation(self, scaled_score):
        current_max_scaled_score = len(self.island_expectation) - 1
        if scaled_score > current_max_scaled_score:
            # index is the scaled_score
            for index in range(current_max_scaled_score + 1, scaled_score + 1):
                temp = 0.0
                # i is the number of tags in the added window
                i = self.min_tags_in_window
                while (int(round(index - self.window_score[i] / self.bin_size)) >= 0):
                    # while ( (index - self.window_scaled_score[i])>=0):
                    temp += self.poisson_value[i] * self.island_expectation[
                        int(round(index - self.window_score[i] / self.bin_size))]
                    # temp += self.poisson_value[i]* self.island_expectation[index - self.window_scaled_score[i]];
                    i += 1
                temp *= self.gap_contribution
                self.island_expectation.append(temp)
                # print index, temp, self.island_expectation[index];
        return self.island_expectation[scaled_score]

    def generate_cumulative_dist(self, outfile=""):
        """
        Generate cumulative distribution: a list of tuples (bins, hist).
        """
        self.cumulative = [0] * len(self.island_expectation)
        partial_sum = 0.0
        for index in range(1, len(self.island_expectation) + 1):
            complimentary = len(self.island_expectation) - index
            partial_sum += self.island_expectation[complimentary]  # The end is outside of the index
            self.cumulative[complimentary] = partial_sum

        if outfile != "":
            fixpoint = int(len(self.island_expectation) / 2)
            outf = open(outfile, "w")
            outline = "# Score" + "\t" + "Expect # islands" + "\t" + "Cumulative # Islands" + "\t" + "Asymptotics" + "\n"
            outf.write(outline)
            for index in range(len(self.island_expectation)):
                outline = str(index * self.bin_size) + "\t" + str(self.island_expectation[index]) + "\t" + str(
                    self.cumulative[index]) + "\n"
                # outline = str(index * self.bin_size) + "\t" + str(self.island_expectation[index])+ "\t" +str(self.cumulative[index]) + "\t" + str(self.cumulative[fixpoint] * exp(-self.root*(self.cumulative[index]-self.cumulative[fixpoint]))) + "\n";
                outf.write(outline)
            outf.close()

    def find_island_threshold(self, e_value_threshold):
        """
        average is the average number of tags in a window:
        opt.tag_density * opt.window_size

        This one allows single-window islands.
        Returns the island threshold
        """
        threshold = .0000001 * e_value_threshold
        current_scaled_score = len(self.island_expectation) - 1
        current_expectation = self.island_expectation[-1]
        assert (current_expectation == self.island_expectation[current_scaled_score]);
        interval = int(1 / self.bin_size)
        if len(self.island_expectation) > interval:
            partial_cumu = sum(self.island_expectation[-interval: -1])
        else:
            partial_cumu = sum(self.island_expectation)
        while (partial_cumu > threshold or partial_cumu < 1e-100):
            current_scaled_score += interval
            current_expectation = self.background_island_expectation(current_scaled_score)
            if len(self.island_expectation) > interval:
                partial_cumu = sum(self.island_expectation[-interval: -1])
            else:
                partial_cumu = sum(self.island_expectation)
            # for index in  xrange(len(self.island_expectation)):
            # print  index*self.bin_size, self.island_expectation[index];

        self.generate_cumulative_dist()
        for index in range(len(self.cumulative)):
            if self.cumulative[index] <= e_value_threshold:
                score_threshold = index * self.bin_size
                break
        return score_threshold

    def func(self, x):
        sum = 0.0
        for index in range(self.min_tags_in_window, self.max_index):
            sum += self.gap_contribution * pow(self.poisson_value[index], 1 - x)
        return sum - 1

    # from Mathematical Utility routines, based on numerical recipe
    #   Copyright (C) 1999, Wesley Phoa
    def bracket_root(self, f, interval, max_iterations=50):
        """\
        Given a univariate function f and a tuple interval=(x1,x2),
        return a new tuple (bracket, fnvals) where bracket=(x1,x2)
        brackets a root of f and fnvals=(f(x1),f(x2)).
        """
        GOLDEN = (1 + 5 ** .5) / 2
        (x1, x2) = interval
        if x1 == x2:
            raise Exception("initial interval has zero width")
        elif x2 < x1:
            x1, x2 = x2, x1
        f1, f2 = f(x1), f(x2)
        for j in range(max_iterations):
            while f1 * f2 >= 0:  # not currently bracketed
                if abs(f1) < abs(f2):
                    x1 = x1 + GOLDEN * (x1 - x2)
                else:
                    x2 = x2 + GOLDEN * (x2 - x1)
                f1, f2 = f(x1), f(x2)
            return (x1, x2), (f1, f2)
        raise Exception("too many iterations")

    # based on Numerical Recipes, p. 354
    def bisect_root(self, func, interval, xacc):
        JMAX = 50;
        (x1, x2) = interval
        f = func(x1)
        fmid = func(x2)
        if (f * fmid >= 0.0):
            print("Root must be bracketed for bisection")
        if (f < 0.0):
            dx = x2 - x1
            rtb = x1
        else:
            dx = x1 - x2
            rtb = x2
        for j in range(JMAX):
            dx *= 0.5
            xmid = rtb + dx
            fmid = func(xmid)
            if (fmid <= 0.0): rtb = xmid;
            if (fabs(dx) < xacc or fmid == 0.0): return rtb
        print("Too many bisections")
        return 0.0

    def find_asymptotics_exponent(self, xacc=.00001):
        num = 100
        # for index in xrange(num):
        #   x = index/float(num);
        #   print x, self.func(x);
        input_bracket = (0.1, 1)
        (xresult, yresult) = self.bracket_root(self.func, input_bracket)
        root = self.bisect_root(self.func, xresult, xacc)
        # print "# The exponent is: ", root;
        return root
