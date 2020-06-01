import os
import subprocess
import re
import argparse
from itertools import zip_longest

FILES = ['45210_treat_rep1.bed','36626_treat_rep1.bed','36629_treat_rep1.bed','45206_treat_rep1.bed',
        '36615_treat_rep1.bed','36646_treat_rep1.bed','45213_treat_rep1.bed','36579_treat_rep1.bed',
        '36609_treat_rep1.bed','45216_treat_rep1.bed','36581_treat_rep1.bed','36598_treat_rep1.bed',
        '45220_treat_rep1.bed','36572_treat_rep1.bed','36582_treat_rep1.bed','45212_treat_rep1.bed',
        '45394_treat_rep1.bed','35394_treat_rep1.bed','35400_treat_rep1.bed','45393_treat_rep1.bed',
        '35393_treat_rep1.bed','45383_treat_rep1.bed','35389_treat_rep1.bed','35456_treat_rep1.bed',
        '45415_treat_rep1.bed','35417_treat_rep1.bed','35420_treat_rep1.bed','45384_treat_rep1.bed',
        '35408_treat_rep1.bed','35451_treat_rep1.bed','45404_treat_rep1.bed','35424_treat_rep1.bed',
        '35445_treat_rep1.bed']

gm12878_ctrl_group = set(['45210_treat_rep1.bed','36626_treat_rep1.bed','36629_treat_rep1.bed','45206_treat_rep1.bed',
        '36615_treat_rep1.bed','36646_treat_rep1.bed','45213_treat_rep1.bed','36579_treat_rep1.bed',
        '36609_treat_rep1.bed','45216_treat_rep1.bed','36581_treat_rep1.bed','36598_treat_rep1.bed',
        '45220_treat_rep1.bed','36572_treat_rep1.bed','36582_treat_rep1.bed','45212_treat_rep1.bed'])

sicer_suffix = [
    '-W200-normalized.wig',
    '-W200-G600.scoreisland',
    '-W200-G600-islands-summary',
    '-W200-G600-FDR0.01-island.bed',
    '-W200-G600-FDR0.01-islandfiltered.bed',
    '-W200-G600-FDR0.01-islandfiltered-normalized.wig',
    ]

sicer_df_suffix = [
    '-W200-G600-summary',
    '-W200-G600-E1000-union.island',
    '-W200-G600-decreased-islands-summary-FDR0.01',
    '-W200-G600-increased-islands-summary-FDR0.01'
]

recognicer_suffix = [
    '-W200-normalized.wig',
    '-W200.cgisland',
    '-W200-islands-summary',
    '-W200-FDR0.01-island.bed',
    '-W200-FDR0.01-islandfiltered.bed',
    '-W200-FDR0.01-islandfiltered-normalized.wig'
]

recognicer_df_suffix = [
    '-W200-summary',
    '-W200-union.island',
    '-W200-decreased-islands-summary-FDR0.01',
    '-W200-increased-islands-summary-FDR0.01'
]

def format_err_msg(line1, line2, f):
    file = os.path.basename(f)
    return 'Differing lines in ' + file + ':\n\"' + line1 + '\"   ||  \"' + line2 + '\"'

def isclose(a, b, rel_tol=1e-07, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def compare_files(f1, f2):
    base = os.path.basename(f1)
    file1 = open(f1, 'r')
    file2 = open(f2, 'r')

    for line1, line2 in zip_longest(file1, file2, fillvalue="EOF-EOF-EOF-EOF-EOF"):
        if line1 == "EOF-EOF-EOF-EOF-EOF" or line2 == "EOF-EOF-EOF-EOF-EOF":
            raise ValueError("Unequal line length in " + str(f1))

        line1 = line1.replace('\n', '')
        line2 = line2.replace('\n', '')
        line1_split = line1.split('\t')
        line2_split = line2.split('\t')
        if '' in line1_split:
            line1_split.remove('')
        if '' in line2_split:
            line2_split.remove('')
        if len(line1_split) != len(line2_split):
            raise ValueError(format_err_msg(line1, line2, f1))

        for i in range(len(line1_split)):
            if line1_split[i].isdigit():
                line1_split[i] = int(line1_split[i])
            else:
                try:
                    line1_split[i] = float(line1_split[i])
                except ValueError:
                    pass
            if line2_split[i].isdigit():
                line2_split[i] = int(line2_split[i])
            else:
                try:
                    line2_split[i] = float(line2_split[i])
                except ValueError:
                    pass

        line_equal = True
        for a, b in zip(line1_split, line2_split):
            if type(a) != type(b):
                if (type(a) is float and type(b) is int):
                    line_equal &= isclose(a, float(b), abs_tol=1e-6)
                elif (type(a) is int and type(b) is float):
                    line_equal &= isclose(float(a), b, abs_tol=1e-6)
                else:
                    str_a = str(a) + '(' + str(type(a)) + ')'
                    str_b = str(b) + '(' + str(type(b)) + ')' 
                    raise ValueError('Error: In \"' + base + '\", comparing different types: ' + str_a + ' and ' + str_b +'.')
            if type(a) is str or type(a) is int:
                line_equal &= a == b
            else:
                line_equal &= isclose(a, b, abs_tol=1e-6)

        if not line_equal:
            raise ValueError(format_err_msg(line1, line2, f1))

    return True 

def compare_wig_files(f1, f2):
    file1 = open(f1, 'r')
    file2 = open(f2, 'r')
    equal = True

    for line1, line2 in zip_longest(file1, file2, fillvalue="EOF-EOF-EOF-EOF-EOF"):
        line1 = line1.replace('\n', '')
        line2 = line2.replace('\n', '')

        if line1 != line2 and (line1 == "EOF-EOF-EOF-EOF-EOF" or line2 == "EOF-EOF-EOF-EOF-EOF"):
            raise ValueError("Error: Unequal # of length in \"" + f1 + "\"")

        line_equal = True
        if(re.match("^track", line1)):
            if (re.match("^track", line2)):
                line_equal &= line1 == line2
            else:
                line_equal = False
        elif(re.match("^variableStep", line1)):
            if(re.match("^variableStep", line2)):
                line_equal &= line1 == line2
            else:
                line_equal = False
        else:
            line1_split = line1.split('\t')
            line2_split = line2.split('\t')
            if len(line2_split) == 1:
                line2_split = line2.split(' ')
                while '' in line2_split:
                    line2_split.remove('')

            try:
                line1_split[1] = float(line1_split[1])
                line2_split[1] = float(line2_split[1])
            except:
                print(repr(line1_split))
                print(repr(line2_split))

            line_equal &= line1_split[0] == line2_split[0] and isclose(line1_split[1], line2_split[1], abs_tol=0.01)

        if not line_equal:
            raise ValueError(format_err_msg(line1, line2, f1))

    return True

def compare_sicer(file, output_dir, test_dir):
    base = os.path.splitext(os.path.basename(file))[0]
    result = True
    for suffix in sicer_suffix:
        target = output_dir + "/" + base + suffix
        truth = test_dir + "/" + base + suffix

        try:
            if suffix.endswith('.wig'):
                result &= compare_wig_files(target, truth)
            else:
                result &= compare_files(target, truth)
        except ValueError as err:
            result &= False
            print(err)
            break

    return result

def compare_sicer_df(output_dir, test_dir, file1, file2):
    base1 = os.path.splitext(os.path.basename(file1))[0]
    base2 = os.path.splitext(os.path.basename(file2))[0]
    result = compare_sicer(file1, output_dir, test_dir) and compare_sicer(file2, output_dir, test_dir)
    if not result:
        return result
    else:
        for suffix in sicer_df_suffix:
            if suffix.endswith('-summary'):
                target = output_dir + "/" + base1 + "-and-" + base2 + suffix
                truth = test_dir + "/" + base1 + "-and-" + base2 + suffix
            elif "union.island" in suffix:
                target = output_dir + "/" + base1 + "-vs-" + base2 + suffix
                truth = test_dir + "/" + base1 + "-vs-" + base2 + suffix
            else:
                target = output_dir + "/" + base1 + suffix
                truth = test_dir + "/" + base1 + suffix
            try:
                result &= compare_files(target, truth)
            except ValueError as err:
                result &= False
                print(err)
                break

    return result
    
def compare_recognicer(file, output_dir, test_dir):
    base = os.path.splitext(os.path.basename(file))[0]
    result = True
    for suffix in recognicer_suffix:
        target = output_dir + "/" + base + suffix
        truth = test_dir + "/" + base + suffix

        try:
            if ".wig" in suffix:
                result &= compare_wig_files(target, truth)
            else:
                result &= compare_files(target, truth)
        except ValueError as err:
            result &= False
            print(err)
            break

    return result

def compare_recognicer_df(output_dir, test_dir, file1, file2):
    base1 = os.path.splitext(os.path.basename(file1))[0]
    base2 = os.path.splitext(os.path.basename(file2))[0]
    result = compare_recognicer(file1, output_dir, test_dir) and compare_recognicer(file2, output_dir, test_dir)
    if not result:
        return result
    else:
        for suffix in recognicer_df_suffix:
            if suffix.endswith('-summary'):
                target = output_dir + "/" + base1 + "-and-" + base2 + suffix
                truth = test_dir + "/" + base1 + "-and-" + base2 + suffix
            elif "union.island" in suffix:
                target = output_dir + "/" + base1 + "-vs-" + base2 + suffix
                truth = test_dir + "/" + base1 + "-vs-" + base2 + suffix
            else:
                target = output_dir + "/" + base1 + suffix
                truth = test_dir + "/" + base1 + suffix
            try:
                result &= compare_files(target, truth)
            except ValueError as err:
                result &= False
                print(err)
                break

    return result

def run_sicer(treatment, control, output):
    run = subprocess.Popen(' '.join(['sicer', '-t', treatment, '-c', control, '-s', 'hg38', '-o', output, '--significant_reads', '>', '/dev/null']), stdin=subprocess.PIPE, shell=True)
    run.wait()
    if run.returncode != 0:
        raise Exception("SICER2 execution failure")
    return run.returncode

def run_recognicer(treatment, control, output):
    run = subprocess.Popen(' '.join(['recognicer', '-t', treatment, '-c', control, '-s', 'hg38', '-o', output, '--significant_reads', '>', '/dev/null']), stdin=subprocess.PIPE, shell=True)
    run.wait()
    if run.returncode != 0:
        raise Exception("SICER2 execution failure")
    return run.returncode

def run_sicer_df(t1, t2, c1, c2, output):
    run = subprocess.Popen(' '.join(['sicer_df', '-t', t1, t2, '-c', c1, c2, '-s', 'hg38', '-o', output, '--significant_reads', '>', '/dev/null']), stdin=subprocess.PIPE, shell=True)
    run.wait()
    if run.returncode != 0:
        raise Exception("SICER2 execution failure")
    return run.returncode

def run_recognicer_df(t1, t2, c1, c2, output):
    run = subprocess.Popen(' '.join(['recognicer_df', '-t', t1, t2, '-c', c1, c2, '-s', 'hg38', '-o', output, '--significant_reads', '>', '/dev/null']), stdin=subprocess.PIPE, shell=True)
    run.wait()
    if run.returncode != 0:
        raise Exception("SICER2 execution failure")
    return run.returncode

def df_test(test_type, data_dir, output_dir, test_dir, files, controls):

    faulty = False

    if files:
        t1 = data_dir + '/' + files[0]
        t2 = data_dir + '/' + files[1]
        c1 = data_dir + '/' + controls[0]
        c2 = data_dir + '/' + controls[1]
        if test_type == "sicer":
            print('Testing `sicer_df` with \"' + files[0] + '\" and \"' + files[1] + '\"...')
            run_sicer_df(t1, t2, c1, c2, output_dir)
            passed = compare_sicer_df(output_dir, test_dir, files[0], files[1])
        else:
            print('Testing `recognicer_df` with \"' + files[0] + '\" and \"' + files[1] + '\"...')
            run_recognicer_df(t1, t2, c1, c2, output_dir)
            passed = compare_recognicer_df(output_dir, test_dir, files[0], files[1])

        if passed:
            print("Test passed!")
        else:
            faulty = True

    else:

        for i in range(len(FILES),2):
            f1 = FILES[i]
            f2 = FILES[i+1]

            if f1 in gm12878_ctrl_group:
                c1 = 'GSM733742_GM12878_input.bed'
            else:
                c1 = 'GSM733780_K562_input.bed'

            if f2 in gm12878_ctrl_group:
                c2 = 'GSM733742_GM12878_input.bed'
            else:
                c2 = 'GSM733780_K562_input.bed'

            t1 = data_dir + '/' + f1
            t2 = data_dir + '/' + f1
            c1 = data_dir + '/' + c1
            c2 = data_dir + '/' + c2

            if test_type == "sicer":
                print('Testing `sicer_df` with \"' + f1 + '\" and \"' + f2 + '\"...')
                run_sicer_df(t1, t2, c1, c2, output_dir)
                passed = compare_sicer_df(output_dir, test_dir, f1, f2)
            else:
                print('Testing `recognicer_df` with \"' + f1 + '\" and \"' + f2 + '\"...')
                run_recognicer_df(t1, t2, c1, c2, output_dir)
                passed = compare_recognicer_df(output_dir, test_dir, f1, f2)

            if passed:
                print("Test passed!")
            else:
                faulty = True
                break

    if faulty:
        return False
    else:
        return True

def recognicer_test(data_dir, output_dir, test_dir, files=None, controls=None):

    faulty = False

    if file:
        print('Testing `recognicer` with \"' + file + '\"...')
        run_recognicer(data_dir + '/' + file, data_dir + '/' + control, output_dir)
        passed = compare_recognicer(file, output_dir, test_dir)

        if passed:
            print("Test passed!")
        else:
            faulty = True   

    else:
        for i in range(len(FILES)):
            file = data_dir + "/" + FILES[i]
            if file in gm12878_ctrl_group:
                control = 'GSM733742_GM12878_input.bed'
            else:
                control = 'GSM733780_K562_input.bed'
            control = data_dir + '/' + control

            print('Testing `recognicer` with \"' + file + '\"...')
            run_recognicer(file, control, output_dir)
            passed = compare_recognicer(file, output_dir, test_dir)

            if passed:
                print("Test passed!")
            else:
                faulty = True
                break
    if faulty:
        return False
    else:
        return True

def sicer_test(data_dir, output_dir, test_dir, file=None, control=None):

    faulty = False

    if file:
        print('Testing `sicer` with \"' + file + '\"...')
        run_sicer(data_dir + '/' + file, data_dir + '/' + control, output_dir)
        passed = compare_sicer(file, output_dir, test_dir)

        if passed:
            print("Test passed!")
        else:
            faulty = True        
    else:
        for i in range(len(FILES)):
            file = data_dir + "/" + FILES[i]
            if file in gm12878_ctrl_group:
                control = 'GSM733742_GM12878_input.bed'
            else:
                control = 'GSM733780_K562_input.bed'
            control = data_dir + '/' + control

            print('Testing `sicer` with \"' + file + '\"...')
            run_sicer(file, control, output_dir)
            passed = compare_sicer(file, output_dir, test_dir)

            if passed:
                print("Test passed!")
            else:
                faulty = True
                break
    if faulty:
        return False
    else:
        return True

def get_args():
    parser = argparse.ArgumentParser(description="Test SICER2")
    parser.add_argument("--sicer", action="store_true", help="If set, test SICER")
    parser.add_argument("--recognicer", action="store_true", help="If set, test RECOGNICER")
    parser.add_argument("--df", action="store_true", help="If set, test df mode")

    parser.add_argument("--data_dir", type=str, help="Path to location of data")
    parser.add_argument("--output_dir", type=str, help="Path to output results of SICER")
    parser.add_argument("--test_dir", type=str, help="Path to where ground-truth results are located")

    parser.add_argument('--test_file', type=str, nargs='+', help='Name of the file to test. Enter two files for df-testing. If not set, default to pre-set files')
    parser.add_argument('--control_file', type=str, nargs='+', help='Name of the corresponding control_file to test. If not set, default to pre-set files')

    return parser.parse_args()

if __name__ == "__main__":
    args = get_args()

    if not os.path.isdir(args.data_dir):
        raise ValueError("Directory \"" + args.data_dir + "\"doesn't exist.")

    if not os.path.isdir(args.test_dir):
        raise ValueError("Directory \"" + args.test_dir + "\"doesn't exist.")

    test_passed = True

    if args.sicer and not args.df:
        if args.test_file:
            test_passed &= sicer_test(args.data_dir, args.output_dir, args.test_dir, args.test_file[0], args.control_file[0])
        else:
            test_passed &= sicer_test(args.data_dir, args.output_dir, args.test_dir)

    if args.recognicer and not args.df:
        if args.test_file:
            test_passed &= recognicer_test(args.data_dir, args.output_dir, args.test_dir, args.test_file[0], args.control_file[0])
        else:
            test_passed &= recognicer_test(args.data_dir, args.output_dir, args.test_dir)

    if args.df:
        if args.sicer:
            test_passed &= df_test("sicer", args.data_dir, args.output_dir, args.test_dir, args.test_file, args.control_file)
        if args.recognicer:
            test_passed &= df_test("recognicer", args.data_dir, args.output_dir, args.test_dir, args.test_file, args.control_file)

    if test_passed:
        print("All tests passed!")
        exit(0)
    else:
        print("Test failed!")
        exit(1)
