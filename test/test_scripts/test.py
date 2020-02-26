import os
import subprocess 
import time
import compare

core_count = 12

test_dir =  '/nv/vol190/zanglab/jy2ma/sicer2'
data_path = test_dir+'/data'
new_sicer_result_path = test_dir+'/sicer_result/new'
old_sicer_result_path = test_dir+'/sicer_result/old'
new_recog_result_path = test_dir+'/recog_result/new'
old_recog_result_path = test_dir+'/recog_result/old'
old_sicer_path = test_dir+'/SICER_1'

FILES = ["36629_treat_rep1.bed"]


"""
FILES = ['45210_treat_rep1.bed','36626_treat_rep1.bed','36629_treat_rep1.bed','45206_treat_rep1.bed',
        '36615_treat_rep1.bed','36646_treat_rep1.bed','45213_treat_rep1.bed','36579_treat_rep1.bed',
        '36609_treat_rep1.bed','45216_treat_rep1.bed','36581_treat_rep1.bed','36598_treat_rep1.bed',
        '45220_treat_rep1.bed','36572_treat_rep1.bed','36582_treat_rep1.bed','45212_treat_rep1.bed',
        '45394_treat_rep1.bed','35394_treat_rep1.bed','35400_treat_rep1.bed','45393_treat_rep1.bed',
        '35393_treat_rep1.bed','45383_treat_rep1.bed','35389_treat_rep1.bed','35456_treat_rep1.bed',
        '45415_treat_rep1.bed','35417_treat_rep1.bed','35420_treat_rep1.bed','45384_treat_rep1.bed',
        '35408_treat_rep1.bed','35451_treat_rep1.bed','45404_treat_rep1.bed','35424_treat_rep1.bed',
        '35445_treat_rep1.bed']
"""

"""
gm12878_ctrl_group = set(['45210_treat_rep1.bed','36626_treat_rep1.bed','36629_treat_rep1.bed','45206_treat_rep1.bed',
        '36615_treat_rep1.bed','36646_treat_rep1.bed','45213_treat_rep1.bed','36579_treat_rep1.bed',
        '36609_treat_rep1.bed','45216_treat_rep1.bed','36581_treat_rep1.bed','36598_treat_rep1.bed',
        '45220_treat_rep1.bed','36572_treat_rep1.bed','36582_treat_rep1.bed','45212_treat_rep1.bed'])
"""


def compare_files(f1, f2):
    file1 = open(file1_name, 'r')
    file2 = open(file2_name, 'r')

    for line1, line2 in zip(file1, file2):
        line1_split = line1.replace('\n', '').split('\t')
        line2_split = line2.replace('\n', '').split('\t')
        if len(line1_split) != len(line_split)
            raise ValueError(f'Different # of Columns: \"{line1}\" and \"{line2}\"')

        for i in range(len(line1_split)):
            if line1_split[i].isDigit():
                line1_split[i] = int(line1_split[i])
            else:
                try:
                    line1_split[i] = float(line1_split[i])
                except ValueError:
                    pass
            if line2_split[i].isDigit():
                line2_split[i] = int(line2_split[i])
            else:
                try:
                    line2_split[i] = float(line2_split[i])
                except ValueError:
                    pass

        line_equal = True
        for a, b in zip(line1_split, line2_split):
            if type(a) != type(b):
                raise ValueError(f'Comparing different types: {a}({type(a)}) and {b}({type(b)})')
            if type(a) is str or type(a) is int:
                line_equal &= a == b
            else:
                line_equal &= isclose(a, b, abs_tol=1e-6)

        if not line_equal:
            raise ValueError(f'Differing lines: \"{line1}\" and \"{line2}\"')
    return True 

def compare_wig_files(f1, f2):
    file1 = open(file1_name, 'r')
    file2 = open(file2_name, 'r')
    equal = True
    for line1, line2 in zip(file1, file2):
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
        elif(line1 != ""):
            line1_split = line1.replace('\n', '').split('\t')
            line1_split[1] = float(line1_split[1])
            line2 = line2.replace('\n', '').split('\t')
            line2_split[1] = float(line2_split[1])

            line_equal &= line1_split[0] == line2_split[0] and isclose(line1_split[1], line2_split[1], abs_tol=0.01)

        if not line_equal:
            raise ValueError(f'Differing lines: \"{line1}\" and \"{line2}\"')

    return True

def compare_sicer():
    pass

def compare_sicer_df():
    pass

def compare_recognicer():

def compare_recognicer_df():


def run_sicer(treatment, control, output):
    subprocess.call(['sicer', '-t', treatment, '-c', control, '-s', 'hg38', '-o', output, '--significant_reads', '>', '/dev/null'])

def run_recognicer(treatment, control, output):
    subprocess.call(['recognicer', '-t', treatment, '-c', control, '-s', 'hg38', '-o', output, '--significant_reads', '>', '/dev/null'])

def run_sicer_df(t1, t2, c1, c2, output):
    subprocess.call(['sicer_df', '-t', t1, t2, '-c', c1, c2, '-s', 'hg38', '-o', output, '--significant_reads', '>', '/dev/null'])
  

def run_recognicer_df(t1, t2, c1, c2, output):
    subprocess.call(['recognicer_df', '-t', t1, t2, '-c', c1, c2, '-s', 'hg38', '-o', output, '--significant_reads', '>', '/dev/null'])

def df_test(test_type, data_dir, output_dir, test_dir, file=None):

    faulty = False

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
            print(f"Testing `sicer_df` with \"{f1}\" and \"{f2}\"...")
            run_sicer_df(t1, t2, c1, c2, output_dir)
            passed = compare.check_sicer_df(output_dir, test_dir, f1, f2)
        else:
            print(f"Testing `recognicer_df` of \"{f1}\" and \"{f2}\"...")
            run_recognicer_df(t1, t2, c1, c2, output_dir)
            passed = compare.check_recognicer_df(output_dir, test_dir, f1, f2)

        if passed:
            print("Test passed!")

        else passed:
            faulty = True
            break

    if faulty:
        return False
    else:
        return True

def recognicer_test(data_dir, output_dir, test_dir, file=None):

    faulty = False

    for i in range(len(FILES)):
        file = FILES[i]
        if file in gm12878_ctrl_group:
            control = 'GSM733742_GM12878_input.bed'
        else:
            control = 'GSM733780_K562_input.bed'

        print(f"Testing `recognicer` with \"{file}\"...")
        run_recognicer(file, control, output_dir)
        passed = compare.check_recognicer(output_dir, test_dir, file)

        if passed:
            print("Test passed!")
        else:
            faulty = True
            break
    if faulty:
        return False
    else:
        return True

def sicer_test(data_dir, output_dir, test_dir, file=None):

    faulty = False

    for i in range(len(FILES)):
        file = FILES[i]
        if file in gm12878_ctrl_group:
            control = 'GSM733742_GM12878_input.bed'
        else:
            control = 'GSM733780_K562_input.bed'

        print(f"Testing `sicer` with \"{file}\"...")
        run_sicer(file, control, output_dir)
        passed = compare.check_sicer(output_dir, test_dir, file)

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

    parser.add_argument('--test_file', type=str, nargs='+', help='Name of the files to test. If not set, default to pre-set files')
    parser.add_argument('--control_file', type=str, nargs='+', help='Name of the corresponding control_file to test. If not set, default to pre-set files')

if __name__ == "__main__":
    args = parser.parse_args()

    if not os.path.isdir(args.data_dir):
        raise ValueError(f"Directory `{args.data_dir}` doesn't exist.")

    if not os.path.isdir(args.output_dir):
        raise ValueError(f"Directory `{args.output_dir}` doesn't exist.")

    if not os.path.isdir(args.test_dir):
        raise ValueError(f"Directory `{args.test_dir}` doesn't exist.")

    test_passed = 0

    if test_sicer and not test_df:
        test_passed &= sicer_test(args.data_dir, args.output_dir, args.test_dir)

    if test_recognicer and not test_df:
        test_passed &= recognicer_test(args.data_dir, args.output_dir, args.test_dir)

    if test_df:
        if test_sicer:
            test_passed &= df_test("sicer", args.data_dir, args.output_dir, args.test_dir)
        else:
            test_passed &= df_test("recognicer", args.data_dir, args.output_dir, args.test_dir)

    if test_passed:
        print("All tests passed!")
    else:
        print("Test failed!")