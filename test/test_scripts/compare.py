import re
import os
import sys
import subprocess

output_files_suffix = ['-W200-G600.scoreisland', '-W200-G600-FDR0.01-island.bed',
                       '-W200-G600-FDR0.01-islandfiltered.bed', '-W200-G600-FDR0.01-islandfiltered-normalized.wig',
                       '-W200-G600-islands-summary', '-W200-normalized.wig']

df_output_file_suffix = ['-W200-G600-summary', '-W200-G600-E1000-union.island',
                         '-W200-G600-decreased-islands-summary-FDR0.01', '-W200-G600-increased-islands-summary-FDR0.01']


recog_output_files_suffix = ['-W200.cgisland', '-W200-FDR0.01-island.bed',
                             '-W200-FDR0.01-islandfiltered.bed', '-W200-FDR0.01-islandfiltered-normalized.wig',
                             '-W200-islands-summary', '-W200-normalized.wig']
recogdf_output_files_suffix = ['-W200-summary', '-W200-union.island',
                               '-W200-decreased-islands-summary-FDR0.01', '-W200-increased-islands-summary-FDR0.01']


def isclose(a, b, rel_tol=1e-07, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def check_13columns(file1_name, file2_name):
    file1 = open(file1_name, 'r')
    file1.readline()
    file2 = open(file2_name, 'r')
    file2.readline()
    equal = True
    for line1, line2 in zip(file1, file2):
        line1 = re.split('\t', line1)
        line1[4] = float(line1[4])
        line1[6] = float(line1[6])
        line1[7] = float(line1[7])
        line1[8] = float(line1[8])
        line1[9] = float(line1[9])
        line1[10] = float(line1[10])
        line1[11] = float(line1[11])
        line1[12] = float(line1[12].replace('\n', ''))
        line2 = re.split('\t', line2)
        line2[4] = float(line2[4])
        line2[6] = float(line2[6])
        line2[7] = float(line2[7])
        line2[8] = float(line2[8])
        line2[9] = float(line2[9])
        line2[10] = float(line2[10])
        line2[11] = float(line2[11])
        line2[12] = float(line2[12].replace('\n', ''))

        test4 = False
        test6 = False
        test7 = False
        test8 = False
        test9 = False
        test10 = False
        test11 = False
        test12 = False

        if(line1[4] == 0 or line2[4] == 0):
            test4 = (line1[4] == line2[4])
        else:
            test4 = isclose(line1[4], line2[4])

        if(line1[6] == 0 or line2[6] == 0):
            test6 = (line1[6] == line2[6])
        else:
            test6 = isclose(line1[6], line2[6])

        if(line1[7] == 0 or line2[7] == 0):
            test7 = (line1[7] == line2[7])
        else:
            test7 = isclose(line1[7], line2[7])

        test8 = isclose(line1[8], line2[8])

        test9 = isclose(line1[9], line2[9])

        if(line1[10] == 0 or line2[10] == 0):
            test10 = (line1[10] == line2[10])
        else:
            test10 = isclose(line1[10], line2[10])

        if(line1[11] == 0 or line2[11] == 0):
            test11 = (line1[11] == line2[11])
        else:
            test11 = isclose(line1[11], line2[11])

        if(line1[12] == 0 or line2[12] == 0):
            test12 = (line1[12] == line2[12])
        else:
            test12 = isclose(line1[12], line2[12])

        equal = ((line1[0] == line2[0]) and (line1[1] == line2[1]) and (line1[2] == line2[2]) and (line1[3] == line2[3]) and (line1[5] == line2[5])
                 and test4 and test6 and test7 and test8 and test9 and test10 and test11 and test12)

        if(equal == False):
            print(line1, line2)
            print(test4, test6, test7, test8, test9, test10, test11, test12)
            break
    # print(equal)
    return equal


def check_unionisland(file1_name, file2_name):
    file1 = open(file1_name, 'r')
    file2 = open(file2_name, 'r')
    equal = True
    for line1, line2 in zip(file1, file2):
        line1 = re.split('\t', line1)
        line2 = re.split('\t', line2)

        equal = ((line1[0] == line2[0]) and (
            line1[1] == line2[1]) and (line1[2] == line2[2]))

        if (equal == False):
            print(line1, line2)
            break
    return equal


def check_islandsummary(file1_name, file2_name):
    file1 = open(file1_name, 'r')
    file2 = open(file2_name, 'r')
    equal = True
    for line1, line2 in zip(file1, file2):
        line1 = re.split('\t', line1)
        line1[5] = float(line1[5])
        line1[6] = float(line1[6])
        line1[7] = float(line1[7].replace('\n', ''))
        line2 = re.split('\t', line2)
        line2[5] = float(line2[5])
        line2[6] = float(line2[6])
        line2[7] = float(line2[7].replace('\n', ''))

        equal = ((line1[0] == line2[0]) and (line1[1] == line2[1]) and (line1[2] == line2[2]) and (line1[3] == line2[3])
                 and (line1[4] == line2[4]) and (isclose(line1[5], line2[5], abs_tol=1e-9))
                 and (isclose(line1[6], line2[6], abs_tol=1e-5)) and (isclose(line1[7], line2[7], abs_tol=1e-7)))

    return equal


def check_WIG(file1_name, file2_name):
    file1 = open(file1_name, 'r')
    file2 = open(file2_name, 'r')
    equal = True
    for line1, line2 in zip(file1, file2):
        if(re.match("^track", line1)):
            equal = (line1 == line2)
        elif(re.match("^variableStep", line1)):
            equal = (line1 == line2)
        elif(line1 != ""):
            line1 = re.split('\t', line1)
            line1[1] = float(line1[1].replace('\n', ''))
            line2 = re.split('\t', line2)
            line2[1] = float(line2[1].replace('\n', ''))

            equal == (line1[0] == line2[0]) and (
                isclose(line1[1], line2[1], abs_tol=0.01))

    return equal


def check_filteredbed(file1_name, file2_name):
    file1 = open(file1_name, 'r')
    file2 = open(file2_name, 'r')
    equal = True
    for line1, line2 in zip(file1, file2):
        line1 = re.split('\t', line1)
        line2 = re.split('\t', line2)

        equal = ((line1[0] == line2[0]) and (line1[1] == line2[1]) and (
            line1[2] == line2[2]) and (line1[5] == line2[5]))

    return equal


def check_islandbed(file1_name, file2_name):
    file1 = open(file1_name, 'r')
    file2 = open(file2_name, 'r')
    equal = True
    for line1, line2 in zip(file1, file2):
        line1 = re.split('\t', line1)
        line1[3] = line1[3].replace('\n', '')
        line2 = re.split('\t', line2)
        line2[3] = line2[3].replace('\n', '')

        equal = ((line1[0] == line2[0]) and (line1[1] == line2[1]) and (
            line1[2] == line2[2]) and (line1[3] == line2[3]))

    return equal


def check_scoreisland(file1_name, file2_name):
    file1 = open(file1_name, 'r')
    file2 = open(file2_name, 'r')
    equal = True
    for line1, line2 in zip(file1, file2):
        line1 = re.split('\t', line1)
        line1[3] = float(line1[3].replace('\n', ''))
        line2 = re.split('\t', line2)
        line2[3] = float(line2[3].replace('\n', ''))

        equal = (line1[0] == line2[0]) and (line1[1] == line2[1]) and (
            line1[2] == line2[2]) and (isclose(line1[3], line2[3], abs_tol=1e-6))

    return equal


def compare(new_path, old_path, file1):
    file1_name = file1.replace('bed', '')
    new_path += '/'
    old_path += '/'
    #subprocess.call(['module', 'load', 'bedtools'])
    try:
        print("calling bedtools")
        subprocess.call('bedtools sort -i %s > %s' % (old_path+file1_name +
                                                      output_files_suffix[0], (old_path+file1_name+output_files_suffix[0]+'.s')), shell=True)
        subprocess.call(['mv', (old_path+file1_name+output_files_suffix[0] +
                                '.s'), (old_path+file1_name+output_files_suffix[0])])
        print("done with bedtools")
    except:
        print("Error: bedtools not installed")
    #subprocess.call(['sort', '-k', '1', '-V', ('sorted_'+old_path+file1_name+output_files_suffix[0]), '>', old_path+file1_name+output_files_suffix[0]],shell=True)
    #subprocess.call(['rm', ('sorted_'+old_path+file1_name+output_files_suffix[0])])

    # Check file 1
    chk_score_island_1 = check_scoreisland(
        new_path+file1_name+output_files_suffix[0], old_path+file1_name+output_files_suffix[0])
    chk_island_bed_1 = check_islandbed(
        new_path+file1_name+output_files_suffix[1], old_path+file1_name+output_files_suffix[1])
    chk_filtered_bed_1 = check_filteredbed(
        new_path+file1_name+output_files_suffix[2], old_path+file1_name+output_files_suffix[2])
    chk_wig1_1 = check_WIG(
        new_path+file1_name+output_files_suffix[3], old_path+file1_name+output_files_suffix[3])
    chk_island_summary_1 = check_islandsummary(
        new_path+file1_name+output_files_suffix[4], old_path+file1_name+output_files_suffix[4])
    chk_wig2_1 = check_WIG(
        new_path+file1_name+output_files_suffix[5], old_path+file1_name+output_files_suffix[5])

    file1_result = chk_score_island_1 and chk_island_bed_1 and chk_filtered_bed_1 and chk_wig1_1 and chk_island_summary_1 and chk_wig2_1

    return file1_result


def df_compare(new_path, old_path, file1, file2):
    file1_name = file1.replace('.bed', '')
    file2_name = file2.replace('.bed', '')
    new_path += '/'
    old_path += '/'
    #subprocess.call(['module', 'load', 'bedtools'])
    ret_dir = os.getcwd()
    # os.chdir(old_path)
    try:
        print("calling bedtools")
        subprocess.call('bedtools sort -i %s > %s' % (old_path+file1_name +
                                                      output_files_suffix[0], (old_path+file1_name+output_files_suffix[0]+'.s')), shell=True)
        subprocess.call(['mv', (old_path+file1_name+output_files_suffix[0] +
                                '.s'), (old_path+file1_name+output_files_suffix[0])])
        print("done with bedtools")
    except:
        print("Error: bedtools not installed")
    #subprocess.run(['sort', '-k', '1', '-V', (old_path+file1_name+output_files_suffix[0]+'.s'), '>', old_path+file1_name+output_files_suffix[0]],shell=True)

    #subprocess.call(['rm', (old_path+file1_name+output_files_suffix[0]+'.s')])
    # os.chdir(ret_dir)

    # Check file 1
    print("=======Checking file 1=======")
    print("Comparing .scoreisland...")
    chk_score_island_1 = check_scoreisland(
        new_path+file1_name+output_files_suffix[0], old_path+file1_name+output_files_suffix[0])
    print("Comparing island .bed...")
    chk_island_bed_1 = check_islandbed(
        new_path+file1_name+output_files_suffix[1], old_path+file1_name+output_files_suffix[1])
    print("Comparing filtered island .bed...")
    chk_filtered_bed_1 = check_filteredbed(
        new_path+file1_name+output_files_suffix[2], old_path+file1_name+output_files_suffix[2])
    print("Comparing .wig...")
    chk_wig1_1 = check_WIG(
        new_path+file1_name+output_files_suffix[3], old_path+file1_name+output_files_suffix[3])
    print("Comparing island summary...")
    chk_island_summary_1 = check_islandsummary(
        new_path+file1_name+output_files_suffix[4], old_path+file1_name+output_files_suffix[4])
    print("Comparing normalized .wig...")
    chk_wig2_1 = check_WIG(
        new_path+file1_name+output_files_suffix[5], old_path+file1_name+output_files_suffix[5])

    file1_result = chk_score_island_1 and chk_island_bed_1 and chk_filtered_bed_1 and chk_wig1_1 and chk_island_summary_1 and chk_wig2_1
    print("File 1 result:", file1_result)

    #subprocess.call(['module', 'load', 'bedtools'])
    try:
        print("calling bedtools")
        subprocess.call('bedtools sort -i %s > %s' % (old_path+file2_name +
                                                      output_files_suffix[0], (old_path+file2_name+output_files_suffix[0]+'.s')), shell=True)
        subprocess.call(['mv', (old_path+file2_name+output_files_suffix[0] +
                                '.s'), (old_path+file2_name+output_files_suffix[0])])
        print("done with bedtools")
    except:
        print("Error: bedtools not installed")
    #subprocess.run(['sort', '-k', '1', '-V', (old_path+file2_name+output_files_suffix[0]+'.s'), '>', old_path+file2_name+output_files_suffix[0]],shell=True)
    #subprocess.call(['rm', (old_path+file1_name+output_files_suffix[0]+'.s')])

    # Check file 2
    print("=======Checking file 2=======")
    print("Comparing .scoreisland...")
    chk_score_island_2 = check_scoreisland(
        new_path+file2_name+output_files_suffix[0], old_path+file2_name+output_files_suffix[0])
    print("Comparing island .bed...")
    chk_island_bed_2 = check_islandbed(
        new_path+file2_name+output_files_suffix[1], old_path+file2_name+output_files_suffix[1])
    print("Comparing filtered island .bed...")
    chk_filtered_bed_2 = check_filteredbed(
        new_path+file2_name+output_files_suffix[2], old_path+file2_name+output_files_suffix[2])
    print("Comparing .wig...")
    chk_wig1_2 = check_WIG(
        new_path+file2_name+output_files_suffix[3], old_path+file2_name+output_files_suffix[3])
    print("Comparing island summary...")
    chk_island_summary_2 = check_islandsummary(
        new_path+file2_name+output_files_suffix[4], old_path+file2_name+output_files_suffix[4])
    print("Comparing normalized .wig...")
    chk_wig2_2 = check_WIG(
        new_path+file2_name+output_files_suffix[5], old_path+file2_name+output_files_suffix[5])

    file2_result = chk_score_island_2 and chk_island_bed_2 and chk_filtered_bed_2 and chk_wig1_2 and chk_island_summary_2 and chk_wig2_2
    print("File 2 result:", file2_result)

    # Check df result
    summary_name = file1_name+'-and-'+file2_name+df_output_file_suffix[0]
    union_island_file_name = file1_name+'-vs-' + \
        file2_name+df_output_file_suffix[1]
    old_union_island_name = file1_name+'-vs-' + \
        file2_name+'-W200-G600-E-union.island'
    decreased_fname = file1_name+df_output_file_suffix[2]
    increased_fname = file1_name+df_output_file_suffix[3]

    print("Comparing DF outputs...")
    chk_unionisland = check_unionisland(
        new_path+union_island_file_name, old_path+old_union_island_name)
    chk_summary = check_13columns(new_path+summary_name, old_path+summary_name)
    chk_decreased = check_13columns(
        new_path+decreased_fname, old_path+decreased_fname)
    chk_increased = check_13columns(
        new_path+increased_fname, old_path+increased_fname)

    df_result = chk_unionisland and chk_summary and chk_decreased and chk_increased

    return (file1_result and file2_result and df_result)


def recog_df_compare(new_path, old_path, file1, file2):
    file1_name = file1.replace('.bed', '')
    file2_name = file2.replace('.bed', '')
    new_path += '/'
    old_path += '/'

    try:
        print("calling bedtools")
        subprocess.call('bedtools sort -i %s > %s' % (old_path+file1_name +
                                                      recog_output_files_suffix[0], (old_path+file1_name+recog_output_files_suffix[0]+'.s')), shell=True)
        subprocess.call(['mv', (old_path+file1_name+recog_output_files_suffix[0] +
                                '.s'), (old_path+file1_name+recog_output_files_suffix[0])])
        print("done with bedtools")
    except:
        print("Error: bedtools not installed")

    # Check file 1
    print("=======Checking file 1=======")
    print("Comparing .scoreisland...")
    chk_score_island_1 = check_scoreisland(
        new_path+file1_name+recog_output_files_suffix[0], old_path+file1_name+recog_output_files_suffix[0])
    print("Comparing island .bed...")
    chk_island_bed_1 = check_islandbed(
        new_path+file1_name+recog_output_files_suffix[1], old_path+file1_name+recog_output_files_suffix[1])
    print("Comparing filtered island .bed...")
    chk_filtered_bed_1 = check_filteredbed(
        new_path+file1_name+recog_output_files_suffix[2], old_path+file1_name+recog_output_files_suffix[2])
    print("Comparing .wig...")
    chk_wig1_1 = check_WIG(
        new_path+file1_name+recog_output_files_suffix[3], old_path+file1_name+recog_output_files_suffix[3])
    print("Comparing island summary...")
    chk_island_summary_1 = check_islandsummary(
        new_path+file1_name+recog_output_files_suffix[4], old_path+file1_name+recog_output_files_suffix[4])
    print("Comparing normalized .wig...")
    chk_wig2_1 = check_WIG(
        new_path+file1_name+recog_output_files_suffix[5], old_path+file1_name+recog_output_files_suffix[5])

    file1_result = chk_score_island_1 and chk_island_bed_1 and chk_filtered_bed_1 and chk_wig1_1 and chk_island_summary_1 and chk_wig2_1
    print("File 1 result:", file1_result)

    #subprocess.call(['module', 'load', 'bedtools'])
    try:
        print("calling bedtools")
        subprocess.call('bedtools sort -i %s > %s' % (old_path+file2_name +
                                                      recog_output_files_suffix[0], (old_path+file2_name+recog_output_files_suffix[0]+'.s')), shell=True)
        subprocess.call(['mv', (old_path+file2_name+recog_output_files_suffix[0] +
                                '.s'), (old_path+file2_name+recog_output_files_suffix[0])])
        print("done with bedtools")
    except:
        print("Error: bedtools not installed")
    #subprocess.run(['sort', '-k', '1', '-V', (old_path+file2_name+recog_output_files_suffix[0]+'.s'), '>', old_path+file2_name+recog_output_files_suffix[0]],shell=True)
    #subprocess.call(['rm', (old_path+file1_name+recog_output_files_suffix[0]+'.s')])

    # Check file 2
    print("=======Checking file 2=======")
    print("Comparing .scoreisland...")
    chk_score_island_2 = check_scoreisland(
        new_path+file2_name+recog_output_files_suffix[0], old_path+file2_name+recog_output_files_suffix[0])
    print("Comparing island .bed...")
    chk_island_bed_2 = check_islandbed(
        new_path+file2_name+recog_output_files_suffix[1], old_path+file2_name+recog_output_files_suffix[1])
    print("Comparing filtered island .bed...")
    chk_filtered_bed_2 = check_filteredbed(
        new_path+file2_name+recog_output_files_suffix[2], old_path+file2_name+recog_output_files_suffix[2])
    print("Comparing .wig...")
    chk_wig1_2 = check_WIG(
        new_path+file2_name+recog_output_files_suffix[3], old_path+file2_name+recog_output_files_suffix[3])
    print("Comparing island summary...")
    chk_island_summary_2 = check_islandsummary(
        new_path+file2_name+recog_output_files_suffix[4], old_path+file2_name+recog_output_files_suffix[4])
    print("Comparing normalized .wig...")
    chk_wig2_2 = check_WIG(
        new_path+file2_name+recog_output_files_suffix[5], old_path+file2_name+recog_output_files_suffix[5])

    file2_result = chk_score_island_2 and chk_island_bed_2 and chk_filtered_bed_2 and chk_wig1_2 and chk_island_summary_2 and chk_wig2_2
    print("File 2 result:", file2_result)

    # Check df result
    summary_name = file1_name+'-and-'+file2_name+recogdf_output_files_suffix[0]
    union_island_file_name = file1_name+'-vs-' + \
        file2_name+recogdf_output_files_suffix[1]
    old_union_island_name = file1_name+'-vs-'+file2_name+'-W200-E-union.island'
    decreased_fname = file1_name+recogdf_output_files_suffix[2]
    increased_fname = file1_name+recogdf_output_files_suffix[3]

    print("Comparing DF outputs...")
    chk_unionisland = check_unionisland(
        new_path+union_island_file_name, old_path+old_union_island_name)
    chk_summary = check_13columns(new_path+summary_name, old_path+summary_name)
    chk_decreased = check_13columns(
        new_path+decreased_fname, old_path+decreased_fname)
    chk_increased = check_13columns(
        new_path+increased_fname, old_path+increased_fname)

    df_result = chk_unionisland and chk_summary and chk_decreased and chk_increased
    print(chk_unionisland)
    print(chk_summary)
    print(chk_decreased)
    print(chk_increased)

    return (file1_result and file2_result and df_result)
