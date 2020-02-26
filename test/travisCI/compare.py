import re
import os
import sys


control_file_1 = "./test/control_1.bed"
treatment_file_1 = "./test/treatment_1.bed"

control_file_2 = "./test/control_2.bed"
treatment_file_2 = "./test/treatment_2.bed"

current_dir = "./test/travisCI/"

sicer_output_files_suffix = ['-W200-G600.scoreisland', '-W200-G600-FDR0.01-island.bed',
                '-W200-G600-FDR0.01-islandfiltered.bed', '-W200-G600-FDR0.01-islandfiltered-normalized.wig',
                '-W200-G600-islands-summary', '-W200-normalized.wig']
sicerdf_output_files_suffix = ['-W200-G600-summary', '-W200-G600-E1000-union.island', '-W200-G600-decreased-islands-summary-FDR0.01', '-W200-G600-increased-islands-summary-FDR0.01']


recog_output_files_suffix = ['-W200.cgisland', '-W200-FDR0.01-island.bed',
                '-W200-FDR0.01-islandfiltered.bed', '-W200-FDR0.01-islandfiltered-normalized.wig',
                '-W200-islands-summary', '-W200-normalized.wig']
recogdf_output_files_suffix = ['-W200-summary', '-W200-union.island', '-W200-decreased-islands-summary-FDR0.01', '-W200-increased-islands-summary-FDR0.01']



def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def check_13columns(file1_name,file2_name):
    file1 = open(file1_name,'r')
    file1.readline()
    file2 = open(file2_name,'r')
    file2.readline()
    equal = True
    for line1, line2 in zip(file1,file2):
        line1 = re.split('\t',line1)
        line1[4]= float(line1[4])
        line1[6]= float(line1[6])
        line1[7]= float(line1[7])
        line1[8]= float(line1[8])
        line1[9]= float(line1[9])
        line1[10]= float(line1[10])
        line1[11]= float(line1[11])
        line1[12]= float(line1[12].replace('\n',''))
        line2 = re.split('\t',line2)
        line2[4]= float(line2[4])
        line2[6]= float(line2[6])
        line2[7]= float(line2[7])
        line2[8]= float(line2[8])
        line2[9]= float(line2[9])
        line2[10]= float(line2[10])
        line2[11]= float(line2[11])
        line2[12]= float(line2[12].replace('\n',''))

        test4=False
        test6=False
        test7=False
        test8=False
        test9=False
        test10=False
        test11=False
        test12=False

        if(line1[4]==0 or line2[4]==0):
            test4=(line1[4]==line2[4])
        else:
            test4=isclose(line1[4],line2[4])

        if(line1[6]==0 or line2[6]==0):
            test6=(line1[6]==line2[6])
        else:
            test6=isclose(line1[6],line2[6])

        if(line1[7]==0 or line2[7]==0):
            test7=(line1[7]==line2[7])
        else:
            test7=isclose(line1[7],line2[7])

        if(line1[8]==0 or line2[8]==0):
            test8=(line1[8]==line2[8])
        else:
            test8=isclose(line1[8],line2[8])

        if(line1[9]==0 or line2[9]==0):
            test9=(line1[9]==line2[9])
        else:
            test9=isclose(line1[9],line2[9])

        if(line1[10]==0 or line2[10]==0):
            test10=(line1[10]==line2[10])
        else:
            test10=isclose(line1[10],line2[10])

        if(line1[11]==0 or line2[11]==0):
            test11=(line1[11]==line2[11])
        else:
            test11=isclose(line1[11],line2[11])

        if(line1[12]==0 or line2[12]==0):
            test12=(line1[12]==line2[12])
        else:
            test12=isclose(line1[12],line2[12])

        equal = ((line1[0]==line2[0]) and (line1[1]==line2[1]) and (line1[2]==line2[2]) and (line1[3]==line2[3]) and (line1[5]==line2[5])
                and test4 and test6 and test7 and test8 and test9 and test10 and test11 and test12)

        if(equal == False):
            #print(line1,line2)
            #print(test4, test6, test7, test8, test9, test10, test11, test12)
            break
    #print(equal)
    return equal


def check_unionisland(file1_name,file2_name):
    file1 = open(file1_name,'r')
    file2 = open(file2_name,'r')
    equal = True
    for line1, line2 in zip(file1,file2):
        line1 = re.split('\t',line1)
        line2 = re.split('\t',line2)

        equal = ((line1[0]==line2[0]) and (line1[1]==line2[1]) and (line1[2]==line2[2]))

        if (equal == False):
            print(line1,line2)
            break
    print(equal)
    return equal


def check_islandsummary (file1_name, file2_name):
    file1 = open(file1_name,'r')
    file2 = open(file2_name,'r')
    equal = True
    for line1, line2 in zip(file1,file2):
        line1 = re.split('\t',line1)
        line1[5]= float(line1[5])
        line1[6]= float(line1[6])
        line1[7]= float(line1[7].replace('\n',''))
        line2 = re.split('\t',line2)
        line2[5]= float(line2[5])
        line2[6]= float(line2[6])
        line2[7]= float(line2[7].replace('\n',''))


        equal = ((line1[0]==line2[0]) and (line1[1]==line2[1]) and (line1[2]==line2[2]) and (line1[3]==line2[3])
                and (line1[4]==line2[4]) and (isclose(line1[5], line2[5], abs_tol=1e-9))
                and (isclose(line1[6], line2[6], abs_tol=1e-5)) and (isclose(line1[7], line2[7], abs_tol=1e-7)))

    return equal


def check_WIG (file1_name, file2_name):
    file1 = open(file1_name,'r')
    file2 = open(file2_name,'r')
    equal = True
    for line1, line2 in zip(file1,file2):
        if(re.match("^track",line1)):
            equal = (line1 == line2)
        elif(re.match("^variableStep",line1)):
            equal = (line1 == line2)
        elif(line1!=""):
            line1 = re.split('\t',line1)
            line1[1]= float(line1[1].replace('\n',''))
            line2 = re.split('\t',line2)
            line2[1]= float(line2[1].replace('\n',''))

            equal == (line1[0]==line2[0]) and (isclose(line1[1], line2[1], abs_tol=0.01))

    return equal

def check_filteredbed(file1_name,file2_name):
    file1 = open(file1_name,'r')
    file2 = open(file2_name,'r')
    equal = True
    for line1, line2 in zip(file1,file2):
        line1 = re.split('\t',line1)
        line2 = re.split('\t',line2)

        equal = ((line1[0]==line2[0]) and (line1[1]==line2[1]) and (line1[2]==line2[2]) and (line1[5]==line2[5]))

    return equal

def check_islandbed(file1_name, file2_name):
    file1 = open(file1_name,'r')
    file2 = open(file2_name,'r')
    equal = True
    for line1, line2 in zip(file1,file2):
        line1 = re.split('\t',line1)
        line1[3]= line1[3].replace('\n','')
        line2 = re.split('\t',line2)
        line2[3]= line2[3].replace('\n','')

        equal = ((line1[0]==line2[0]) and (line1[1]==line2[1]) and (line1[2]==line2[2]) and (line1[3]==line2[3]))

    return equal


def check_scoreisland (file1_name, file2_name):
    file1 = open(file1_name,'r')
    file2 = open(file2_name,'r')
    equal = True
    for line1, line2 in zip(file1,file2):
        line1 = re.split('\t',line1)
        line1[3]= float(line1[3].replace('\n',''))
        line2 = re.split('\t',line2)
        line2[3]= float(line2[3].replace('\n',''))

        equal = (line1[0]==line2[0]) and (line1[1]==line2[1]) and (line1[2]==line2[2]) and (isclose(line1[3], line2[3], abs_tol=1e-6))

    return equal


def compare_sicer(f):
    if f == 1:
        tf = os.path.basename(treatment_file_1).replace(".bed",'')
    elif f == 2:
        tf = os.path.basename(treatment_file_2).replace(".bed",'')

    chk_score_island = check_scoreisland("./test/travisCI/test_output/" + tf+sicer_output_files_suffix[0], current_dir+'expected_output/'+tf+sicer_output_files_suffix[0])
    chk_island_bed = check_islandbed("./test/travisCI/test_output/" + tf+sicer_output_files_suffix[1], current_dir+'expected_output/'+tf+sicer_output_files_suffix[1])
    chk_filtered_bed = check_filteredbed("./test/travisCI/test_output/" + tf+sicer_output_files_suffix[2], current_dir+'expected_output/'+tf+sicer_output_files_suffix[2])
    chk_wig1 = check_WIG("./test/travisCI/test_output/" + tf+sicer_output_files_suffix[3], current_dir+'expected_output/'+tf+sicer_output_files_suffix[3])
    chk_island_summary = check_islandsummary("./test/travisCI/test_output/" + tf+sicer_output_files_suffix[4], current_dir+'expected_output/'+tf+sicer_output_files_suffix[4])
    chk_wig2 = check_WIG("./test/travisCI/test_output/" + tf+sicer_output_files_suffix[5], current_dir+'expected_output/'+tf+sicer_output_files_suffix[5])

    final_result = (chk_score_island and chk_island_bed and chk_filtered_bed and chk_wig1 and chk_island_summary and chk_wig2)
    return final_result

def compare_sicer_df():
    #Check for df execution
    f1_name = os.path.basename(treatment_file_1).replace(".bed",'')
    f2_name = os.path.basename(treatment_file_2).replace(".bed",'')

    summary_name = f1_name+'-and-'+f2_name+sicerdf_output_files_suffix[0]
    union_island_file_name = f1_name+'-vs-'+f2_name+sicerdf_output_files_suffix[1]
    decreased_fname = f1_name+sicerdf_output_files_suffix[2]
    increased_fname = f1_name+sicerdf_output_files_suffix[3]

    chk_unionisland = check_unionisland("./test/travisCI/test_output/" + union_island_file_name, current_dir+'expected_output/'+union_island_file_name)
    chk_summary = check_13columns("./test/travisCI/test_output/" + summary_name, current_dir+'expected_output/'+summary_name)
    chk_decreased = check_13columns("./test/travisCI/test_output/" + decreased_fname, current_dir+'expected_output/'+decreased_fname)
    chk_increased = check_13columns("./test/travisCI/test_output/" + increased_fname, current_dir+'expected_output/'+increased_fname)

    final_result = (chk_unionisland and chk_summary and chk_decreased and chk_increased)

    return final_result

def compare_recog(f):
    if f == 1:
        tf = os.path.basename(treatment_file_1).replace(".bed",'')
    elif f == 2:
        tf = os.path.basename(treatment_file_2).replace(".bed",'')

    chk_score_island = check_scoreisland("./test/travisCI/test_output/" + tf+recog_output_files_suffix[0], current_dir+'expected_output/'+tf+recog_output_files_suffix[0])
    chk_island_bed = check_islandbed("./test/travisCI/test_output/" + tf+recog_output_files_suffix[1], current_dir+'expected_output/'+tf+recog_output_files_suffix[1])
    chk_filtered_bed = check_filteredbed("./test/travisCI/test_output/" + tf+recog_output_files_suffix[2], current_dir+'expected_output/'+tf+recog_output_files_suffix[2])
    chk_wig1 = check_WIG("./test/travisCI/test_output/" + tf+recog_output_files_suffix[3], current_dir+'expected_output/'+tf+recog_output_files_suffix[3])
    chk_island_summary = check_islandsummary("./test/travisCI/test_output/" + tf+recog_output_files_suffix[4], current_dir+'expected_output/'+tf+recog_output_files_suffix[4])
    chk_wig2 = check_WIG("./test/travisCI/test_output/" + tf+recog_output_files_suffix[5], current_dir+'expected_output/'+tf+'-W200-cgnormalized.wig')

    final_result = (chk_score_island and chk_island_bed and chk_filtered_bed and chk_wig1 and chk_island_summary and chk_wig2)
    return final_result

def compare_recog_df():
    #Check for df execution
    f1_name = os.path.basename(treatment_file_1).replace(".bed",'')
    f2_name = os.path.basename(treatment_file_2).replace(".bed",'')

    summary_name = f1_name+'-and-'+f2_name+recogdf_output_files_suffix[0]
    union_island_file_name = f1_name+'-vs-'+f2_name+recogdf_output_files_suffix[1]
    decreased_fname = f1_name+recogdf_output_files_suffix[2]
    increased_fname = f1_name+recogdf_output_files_suffix[3]

    chk_unionisland = check_unionisland("./test/travisCI/test_output/" + union_island_file_name, current_dir+'expected_output/'+union_island_file_name)
    chk_summary = check_13columns("./test/travisCI/test_outpu/t" + summary_name, current_dir+'expected_output/'+summary_name)
    chk_decreased = check_13columns("./test/travisCI/test_outpu/t" + decreased_fname, current_dir+'expected_output/'+decreased_fname)
    chk_increased = check_13columns("./test/travisCI/test_output/" + increased_fname, current_dir+'expected_output/'+increased_fname)

    final_result = (chk_unionisland and chk_summary and chk_decreased and chk_increased)

    return final_result

if __name__ == "__main__":
    sicer_result_1 = compare_sicer(1)
    sicer_result_2 = compare_sicer(2)
    sicer_df_result = compare_sicer_df()
    sicer = (sicer_result_1 and sicer_result_2 and sicer_df_result)
    
    recog_result_1 = compare_recog(1)
    recog_result_2 = compare_recog(2)
    recog_df_result = compare_recog_df()
    recog = (recog_result_1 and recog_result_2 and recog_df_result)
    
    if sicer and recog:
        sys.exit(0)
    else:
        sys.exit(1)
