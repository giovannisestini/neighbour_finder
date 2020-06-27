#!/usr/bin/env python

#This script allows to identify the combinations and the frequencies of SalI boxes.
#The input file has to be in the .txt file format. Each line has to be occupied by one BAC. Each repeat has to be separated by /tab. Look to the following example of a file containing 2 BACs:
#
#3(EEV)111	3(EEV)001	3(EEV)111	3(EEV)001
#2(EV)011	3(ELV)011	2(EZ)011	2(EE)011
#
#In the example, the file contains the barcode of 2 BACs, each of them consist of 4 repeats. Each repeat is separated from the previous one by tab.

import argparse, random, os, json, sys
from collections import Counter

def sal_extractor(in_file, out_file, random_shuffle, typ_analysis) : # it extracts the portion of the barcoding corresponding to SalI boxes or SNPs/promoters
    with open(in_file, 'r') as f:
        if typ_analysis == 0 and random_shuffle == 0: #salI box analysis without randomization
            for element in f :
                x = element.split('\t')
                for sal in x : # it extracts the portion of the code corresponding to SalI boxes, which is found between round brackets
                    y = sal.split('(') 
                    z = y[1].split(')')
                    out_file.write(z[0] + '\t') # Just the portion of the barcode corresponding to SalI boxes is written in the temporary file
                out_file.write('\n')
        elif typ_analysis == 0 and random_shuffle == 1 : #SalI box analysis with randomization
            random_list = []
            for element in f :
                x = element.split('\t')
                for sal in x : #it extracts the portion of the code corresponding to SalI boxes, which is found between round brackets
                    y = sal.split('(')
                    z = y[1].split(')')
                    random_list.append(z[0]) # SalI boxes of the BAC under analysis are randomized
                random.shuffle(random_list) # All SalI boxes from all the BACs are randomized
            value = 0 
            sal_value = 0 
            for x in range (0, len(random_list)) : # Just the portion of the barcode corresponding to SalI boxes is written in the temporary file
                out_file.write(random_list[sal_value] + '\t')
                value += 1
                sal_value += 1 
                if value == 8 : # Up to 8 SalI boxes are written on the same line before newline
                    out_file.write('\n')
                    value = 0 
            out_file.write('\n')
        elif typ_analysis == 1 and random_shuffle == 0 : #SNP analysis without randomization
            for element in f :
                x = element.split('\t')
                for sal in x : # it extracts the portion of the code corresponding to SNPs/promoters, which is found outside round brackets
                    y = sal.split('(')
                    z = y[1].split(')')
                    out_file.write(y[0] + z[1] + '\t')
                out_file.write('\n')
        elif typ_analysis == 1 and random_shuffle == 1 : #SNP analysis with randomization
            random_list = []
            for element in f :
                x = element.split('\t')
                for sal in x : # it extracts the portion of the code corresponding to SNPs/promoters, which is found outside round brackets
                    y = sal.split('(')
                    z = y[1].split(')')
                    w = z[1].split('\n')
                    random_list.append(y[0] + w[0]) # The SNPs of the BAC under analysis are randomized
                random.shuffle(random_list) # All SNPs from all the BACs are randomized
            value = 0
            sal_value = 0
            for x in range (0, len(random_list)) : # Just the portion of the barcode corresponding to SNPs is written in the temporary file
                out_file.write(random_list[sal_value] + '\t')
                value += 1
                sal_value += 1
                if value == 8 : # Up to 8 SNPs are written on the same line before newline
                    out_file.write('\n')
                    value = 0
            out_file.write('\n')
        elif typ_analysis == 2 and random_shuffle == 0 : # AvaI + CAT
            for element in f :
                x = element.split('\t')
                for sal in x :
                    y = sal.split('(')
                    z = y[1].split(')')
                    w = y[0] + z[1]
                    try :
                        out_file.write(w[1] + w[2] + '\t')
                    except IndexError :
                        pass
                out_file.write('\n')
        elif typ_analysis == 2 and random_shuffle == 1 :
            random_list = []
            for element in f :
                x = element.split('\t')
                for sal in x :
                    y = sal.split('(')
                    z = y[1].split(')')
                    w = y[0] + z[1]
                    try :
                        k = w[1] + w[2]
                    except IndexError :
                        pass
                    random_list.append(k)
                random.shuffle(random_list)
            value = 0
            sal_value = 0
            for x in range (0, len(random_list)) : # Just the portion of the barcode corresponding to SNPs is written in the temporary file
                out_file.write(random_list[sal_value] + '\t')
                value += 1
                sal_value += 1
                if value == 8 : # Up to 8 features are written on the same line before newline
                    out_file.write('\n')
                    value = 0
            out_file.write('\n')
        elif typ_analysis == 3 and random_shuffle == 0 : # whole repeat analysis
            for element in f :
                    out_file.write(element)
        elif typ_analysis == 3 and random_shuffle == 1 : # whole repeat anlysis random
            random_list = []
            for element in f :
                x = element.split('\n')
                y = x[0].split('\t')
                for sal in y :
                    random_list.append(sal)
                random.shuffle(random_list)
            value = 0
            sal_value = 0   
            for x in range (0, len(random_list)) :
                out_file.write(random_list[sal_value] + '\t')
                value += 1
                sal_value += 1
                if value == 8 :
                    out_file.write('\n')
                    value = 0
            out_file.write('\n')
        elif typ_analysis == 4 and random_shuffle == 0 : # ava/cat/3'ets analysis
            for element in f :
                x = element.split('\t')
                for sal in x :
                    y = sal.split('(')
                    z = y[1].split(')')
                    try :
                        w = z[1][0] + z[1][1] + z[1][2]
                        out_file.write(w + '\t')
                    except IndexError :
                        pass
                out_file.write('\n')
        elif typ_analysis == 4 and random_shuffle == 1 : # ava/cat/3ets analysis random
            random_list = []
            for element in f :
                x = element.split('\t')
                for sal in x :
                    y = sal.split('(')
                    z = y[1].split(')')
                    try :
                        w = z[1][0] + z[1][1] + z[1][2]
                        random_list.append(w)
                    except IndexError :
                        pass
                random.shuffle(random_list)
            value = 0
            sal_value = 0
            for x in range (0, len(random_list)) : # Just the portion of the barcode corresponding to SNPs is written in the temporary file
                out_file.write(random_list[sal_value] + '\t')
                value += 1
                sal_value += 1
                if value == 8 : # Up to 8 features are written on the same line before newline
                    out_file.write('\n')
                    value = 0
            out_file.write('\n')    
        elif random_shuffle == 0 : # ava + 3ets anlysis (5) cat +3ets analysis (6)
            if typ_analysis == 5 or typ_analysis == 6 : 
                for element in f :
                    x = element.split('\t')
                    for sal in x :
                        y = sal.split('(')
                        z = y[1].split(')')
                        w = y[0] + z[1]
                        if typ_analysis == 5 :
                            try :
                                out_file.write(w[1] + w[3] + '\t')
                            except IndexError :
                                pass
                        elif typ_analysis == 6 :
                            try :
                                out_file.write(w[2] + w[3] + '\t')
                            except IndexError :
                                pass
                        else :
                            print('ERROR')
                    out_file.write('\n')
            else :
                print('ERROR --analysis_type not correct')
        elif random_shuffle == 1 : # ava + 3ets anlysis (5), cat +3ets analysis (6) random
            if typ_analysis == 5 or typ_analysis == 6 :
                print('ok random')
                random_list = []
                for element in f :
                    x = element.split('\t')
                    for sal in x :
                        y = sal.split('(')
                        z = y[1].split(')')
                        w = y[0] + z[1]
                        if typ_analysis == 5 :
                            try :
                                k = w[1] + w[3]
                                random_list.append(k)
                                #print(random_list)
                            except IndexError :
                                pass
                        elif typ_analysis == 6 :
                            try :
                                k = w[2] + w[3]
                                random_list.append(k)
                            except IndexError :
                                pass
                        else :
                            print('ERROR')
                    random.shuffle(random_list)
                value = 0
                sal_value = 0
                for x in range (0, len(random_list)) :
                    out_file.write(random_list[sal_value] + '\t')
                    value += 1
                    sal_value += 1
                    if value == 8 :
                        out_file.write('\n')
                        value = 0
                out_file.write('\n')
            else :
                print('ERROR --analysis_type not correct')

def montecarlo_extractor(in_file, out_file, typ_analysis) :
    random_list = []
    count_lol = 0
    with open(in_file, 'r') as f:
        for element in f :
            count_lol += 1
            x = element.split('\t')
            for sal in x : # it extracts the portion of the code corresponding to SalI boxes, which is found between round brackets
                y = sal.split('(') 
                z = y[1].split(')')
                w = y[0] + z[1]
                if typ_analysis == 0 : # SalI boxes
                    random_list.append(z[0])
                else :
                    if typ_analysis == 1 : # promoters + all SNPs
                        q = y[0] + z[1]
                        random_list.append(q)
                    elif typ_analysis == 2 : #AVA - CAT
                        q = w[1] + w[2]
                        random_list.append(q)
                    elif typ_analysis == 3 : # total repeat
                        q = w[0] + '(' + y[1]
                        random_list.append(q)
                    elif typ_analysis == 4 : # AVA + CAT + 3ETS
                        q = w[1] + w[2] + w[3]
                        random_list.append(q)
                    elif typ_analysis == 5 : # AVA + 3ETS
                        q = w[1] + w[3]
                        random_list.append(q)
                    elif typ_analysis == 6 : # CAT + 3ETS
                        q = w[2] + w[3]
                        random_list.append(q)
                    else :
                        print('ERROR. This option is not available')
            random.shuffle(random_list)
            random_list = map(lambda s: s.strip(), random_list)
        value = 0 
        sal_value = 0 
        for x in range (0, len(random_list)) : # Just the portion of the barcode corresponding to SalI boxes is written in the temporary file
            out_file.write(random_list[sal_value] + '\t')
            value += 1
            sal_value += 1 
            if value == 8 : # Up to 8 SalI boxes are written on the same line before newline
                out_file.write('\n')
                value = 0 
        out_file.write('\n')
        
        
def combi_maker(in_file_2, combi_list_in, value_round, count_in) : # it creates the SalI boxes / SNPs combinations
    with open(in_file_2, 'r') as f :
        for element in f :
            x = element.split('\t')
            sal1 = 0
            sal2 = 1 + value_round
            list_fused = []
            for sal in range(int(len(x)) - count_in) : # the combinations of SalI boxes / SNPs are created
                fused = str(x[sal1]) + '/' + str(x[sal2])
                list_fused.append(fused)
                sal1 += 1
                sal2 += 1
            for combi in list_fused : # the combinations are appended to a new list
                combi_list_in.append(combi)
    return(combi_list_in)

def recursive_count(lst, item, out): # it counts how many time a SalI box/SNP combination recurs in the pool of SalI box/SNP combinations
    count = 0
    for elem in lst:
        if elem == item:
            count += 1
        elif type(elem) in (list, dict, set, tuple):  # it outputs the combination and the number of times it was found in the pool
            count += recursive_count(elem, item)
    out.write(str(item) + ' ' + str(count) + '\n')
    return item, count

def features_counter(temporary_file, features_count) :
    if features_count == 0 :
        pass
    elif features_count == 1 :
        file_name = os.path.splitext(os.path.basename(temporary_file))[0]
        tmp_directory = os.path.dirname(os.path.realpath(temporary_file))
        tmp_locus = os.path.join(tmp_directory, file_name + '_counts' + '.txt')
        with open(temporary_file, 'r') as f :
            main_list = []
            for element in f :
                x = element.split('\t')
                x[-1] = x[-1].strip()
                y = [x for x in x if len(x.strip()) > 0]
                if not y :
                    pass    
                else :
                    for element in y :
                        main_list.append(element)
        the_count = {i:main_list.count(i) for i in main_list}
        with open(tmp_locus, 'w') as f:
            f.write(json.dumps(the_count))
    else :
        print('ERROR: invalid --features_count option')
    
parser = argparse.ArgumentParser(description='Neighbour_finder')

parser.add_argument('-i', '--input', help = 'Input file to analyze', type = str, required = True )
parser.add_argument('-t', '--temporary', help = 'temporary_file file', type = str, required = True)
parser.add_argument('-tr', '--temporary_remove', help = 'enables temporary file removal, 0 = remove, 1 = maintain. (default 0)', type = int, default = 0)
parser.add_argument('-o', '--out_combinations', help = 'combi_file', type = str, required = True)
parser.add_argument('-a', '--analysis', help = 'Type of analysis: 0 = SalI boxes, 1 = promoters + all SNPs, 2 = AvaI + CAT, 3 = total repeat, 4 = AvaI + CAT + 3ETS, 5 = AvaI + 3ETS, 6 = CAT + 3ETS', type = int, required = True)
parser.add_argument('-r', '--random_in', help = 'Enable data randomization: 0 = no randomization, 1 = randomization. (default 0)', type = int, default = 0)
parser.add_argument('-c', '--features_count', help = 'Enable features recurrence: 0 = not enabled, 1 = enabled. (dafault 0)', type = int, default = 0)
parser.add_argument('-m', '--random_montecarlo', help = 'Run a montecarlo simulation', type = int, required = False, default = 0)

args = parser.parse_args()

input_file = args.input
temporary_file = args.temporary
temporary_remove = args.temporary_remove
out_combi_file = args.out_combinations
typ_analysis = args.analysis
answer = args.random_in
features_count = args.features_count
random_montecarlo = args.random_montecarlo

output_file_write = open(temporary_file, 'w')
out_combi_file_write = open(out_combi_file, 'a')

if random_montecarlo == 0 :
    print('Standard analysis initiated')
    sal_extractor(in_file = input_file, out_file = output_file_write, random_shuffle = answer, typ_analysis = typ_analysis)
elif random_montecarlo == 1 :
    print('Montecarlo analysis initiated')
    list10k = []
    list20k = []
    list30k = []
    list_10k_total = []
    list_20k_total = []
    list_30k_total = []
    for z in range(0,1000) : # the randomization is carried out 1000 times
        output_file_write = open(temporary_file, 'w')
        montecarlo_extractor(in_file = input_file, out_file = output_file_write, typ_analysis = typ_analysis)
        output_file_write.close()
        value_round = 0
        count_in = 2
        for x in range(0,3) : # the analysis is performed for 10, 20 and 30 kb, so 3 cycles are needed
            combi_list = []
            combi_list_in = combi_maker(temporary_file, combi_list, value_round, count_in)
            if value_round == 0 : # 10kb
                counter_dictionary = Counter(combi_list_in)
                dictList = counter_dictionary.items()
                list10k.append(dictList)
            elif value_round == 1 : # 20 kb
                counter_dictionary = Counter(combi_list_in)
                dictList = counter_dictionary.items()
                list20k.append(dictList)
            elif value_round == 2 : # 30 kb
                counter_dictionary = Counter(combi_list_in)
                dictList = counter_dictionary.items()
                list30k.append(dictList)
            value_round += 1
            count_in += 1
    for element in list10k :
        for combination in element :
            list_10k_total.append(combination[0])
    for element in list20k :
        for combination in element :
            list_20k_total.append(combination[0])
    for element in list30k :
        for combination in element :
            list_30k_total.append(combination[0])
    
    counter_dictionary_10k = Counter(list_10k_total)
    dictList_10k = counter_dictionary_10k.items()
    counter_dictionary_20k = Counter(list_20k_total)
    dictList_20k = counter_dictionary_20k.items()
    counter_dictionary_30k = Counter(list_30k_total)
    dictList_30k = counter_dictionary_30k.items()
    
    out_combi_file_write.write('Distance ' + '10 kb' + '\n')
    for h in dictList_10k :
        total_recurrence = 0
        list_of_numbers = []
        for element in list10k :
            for subelement in element :
                if h[0] == subelement[0] :
                    list_of_numbers.append(int(subelement[1]))
            total_recurrence = sum(list_of_numbers)
        average = total_recurrence/float(1000) # for the aritmentic average substitute the number with float(h[1])
        out_combi_file_write.write(str(h[0]) + ' ' + str(average) + '\n')
    
    out_combi_file_write.write('Distance ' + '20 kb' + '\n')
    for h in dictList_20k :
        total_recurrence = 0
        list_of_numbers = []
        for element in list20k :
            for subelement in element :
                if h[0] == subelement[0] :
                    list_of_numbers.append(int(subelement[1]))
            total_recurrence = sum(list_of_numbers)
        average = total_recurrence/float(1000) # for the aritmentic average substitute the number with float(h[1])
        out_combi_file_write.write(str(h[0]) + ' ' + str(average) + '\n')
    
    out_combi_file_write.write('Distance ' + '30 kb' + '\n')
    for h in dictList_30k :
        total_recurrence = 0
        list_of_numbers = []
        for element in list30k :
            for subelement in element :
                if h[0] == subelement[0] :
                    list_of_numbers.append(int(subelement[1]))
            total_recurrence = sum(list_of_numbers)
        average = total_recurrence/float(1000) # for the aritmentic average substitute the number with float(h[1])
        out_combi_file_write.write(str(h[0]) + ' ' + str(average) + '\n')

else :
    print('ERROR. Random Montecarlo option must be 0 (deactivated) or 1 (activated)')

output_file_write.close()
combi_list = []

value_round = 0
count_in = 2
if random_montecarlo == 0 :
    for x in range(0,3) : # the process is repeated three times to see the association of SalI boxes/SNPs at 10, 20 and 30 kb of distance
        out_combi_file_write.write('Distance ' + str(count_in - 1) + '0 kb' + '\n')
        combi_list_in = combi_maker(temporary_file, combi_list, value_round, count_in)
        hello = dict.fromkeys(combi_list_in)
        combi_list2 = list(hello)
        for item in combi_list2 :
            item, count = recursive_count(combi_list, item, out_combi_file_write)
        hello.clear()
        del combi_list2[:]
        del combi_list_in[:] 
        value_round += 1
        count_in += 1
elif random_montecarlo != 1 :
    for x in range(0,3) :
        out_combi_file_write.write('Distance ' + str(count_in - 1) + '0 kb' + '\n')

features_counter(temporary_file = temporary_file, features_count = features_count)

if temporary_remove == 0 :
    os.remove(temporary_file)
elif temporary_remove == 1 :
    output_file_write.close()
else :
    print('ERROR: invalid --temporary_remove option')