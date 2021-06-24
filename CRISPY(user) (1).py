#CRIS.py (fastqjoin.colab)
from __future__ import division
#‘/’ 연산자를 모듈 전체에서 실제 나누기를 의미하도록 변경
import os
# 운영체제를 제어하는 모듈을 사용 (Windows로 파일과 폴더를 만들고 복사 가능)
import glob
# 파일들의 리스트를 뽑을 때 사용 (현재 경로의 모든 파일을 뽑기 가능)
# CRISPY의 핵심 기능임
import sys
# 파이썬 인터프리터를 제어할 수 있음 (경로 제어)
import csv
# csv 모듈은 csv 형식의 표 형식 데이터를 읽고 쓰는 클래스를 구현
from collections import Counter
# 컨테이너에 동일한 값의 자료가 몇 개인지를 파악하는데 사용하는 객체
from collections import OrderedDict
# 입력된 item 들의 순서를 기억하는 Dictionary 클래스
import pandas as pd 
# 데이터 분석을 위한 전용 라이브러리임 
# 왜 pandas와 csv 모듈을 같이 쓸까??


#CRIS.py
#
#Modify parameters below in the get_parameters() section.

def get_parameters(): #사용자가 입력한 정보를 반환한다.
    Data_Frame = pd.read_excel("/content/drive/MyDrive/workspace/sample1.xlsx") 
    ID = str.upper(Data_Frame.loc[0,'information'])
    ref_seq = str.upper(Data_Frame.loc[1,'information'])
    seq_start = str.upper(Data_Frame.loc[2,'information'])
    seq_end = str.upper(Data_Frame.loc[3,'information'])
    fastq_files = '/content/drive/MyDrive/workspace/*.fastqjoin'
    sgRNA_list = [Data_Frame.columns[2],str.upper(Data_Frame.loc[4, 'sg1'])]
    return ID,ref_seq,seq_start,seq_end,fastq_files,sgRNA_list #입력한 정보들을 반환
    

def pairwise(iterable):
    #Make an ordered dictionary from items in in test list
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable) #값을 차례대로 꺼낼 수 있는 객체로 만들어줌
    return zip(a, a) #튜플형태로 묶어줌


def make_project_directory(save_dir):
    #Create a project directory to store files in
    cwd = os.getcwd()   #Get current working directory
    try:
        os.stat(save_dir) #경로에 해당하는 정보를 불러옴
    except:
        os.mkdir(save_dir) #또는 해당 폴더를 만듬

def write_to_file(record_entry,f):
    #Writes results to the results_counter.txt file.  Method of inserting space between each well to make the .txt file easier to read.  Based on # of items in well list
    f.write("\n")
    for record in record_entry:         #record entry is a giant list (master_Record) that contains well-lists as elements
        for line in record:
            if len(record) > 2:
                f.write(str(line) + '\n')
            else:
                pass
        if len(line) > 2:
            f.write('\n\n') 
        else:
            pass        

def make_counter(indel_size_list, current_fastq_file, dict_Counters, c_Counter, SNP_test, raw_wt_counter):
    top_common = 12             #Number of top found reads to write to results_counter .txt file
    temp_counter = Counter(indel_size_list).most_common(top_common)    #Count top indels present in fastq file
    temp_dict =OrderedDict()                                           
    #If there are not a total of at least 'top_common' common reads from the summary file, fill the spaces with NA, NA
    if len(Counter(indel_size_list).most_common(top_common)) <top_common:
            for i in range (0,top_common - len(Counter(indel_size_list).most_common(top_common))):
                temp_counter.append(("NA","NA"))
            else:
                pass
    temp_dict['Name']=str(current_fastq_file)              #Fill in Name column of CSV file with the file fastq name
    temp_dict['Sample']=''                                 #Makes a blank column for user to fill in sample names
    counter =1
    
    for k in dict_Counters:
        try:
            temp_dict[k] = str( str(dict_Counters[k]) + ' (' + str(format((dict_Counters[k] / c_Counter*100), '.1f')) + '%)')
            temp_dict[str('%'+ k)] = str(format((dict_Counters[k] / c_Counter)*100, '.1f'))
            temp_dict['Total'] = c_Counter
            temp_dict['SNP_test'] = SNP_test
            temp_dict['raw_wt_counter'] = raw_wt_counter
            temp_dict['Total_indel'] = str((format(((c_Counter - indel_size_list.count(0))/c_Counter)*100,'.1f')))
        except ZeroDivisionError:
            pass
        
    for k,g in temp_counter:        #Fill in top indels sizes and amount
        temp_dict['#'+ str(counter)+ '-Indel']=k
        try:
          # temp_dict['z%Indel'+str(counter)]= format(g/c_Counter *100, '.1f')
           temp_dict['#'+ str(counter)+'-Reads(%)'] = str(g) + ' ('+ str(format(g/c_Counter*100, '.1f')) + '%)'
        except TypeError:
            pass
        counter+=1
    return temp_dict
    


def search_fastq(ID,ref_seq,seq_start,seq_end,fastq_files,test_list):
    #fastq 파일을 처리하고 그 안에 담긴 내용을 읽는다.
    test_dict=OrderedDict()                  #검색 중인 각 항목의 순서를 저장하기 위한 변수         
    dict_Counters = OrderedDict()            #검색중인 사전의 각 항목 수를 저장하기 위한 변수 (예 : current_fastq_file)
    master_distance_and_count_summary =[]    #리스트
    save_dir = os.getcwd()+"/"+ str(ID) + "/"  #현재 자신의 디렉터리 위치를 돌려준다. 
    # a = os.getcwd() a : 'C:\\Users\\joo\\Python_study\\mod'
    fastq_counter = 0               #fastq_files의 리드수를 세기 위한 카운터 변수
    master_Record = []             #Master list, Master list contains lists
    indel_size_list = []             #list of indel sizes found in the fastq.  For each read, one indel size is recorded
    csv_summary_df = pd.DataFrame() #DataFrame 형을 만들고 초기화 
    master_distance_and_count_summary =[]

    for x,y in pairwise(test_list):          #Create an ordered dictionary of items in test_list
        test_dict[x] = y
    save_dir = os.getcwd()+"/"+ str(ID) + "/" #계산 후 현재 위치를 저장 
    print("Working directory: {}".format(str(os.getcwd())))#현재 작동되는 디렉토리를 출력
    make_project_directory(save_dir)
    file_name = save_dir+'results_counter ' + ID + '.txt'
    f = open(file_name, "w") #파일 이름을 다음과 같이 지정하고 저장 (쓰기 모드)
    wt_distance = ref_seq.find(seq_end)+len(seq_end) - ref_seq.find(seq_start) #Expected size of WT read, the distance between the two anchor points made from seq_start and seq_end
    f.write(ID + '\n') #문서에 ID 쓰고
    f.write(str("seq_start: "+seq_start+'\n')) #시작서열과   
    f.write(str("seq_end: "+seq_end+'\n')) #끝서열을 씀
    f.write("Test_Sequences: \n") #그리고 나서 
    for key, value in test_dict.items():  #Go through the test_dict and write each item that is being searched for
        f.write(str(key)+": "+value+'\n') #딕셔너리에 존재하는 요소들을 하나씩 뽑아내서 작성
        print(key, value) 
        dict_Counters[str(key)]=0 #딕셔너리에 존재하는 가이드 RNA를 dict_Conters에 딕셔너리 형으로 저장
        
    print('Expected WT distance: {}'.format(wt_distance))     #The expected distance between  seq_start and seq_end if the DNA is WT/ REF
    if wt_distance < 0: #이후 계산할 서열을 계수하고, 0 보다 작을 시 에러값 출력
        f.write("\n\n WARNING: THIS IS NOT GOING TO GIVE YOU THE FULL DATA. YOUR EXPECTED WT DISTANCE IS LESS THAN 0, it is: {}\n  Check your seq_start and seq_end again\n".format(wt_distance))

    else:
        pass

    print("Program Running")

    for each_fastq_file in glob.glob(fastq_files):   #지원하는 모든 파일을 리스트 형식으로 반환(*와 ?만 가능)해서 for문 돌림
        c_Counter = 0                     #Reset control counter to 0, this counter counts how many times both seq_start and seq_end are found in a line.
        start_counter = 0                  #How many times seq_start is found in a fastq file, used for SNP check
        end_counter = 0                    #How many times seq_end is found in a fastq file, used for SNP check
        raw_wt_counter=0                   #Calculate how many times first item in test_list is found in fastqfile, a check for SNPs
        for items in test_dict:            #RESET DICT COUNTERS TO 0
            dict_Counters[str(items)]=0
        indel_size_list = []                          #List of all the indel sizes found in each file
        fastq_name = str(each_fastq_file)                                #Well name that contains the clone being screened
        line_list =[]                                  #List of all the lines found in a fastq file.  These lines must contain both seq_start and seq_end.  Used to find most highly occuring sequences
        current_fastq_file = open(str(each_fastq_file), "r")     
        for line in current_fastq_file:   
            if list(test_dict.items())[0][1] in line:           #각각의 fastq 파일을 검사하여 첫번째 가이드 RNA가 있는 횟수를 센다.
                raw_wt_counter+=1
            if line.find(seq_start)>0 and line.find(seq_end)>0:  
                c_Counter += 1
                start_counter +=1     #Count # of times seq_start is found
                end_counter +=1       #Count # of times seq_end is found
                read_start = line.find(seq_start)
                read_end = line.find(seq_end)+len(seq_end) #실제 파일과 야생형 파일을 비교해서 indel_size를 계산
                indel_size = line.find(seq_end)+len(seq_end) - line.find(seq_start) - wt_distance 
                indel_size_list.append(indel_size) #이를 indel_size list에 append 시킴
                line_list.append(line[read_start:(read_end)])
                for item in test_dict:
                    if test_dict[item] in line:
                        dict_Counters[item] +=1 #line안에 sgRNA가 있으면 개수를 센다.
                    else:
                        pass
            elif line.find(seq_start)>0 and line.find(seq_end)<0:        #If seq_start is found and seq_end is not found, for SNPT test
                start_counter +=1
            elif line.find(seq_end)>0 and line.find(seq_start)<0:        #If seq_end is found and seq_start is not found, for SNP test
                end_counter+=1
            else:
                pass
        current_fastq_file.close()               
        try:
            SNP_test = format(start_counter / end_counter, '.2f')        #Compare counts of seq_start and seq_end for SNP test
        except ZeroDivisionError:
            pass
        try:
            raw_wt_counter = str(str(raw_wt_counter) +  ' ('+ format((raw_wt_counter/ dict_Counters[list(test_dict.items())[0][0]]), '.1f') +')') #Calculate raw_wt_counter
        except ZeroDivisionError:
            pass
       
        if c_Counter == 0 : #counter가 아예 없을 경우에는 패스 함
            pass
        elif c_Counter > 10 :  #if more than 10 control read counts, record data
            print("{}: Total_reads:{}, {}".format(fastq_name,str(c_Counter).ljust(2), dict_Counters.items()))
            fastq_counter += 1
            test_list_string=str(" Testing: ")
            for k,v in dict_Counters.items(): # 딕셔너리에 존재하는 sgRNA서열과 서열의 리드수를 읽는다.  
                test_list_string=test_list_string+"({}:{}), ".format(k,v)
            temp = Counter(line_list).most_common(12) #최빈값 구하기 함수
            #summary_line is a list, format:  Miller-Plate13-C01 TOTAL:2072 OrderedDict([('g3', 2010), ('Block_Only', 0), ('Mod_Only', 2), ('Block_Mod', 0), ('Full', 0)])       [(0, 2070), (-1, 2)]
            summary_line = ([str(fastq_name) + " TOTAL:" + str(c_Counter)+" "+test_list_string+"     "+"Top_reads:"+ str(Counter(indel_size_list).most_common(12))])
            for k,v in temp:          #Append the top found DNA sequences to summary_line
                summary_line.append('{} , {}'.format(k,v))
            master_Record.append(summary_line)
            master_distance_and_count_summary.append(make_counter(indel_size_list,str(fastq_name), dict_Counters, c_Counter,SNP_test,raw_wt_counter))
        else:
             pass


    print("SUMMARY")
    make_project_directory(ID)
    #print master_distance_and_count_summary
    pd_columns = ['Name','Sample','Total', 'Total_indel', '#1-Indel','#1-Reads(%)','#2-Indel','#2-Reads(%)','#3-Indel','#3-Reads(%)','#4-Indel','#4-Reads(%)','#5-Indel','#5-Reads(%)',
             '#6-Indel','#6-Reads(%)','#7-Indel','#7-Reads(%)','#8-Indel','#8-Reads(%)', 'SNP_test', 'raw_wt_counter']

    flip_dict = list(test_dict.items())   #Need to flip order of items in dictionary.  That way when they are inserted into the excel list, the order will come out correct
    flip_dict.reverse()
    flip_dict=OrderedDict(flip_dict)
    for k, v in flip_dict.items():        #Insert the items from test_dict into position 3 for the column output, after 'Total'
        pd_columns.insert(3,k)
    for k, v in test_dict.items():        #Insert the % values for test_dic at the end
        pd_columns.append(str('%'+k))
    csv_summary_df = pd.DataFrame.from_records(master_distance_and_count_summary, index='Name', columns=pd_columns)
    csv_summary_df=csv_summary_df.sort_index()
    csv_summary_df = csv_summary_df[pd.notnull(csv_summary_df['Total'])]
    try:
        csv_summary_df.to_csv(str(save_dir+ID)+'.csv')    #Filename to save csv as
    except (IOError):
        print('ERROR.  Script did not execute properly.')
        print('The requested .csv file {} is either open or you do not have access to it.  If open, please close file and rerun program').format(str(save_dir+ID)+'.csv')
        
    master_Record = sorted(master_Record)
    print("Total wells with product:", fastq_counter)
    write_to_file(master_Record,f)
    f.close()
    

def main():
    ID = '' # DNA의 이름 변수 초기화
    ref_seq = '' # 전체 서열 변수 초기화
    seq_start = '' # 전체 서열이 사작되는 서열 초기화
    seq_end = '' # 전체 서열이 끝나는 서열 초기화
    fastq_files = '' # 분석하고자 하는 파일 초기화 (이미 fastq 파일만 타겟일까?)
    test_list = [] # 리스트 초기화
    print("CRIS.py \nMain program") #문자열 출력
    print(f"현재 디렉터리 위치를 확인하시고 그 장소에 검사 대상 파일을 넣으세요\n")
    print(f"현재 디렉터리 위치 : {os.getcwd()}\n 확인 후 엔터")
    input()
    ID, ref_seq, seq_start, seq_end, fastq_files, test_list = get_parameters() # 다음 함수 이용해서 정보를 받음
    search_fastq(ID, ref_seq, seq_start, seq_end, fastq_files, test_list) # 받은 정보를 토대로 search_fastq 함수를 실행
    print("Done") # 완료시 Done 출력하고 프로그램 종료
if __name__== "__main__":
    """
    해당 모듈이 임포트된 경우가 아니라 인터프리터에서 직접 실행된
    경우에만, if 문 이하의 코드를 돌리라는 명령
    
    인터프리터에서 직접 실행하면, __name__ 변수에 __main__이 담겨서 프린트 됨
    """
    main() #메인 함수의 실행 (가장 먼저 실행)