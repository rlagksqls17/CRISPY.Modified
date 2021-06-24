import gzip
import glob
import re
import pdb

    #DNA의 상보적인 서열을 작성해주는 함수 (‘3 -> '5 방향으로 작성)
def Reverse_Complement(dna_seq):
    Com_seq = []
    for nuc in dna_seq:
        if nuc == 'A':
            Com_seq.append('T')
        elif nuc == 'T':
            Com_seq.append('A')
        elif nuc == 'C':
            Com_seq.append('G')
        elif nuc == 'G':
            Com_seq.append('C')
        else:
            pass
    return "".join(Com_seq) #배열을 문자열로 변경 #Com_seq.reverse() #여기서 (‘5->‘3) >> (‘3->‘5)

def Clear_list(R1file, R2file):
    R1file_list = []
    R2file_list = []
    Total_list = []

    for each_fastq_file in glob.glob(R1file):
        R1file_list.append(each_fastq_file)
    for each_fastq_file in glob.glob(R2file):
        R2file_list.append(each_fastq_file)

    R1file_list.sort(key=lambda a: (int(re.search(r'\d+', a.split('_')[1]).group()), len(a)))
    R2file_list.sort(key=lambda a: (int(re.search(r'\d+', a.split('_')[1]).group()), len(a)))

    list_count = 0
    total_list = []
    for k in (R1file_list):
        total_list.append((R1file_list[list_count],R2file_list[list_count]))
        list_count += 1
    
    return total_list

#main

R1file = '/content/drive/MyDrive/workspace/*_R1_001.fastq.gz'
R2file = '/content/drive/MyDrive/workspace/*_R2_001.fastq.gz'
seq_start = 'CCGTGCCATCA'
seq_end = 'AAGAATCGCTC'
sgRNA_seq = 'ATTTGAGAACGTAGATCGG'

total_list = Clear_list(R1file, R2file)

for k in total_list:
    R1file_list = []
    R2file_list = []
    with gzip.open(str(k[0]),"r") as g:
        print(k[0])
        for line in g:
            R1file_list.append(line)
    with gzip.open(str(k[1]),"r") as a:
        print(k[1])
        for line in a:
            R2file_list.append(line)
    
        
    clear_5 = []

    for line_str in R1file_list:
        striped = ((str(line_str).replace('b', '')).replace('\\n','')).strip("'")
        clear_5.append(striped) #R2 파일 정제
# clear_5.append(Reverse_Complement(striped[::-1])) #R2 파일 정제
    clear_3 = []
    count = 0
    total_read = 0
    
    for line_str in R2file_list:
        striped = ((str(line_str).replace('b', '')).replace('\\n','')).strip("'") + clear_5[total_read] # R1 파일 정제와 R2 파일을 이어붙이기 
 #striped = ((str(line_str).replace('b', '')).replace('\\n','')).strip("'") + clear_5[count]
        count+=1
        total_read +=1
        clear_3.append(striped)
    

    real_read = 0
    for line in clear_3:
        if line.find(sgRNA_seq)>0:
            real_read+=1


    print(f"Total reads: {round(count/5)}")
    print(f"indel frequency : {round((((count/5)-real_read)/(count/5))*100,2)}%")

    print("\n")


    

