from multiprocessing import Pool
import time

input_file = 'Homo_sapiens.GRCh38.dna_sm.toplevel.fasta'
list_file = 'coding.mutation.aa.sort.tsv.pylist'
output_file = 'gene.cpu50.fasta'
name_detail = True

def task(name):
    seq_file = {}
    with open(input_file, 'r') as input_fasta:
        for line in input_fasta:
            line = line.strip()
            if line[0] == '>':
                seq_id = line.split()[0]
                seq_file[seq_id] = ''
            else:
                seq_file[seq_id] += line

    input_fasta.close()

    
    list_dict = {}
    with open(list_file, 'r') as list_table:
        for line in list_table:
            if line.strip():
                line = line.strip().split('\t')
                list_dict[line[0]] = [line[1], int(line[2]) - 1, int(line[3]), line[4]]


    list_table.close()

   
    def rev(seq):
        base_trans = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N', 'a':'t', 'c':'g', 't':'a', 'g':'c', 'n':'n'}
        rev_seq = list(reversed(seq))
        rev_seq_list = [base_trans[k] for k in rev_seq]
        rev_seq = ''.join(rev_seq_list)
        return(rev_seq)


    output_fasta = open(output_file, 'w')
    for key,value in list_dict.items():
        if name_detail:
            print('>' + key, '[' + value[0], value[1] + 1, value[2], value[3] + ']', file = output_fasta)
        else:
            print('>' + key, file = output_fasta)

        seq = seq_file['>' + value[0]][value[1]:value[2]]
        if value[3] == '+':
            print(seq, file = output_fasta)
        elif value[3] == '-':
            seq = rev(seq)
            print(seq, file = output_fasta)

    output_fasta.close()

if __name__ == '__main__':
    start_time = time.time()
    pool = Pool(processes=50)  
    for i in range(50):
        pool.apply_async(task, (i,))  
    pool.close()  
    pool.join()  
    end_time = time.time()
    duration = end_time - start_time
    print('main process duration: %.3fs' % duration) 
