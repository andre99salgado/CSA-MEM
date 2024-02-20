from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from random import seed
from random import randint
import sys

file_in = '_viruses_completos.fasta'
file_out = '_viruses/_viruses_rodados_'
seed(1)

def move_first_char_to_end(text, number_of_chars):
    if len(text) > 0:
        first_char = text[0:number_of_chars-1]
        rest_of_text = text[number_of_chars-1:]
        modified_text = rest_of_text + first_char
        print("aqui")
        return modified_text
    else:
        return text
    
for x in range(50):
    file_out = '_viruses/_viruses_rodados_' + str(x) + '.fasta'
    
    array = [randint(2, 10000),randint(2, 10000),randint(2, 10000),randint(2, 10000),randint(2, 10000),randint(2, 10000),randint(2, 10000),randint(2, 10000),randint(2, 10000),randint(2, 10000),randint(2, 10000),randint(2, 10000)]
    #array = [33,29,31,391,29,454,31,28,609,30,30,28,]

    with open(file_out, 'w') as f_out:
        text = ''
        count = 0
        f = open(file_out, "a")
        for seq_record in SeqIO.parse(open(file_in, mode='r'), 'fasta'):
            # remove .id from .description record (remove all before first space)
            seq_record.description = ' '.join(seq_record.description.split()[1:])
            print('SequenceID = ' + seq_record.id)
            print('Description = ' + seq_record.description + '\n')
            print('------')
            textoo = '>' + seq_record.id + '\n' + move_first_char_to_end(seq_record.seq, array[count]) + '\n'
            f.write(str(textoo))
            print('Text2= ' + move_first_char_to_end(seq_record.seq, array[count])[0:40])
            count += 1
        

