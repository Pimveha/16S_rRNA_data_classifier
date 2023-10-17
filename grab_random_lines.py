from random import shuffle

def grab_random_fasta_lines(input_fasta, shrink_fraction, write_file):
    fasta_list = []
    id, seq = "", ""
    with open(input_fasta, "r") as r_file:
        for i, line in enumerate(r_file):
            if i % 2 == 0:
                id = line
            else:
                seq = line
            fasta_list.append(f"{id}{seq}")
    stop_line_count = int(len(fasta_list) * shrink_fraction)/2
    
    shuffle(fasta_list)
    with open(write_file, "w") as w_file:
        for i, line in enumerate(fasta_list):
            w_file.write(line)
            if i >= stop_line_count:
                break
    
    
if __name__ == "__main__":
    input_fasta = "./raw_data/gg_otus_6oct2010/rep_set/gg_97_otus_6oct2010_aligned.fasta"
    shrink_fraction = 50/320
    write_file = "./raw_data/16S_rRNA_aligned.fasta"
    
    grab_random_fasta_lines(input_fasta, shrink_fraction, write_file)
    