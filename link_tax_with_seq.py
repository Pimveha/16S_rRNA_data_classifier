
# data retrieved from: https://greengenes.lbl.gov/Download/OTUs/


def link_tax_with_seq(seq_file, tax_file, output_file):
    with open(seq_file, "r") as r_seq, open(tax_file, "r") as r_tax:
        
        

if __name__ == "__main__":
    aligned_seqs =  "./raw_data/16S_rRNA_aligned.fasta"
    link_tax_with_seq(aligned_seqs)
    