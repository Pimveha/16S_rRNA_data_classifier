
# data retrieved from: https://greengenes.lbl.gov/Download/OTUs/


def link_tax_with_seq(seq_file, tax_file):
    taxa_id_dict = {}
    taxa_hierarchy = ["kindom", "phylum", "class",
                      "order", "family", "genus", "species"]
    id, sequence = "", ""

    with open(seq_file, "r") as r_seq, open(tax_file, "r") as r_tax:
        for i, line in enumerate(r_seq):
            if i % 2 == 0:
                id = line.rstrip()[1:]
            else:
                sequence = line.rstrip()
            taxa_id_dict[id] = {"sequence": sequence, "taxa": ""}

        for line in r_tax:
            id, taxa = line.rstrip().split("\t")
            # taxa = taxa.rstrip()
            if id in taxa_id_dict:
                taxa_list = [taxon[3:] for taxon in taxa.split(";")]
                taxa_id_dict[id]["taxa"] = dict(zip(taxa_hierarchy, taxa_list))
    # print(taxa_id_dict)
    # print(taxa_id_dict.keys())
    # for key in taxa_id_dict:
    #     print(f"{key=}\n{taxa_id_dict[key]['taxa']}\n\n")
    # len sequence: 7682
    return taxa_id_dict


def get_array(taxa_id_dict):
    # check_set = [set(taxa_id_dict[k]['sequence']) for k in taxa_id_dict.keys()]
    # unique_chars = set().union(*check_set)
    # print(unique_chars)
    # 'G', 'M', 'H', 'Y', 'W', 'S', 'A', 'N', 'C', 'R', 'K', 'T', '-', 'B', 'V', 'D'
    # https://en.wikipedia.org/wiki/Nucleotide#Abbreviation_codes_for_degenerate_bases
    
    for 


if __name__ == "__main__":
    aligned_seqs = "./raw_data/16S_rRNA_aligned.fasta"
    tax_file_green_genes = "./raw_data/otu_id_to_greengenes.txt"

    taxa_id_dict = link_tax_with_seq(aligned_seqs, tax_file_green_genes)
    get_array(taxa_id_dict)
