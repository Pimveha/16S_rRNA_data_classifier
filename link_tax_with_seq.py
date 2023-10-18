
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
                # print(id)
                # if not id.isnumeric():
                #     print(id)
                #     break
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
    #     print(f"{key=}")
    # len sequence: 7682


if __name__ == "__main__":
    aligned_seqs = "./raw_data/16S_rRNA_aligned.fasta"
    tax_file_green_genes = "./raw_data/otu_id_to_greengenes.txt"

    link_tax_with_seq(aligned_seqs, tax_file_green_genes)
