
# data retrieved from: https://greengenes.lbl.gov/Download/OTUs/
import numpy as np


def link_tax_with_seq(seq_file, tax_file):
    taxa_id_dict = {}
    taxa_hierarchy = ["kingdom", "phylum", "class",
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
                if any(char.isdigit() for char in taxa):
                    taxa_id_dict.pop(id)
                    continue
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

    # simple implementation if not in "ACGT-", val = 0

    base_dict = {"-": 0.0, "A": 0.2, "C": 0.4, "G": 0.6, "T": 0.8}
    kingdom_set, phylum_set = set(), set()
    for k, v in taxa_id_dict.items():
        float_list = []
        for base in v["sequence"]:
            float_list.append(base_dict.setdefault(base, 1.0))
        # print(float_list)
        # print(f"{k=}, {v['taxa']['kindom']}")
        # print(f"{v['taxa']['kindom']}")
        # kingdom_set.add(v['taxa']['phylum'])
        phylum_set.add(v['taxa']['phylum'])
    # print(kingdom_set) = {'Bacteria', 'Archaea'}
    print(phylum_set)


def count_non_overlapping(seq1, seq2):
    seq1_arr, seq2_arr = np.array(list(seq1)), np.array(list(seq2))
    diff_indices = np.where(seq1_arr == seq2_arr)[0]
    # print(diff_indices)
    # print(diff_indices.size)
    # return diff_indices.size
    return len(diff_indices)


def top_overlap(taxa_id_dict, compare_seq, top=30):
    # top_taxa_dict = {-i: {"count": 0, "percentage": 0, "taxa": ""} for i in range(30)}
    top_taxa_list = [(0, 0) for _ in range(30)]

    for key in taxa_id_dict:
        compare_me_2 = taxa_id_dict[key]["sequence"]
        # print(f"{compare_seq}\n{compare_me_2}")
        count = count_non_overlapping(compare_seq, compare_me_2)
        if count > top_taxa_list[top-1][1]:
            top_taxa_list[top-1] = (key, count)
            top_taxa_list.sort(key=lambda x: x[1], reverse=True)
            # print(count)
            print(top_taxa_list[-1][1])  # wost overlap
            print(top_taxa_list[0][1])  # best overlap
            # print(top_taxa_list)
    return top_taxa_list


def show_taxa_from_ids(taxa_id_dict, id_list):
    for id in id_list:
        # print(taxa_id_dict[id]["taxa"])
        print(taxa_id_dict[id]["taxa"].values())


def most_likely_taxa(taxa_id_dict, id_list):
    # best_option = id_list.copy()
    for taxa_name in taxa_id_dict[id_list[0]]["taxa"]:
        # cur_taxa_list = [t for t in taxa_id_dict[id_list[0]]["taxa"][taxa_name]]

        # for id in id_list:
        all_options_list = [taxa_id_dict[id]["taxa"][taxa_name]
                            for id in id_list]
        best_option = max(set(all_options_list), key=all_options_list.count)
        print("---------------------------------")
        print(f"{taxa_name=}\n{best_option=}")

def check_accuracy():
    

if __name__ == "__main__":
    aligned_seqs = "./raw_data/16S_rRNA_aligned.fasta"
    tax_file_green_genes = "./raw_data/otu_id_to_greengenes.txt"

    taxa_id_dict = link_tax_with_seq(aligned_seqs, tax_file_green_genes)
    # get_array(taxa_id_dict)
    # percent_non_overlapping("ACTTTCACCGAGA", "TCTTTCCCCGAGA")
    # top_taxa_list = top_overlap(taxa_id_dict, comp_seq)
    # print(taxa_id_dict.keys())
    # for k in taxa_id_dict:
    #     print(f"{k}:\t{taxa_id_dict[k]['taxa']['kingdom']}")
    top_taxa_list = top_overlap(
        taxa_id_dict, taxa_id_dict["574960"]["sequence"])
    # 106803
    top_30_ids = [item for item, _ in top_taxa_list]
    # show_taxa_from_ids(taxa_id_dict, top_30_ids)
    print("true taxa:", taxa_id_dict["574960"]["taxa"])
    most_likely_taxa(taxa_id_dict, top_30_ids[1:])

# way to measure acuracy:
# check for all sequences, how deep into the hierarchy does it go,
# check if guessed taxa corresponds to the true taxa
# if not, how deep into the hierarchy does it go
