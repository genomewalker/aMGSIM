# # Get codons that can get damaged
# # C -> T
import Bio.Data.CodonTable
import pandas as pd
from itertools import product
import numpy as np
import json


def get_codon_table():
    codons_table = Bio.Data.CodonTable.standard_dna_table.forward_table
    codons_stop = Bio.Data.CodonTable.standard_dna_table.stop_codons
    codons_stop = {el: "*" for el in codons_stop}
    codons_table = {**codons_table, **codons_stop}
    return codons_table


codons_table = get_codon_table()

filtered_dict = {k: v for (k, v) in codons_table.items() if "C" in k}
codon_list = list(filtered_dict.keys())

d = {
    "A": "A",
    "T": "T",
    "G": "G",
    "C": ["T", "C"],
}

d_rev = {
    "T": "T",
    "C": "C",
    "A": "A",
    "G": ["G", "A"],
}

codons_damage = {}
for k in codon_list:
    c_list = list(map("".join, product(*map(d.get, k))))
    c_list.remove(str(k))
    codons_damage[k] = {
        "aa": codons_table[str(k)],
        "mods": {el: codons_table[str(el)] for el in c_list},
    }


# A -> G
filtered_dict_rev = {k: v for (k, v) in codons_table.items() if "G" in k}
codon_list_rev = list(filtered_dict_rev.keys())
codons_damage_rev = {}
codons_damage_rev_l = ()
i = 0
for k in codon_list_rev:
    c_list = list(map("".join, product(*map(d_rev.get, k))))
    c_list.remove(str(k))
    codons_damage_rev[k] = {
        "aa": codons_table[str(k)],
        "mods": {el: codons_table[str(el)] for el in c_list},
    }
print(json.dumps(codons_damage_rev, indent=4))

n = 0
codons_damage_df = pd.DataFrame(columns=["codon", "aa", "codon_d", "aa_d"])
for k in codons_damage:
    for i in codons_damage[k]["mods"]:
        codons_damage_df.loc[n] = [
            k,
            codons_damage[k]["aa"],
            i,
            codons_damage[k]["mods"][i],
        ]
        n += 1

conditions = [
    (codons_damage_df["aa"] != codons_damage_df["aa_d"]),
    (codons_damage_df["aa"] == codons_damage_df["aa_d"]),
    (codons_damage_df["aa_d"] == "*"),
]
choices = ["nonsyn", "syn", "stop"]

codons_damage_df["class"] = np.select(conditions, choices)
codons_damage_df[["aa_d"]].groupby("aa_d").size()


n = 0
codons_damage_rev_df = pd.DataFrame(columns=["codon", "aa", "codon_d", "aa_d"])
for k in codons_damage_rev:
    for i in codons_damage_rev[k]["mods"]:
        codons_damage_rev_df.loc[n] = [
            k,
            codons_damage_rev[k]["aa"],
            i,
            codons_damage_rev[k]["mods"][i],
        ]
        n += 1
conditions = [
    (codons_damage_rev_df["aa_d"] == "*"),
    (codons_damage_rev_df["aa"] != codons_damage_rev_df["aa_d"]),
    (codons_damage_rev_df["aa"] == codons_damage_rev_df["aa_d"]),
]
choices = ["stop", "nonsyn", "syn"]
codons_damage_rev_df["class"] = np.select(conditions, choices)
codons_damage_rev_df[["class"]].groupby("class").size()

codons_damage_rev_df.to_csv(
    path_or_buf="/vol/cloud/antonio/projects/anc-com-sim/sandbox/test_data/SPAAM/codons_damage_rev_df",
    sep="\t",
    index=False,
)

codons_damage_df.to_csv(
    path_or_buf="/vol/cloud/antonio/projects/anc-com-sim/sandbox/test_data/SPAAM/codons_damage_df",
    sep="\t",
    index=False,
)

df = pd.read_csv("/tmp/text.tsv", sep="\t")
read_files = {}

d = df.groupby(["comm", "file_type"])[["pair", "file"]].apply(lambda x: dict(x.values))
for comm in df["comm"].unique():
    read_files[str(comm)] = {}
    for index, val in d.items():
        if str(comm) == str(index[0]):
            read_files[str(comm)][str(index[1])] = val
