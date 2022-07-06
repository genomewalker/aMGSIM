from Bio import SeqIO, Seq
from bisect import bisect_left
import re
import pyranges as pr
import Bio.Data.CodonTable
import pandas as pd
import logging
import gzip
from mimetypes import guess_type
from functools import partial
from pathlib import Path
from pandarallel import pandarallel


def get_codon_table():
    codons_table = Bio.Data.CodonTable.standard_dna_table.forward_table
    codons_stop = Bio.Data.CodonTable.standard_dna_table.stop_codons
    codons_stop = {el: "*" for el in codons_stop}
    codons_table = {**codons_table, **codons_stop}
    return codons_table


def get_headers(records, pattern):
    """
    Get information of the reads encoded in the headers
    """
    for i, record in records:
        m1 = re.match(pattern, record.id)
        reads = {
            "Chromosome": m1.group(2),
            "Start": int(m1.group(6)),
            "End": m1.group(7),
            "Name": m1.group(0),
            "Strand": m1.group(5),
            "Start_read": int(m1.group(6)),
            "End_read": int(m1.group(7)),
            "Strand_read": m1.group(5),
            "type": m1.group(4),
            "read_length": m1.group(8),
            "damage": m1.group(9),
        }
        yield reads


def get_prodigal_coords(x):
    """
    Get prodigal coordinates from headers
    """
    # get coords from prodigal fasta header
    s = re.split("\#|\s", x.replace(" ", ""))
    coords = [int(i) for i in s[1:4]]
    return pd.Series(coords)


def fasta_to_dataframe(
    infile, header_sep=None, key="name", seqkey="sequence", feature="CDS"
):
    """Get fasta proteins into dataframe"""
    recs = SeqIO.parse(infile, "fasta")
    keys = [key, seqkey, "description"]
    data = [(r.name, str(r.seq), str(r.description)) for r in recs]
    df = pd.DataFrame(data, columns=(keys))
    df["type"] = feature
    # fix bad names
    if header_sep not in ["", None]:
        df[key] = df[key].parallel_apply(lambda x: x.split(header_sep)[0], 1)
    # df[key] = df[key].str.replace('|','_')
    return df


def get_gene_coordinates(fna):
    """
    Get gene coordinates from the genes
    """
    df_nt = fasta_to_dataframe(fna)

    # df_aa = fasta_to_dataframe(faa)
    # df_genes = df_nt.merge(df_aa, on=["name", "description", "type"])
    df_nt[["Start", "End", "Strand"]] = df_nt.description.parallel_apply(
        get_prodigal_coords, 1
    )
    # df_genes = df_genes.rename(columns={"sequence_x": "nt_seq", "sequence_y": "aa_seq"})
    df_nt["Strand"] = df_nt["Strand"].astype(str).str.replace("-1", "-")
    df_nt["Strand"] = df_nt["Strand"].astype(str).str.replace("1", "+")
    df_nt["Chromosome"] = df_nt["name"].str.extract(r"(\S+)_\d+")
    df_nt.drop(["description", "sequence"], 1)
    df_nt = df_nt[["name", "type", "Chromosome", "Start", "End", "Strand"]]
    df_nt = df_nt.rename(columns={"name": "gene_name"})
    df_nt = pr.PyRanges(df_nt)
    return df_nt


def revcomp(x, var):
    if x["Strand_read"] == "-":
        s = Seq.Seq(x[var])
        s = str(s.reverse_complement())
    else:
        s = x[var]
    return s


def get_intersections(reads, genes, genome, min_len=0):
    """
    Get intersections between reads and genes
    Returns the intersection positions and the sequence contained in the intersection
    """
    # j = df_genes.join(df_reads, report_overlap=True)
    gr = reads.intersect(genes, strandedness=False, how="first")
    # j = df_reads.overlap(df_genes, invert=True)
    gr.intersect_length = gr.lengths()
    # print(j.df[j.df['intersect_length'].astype(int) != j.df['length'].astype(int)])
    gr = pr.PyRanges(gr.df.drop(["Strand"], 1))
    # TODO Optimize retrieval of fasta sequences
    gr.intersect_seq = pr.get_fasta(gr, genome[0])
    gr.intersect_seq = gr.df.parallel_apply(revcomp, var="intersect_seq", axis=1)
    gr = pr.PyRanges(gr.df[gr.df["intersect_length"] >= min_len])
    return gr


def get_nondamaged_seqs(intersections, genome):
    gr = intersections.df[["Chromosome", "Start_read", "End_read"]]
    gr = gr.rename(columns={"Start_read": "Start", "End_read": "End"})
    # gr['Start'] = gr.Start.apply(lambda x: int(x) + 1)
    gr = pr.PyRanges(gr)
    intersections.nondamaged_seq = pr.get_fasta(gr, genome[0])
    nondamaged_seq = intersections.df.parallel_apply(
        revcomp, var="nondamaged_seq", axis=1
    )
    return nondamaged_seq


# Get damaged sequence
def calc_start(x):
    if x["Strand_read"] == "+":
        start = x["Start"] - x["Start_read"]
        return start
    else:
        start = abs(x["End"] - x["End_read"])
        return start


def calc_end(x):
    if x["Strand_read"] == "+":
        end = int(x["read_length"]) - abs(x["End"] - x["End_read"])
        return end
    else:
        end = int(x["intersect_length"]) + abs(x["End"] - x["End_read"])
        return end


def get_damaged_seqs(intersections, deam_file):
    reads_seqs = fasta_to_dataframe(deam_file, feature="read")
    reads_seqs = reads_seqs[["name", "sequence"]]
    reads_seqs = reads_seqs.rename(
        columns={"name": "Name", "sequence": "read_sequence"}
    )
    reads_seq = intersections.df.merge(reads_seqs)

    reads_seq.loc[:, "diff_start"] = reads_seq.parallel_apply(calc_start, axis=1)
    reads_seq.loc[:, "diff_end"] = reads_seq.parallel_apply(calc_end, axis=1)

    damaged_seq = pd.DataFrame(
        reads_seq.parallel_apply(
            lambda x: {
                "Name": x["Name"],
                "damaged_seq": x["read_sequence"][x["diff_start"] : x["diff_end"]],
                "damaged_seq_length": len(
                    x["read_sequence"][x["diff_start"] : x["diff_end"]]
                ),
            },
            axis=1,
        ).tolist()
    )
    damaged_seq = intersections.df.merge(damaged_seq).drop_duplicates()
    return pr.PyRanges(damaged_seq)


def take_closest(positions, position):
    """
    Assumes positions is sorted. Returns closest value to position.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(positions, position)
    if pos == 0:
        return positions[0]
    if pos == len(positions):
        return positions[-1]
    before = positions[pos - 1]
    after = positions[pos]
    if after - position < position - before:
        return after
    else:
        return before


def find_damage(x):
    # Get starting position
    diff_start = int(x["Start"]) - int(x["Start_read"]) + 1
    # Get last position
    rl = int(x["read_length"])
    diff_end = rl - (int(x["End"]) - int(x["End_read"]))
    dam_orig = x["damage"]
    # if diff_start <= number <= diff_end:
    # get damage values
    if x["damage"] == "None":
        damage = "None"
    else:
        damage = x["damage"].split(",")
        c2t = [i for i in damage if int(i) > 0]

        # Convert a2g to fwd coordinates
        a2g = [i for i in damage if int(i) < 0]
        a2g_conv = {i: (int(i) + rl + 1) for i in a2g}
        # filter the positions of the left
        c2t = [i for i in c2t if (diff_start <= int(i) <= diff_end)]
        a2g = [k for k in a2g_conv if diff_start <= int(a2g_conv[k]) <= diff_end]
        damage = c2t + a2g
        if damage:
            damage = (",").join(damage)
        else:
            damage = "None"
    # print("orig:{} :: damage:{} :: start:{} :: end:{}".format(dam_orig,damage, diff_start, diff_end))
    return damage


def get_seqs_inframe(x):
    # print(x.to_dict())
    start = x["Start"]
    end = x["End"]
    start_read = x["Start_read"]
    end_read = x["End_read"]
    start_gene = x["Start_gene"]
    end_gene = x["End_gene"]
    intersect_seq = x["intersect_seq"]
    nondamaged_seq = x["nondamaged_seq"]
    damaged_seq = x["damaged_seq"]
    sequence = x["sequence"]
    strand_read = x["Strand_read"]
    strand_gene = x["Strand_gene"]
    damage = x["damage_intersection"]

    if damage != "None":
        damage = [int(i) for i in damage.split(",")]

        if strand_read == "-" and strand_gene == "+":
            intersect_seq = str(Seq.Seq(intersect_seq).reverse_complement())
            damaged_seq = str(Seq.Seq(damaged_seq).reverse_complement())
            damage = [int(i) * -1 for i in damage]
            # damage = [int(i) * -1 for i in damage]
        elif strand_read == "+" and strand_gene == "-":
            intersect_seq = str(Seq.Seq(intersect_seq).reverse_complement())
            damaged_seq = str(Seq.Seq(damaged_seq).reverse_complement())
            damage = [int(i) * -1 for i in damage]
        # elif strand_read == '-' and strand_gene == '-':
        #     if damage != 'None':
        #         damage = [int(i) * -1 for i in damage ]

        intersect_seq_len = len(intersect_seq)
        damaged_seq_len = len(damaged_seq)
        read_len = len(nondamaged_seq)

        seq_coords = re.search(intersect_seq, sequence).span()
        coord_start = seq_coords[0]
        coord_stop = seq_coords[1]
        gene_len = len(sequence)

        pos = range(coord_start, coord_stop)
        start_codon_pos = [
            i[1] for i in enumerate(range(0, gene_len), 0) if i[0] % 3 == 0
        ]
        start_codon = take_closest(positions=start_codon_pos, position=coord_start)
        stop_codon = take_closest(positions=start_codon_pos, position=coord_stop)

        if start_codon < coord_start:
            start_codon = start_codon_pos[start_codon_pos.index(start_codon) + 1]

        diff_start = start_codon - coord_start
        end_codon_pos = max(
            [
                i[1]
                for i in enumerate(range(diff_start, intersect_seq_len), 0)
                if i[0] % 3 == 0
            ]
        )

        intersect_seq_inframe_nt = intersect_seq[diff_start:end_codon_pos]
        pos = re.search(intersect_seq_inframe_nt, intersect_seq).span()
        damaged_seq_inframe_nt = damaged_seq[pos[0] : pos[1]]

        # Translate
        intersect_seq_inframe_aa = str(Seq.Seq(intersect_seq_inframe_nt).translate())
        damaged_seq_inframe_aa = str(Seq.Seq(damaged_seq_inframe_nt).translate())
        sequence_aa = str(Seq.Seq(sequence).translate())

        # find codons with damage
        # Transform damage positions to positions in the coding fragment

        if strand_read == "-":
            cstart = abs(end_read - end) + diff_start
        else:
            cstart = abs(start_read - start) + diff_start

        c2t = [i - 1 for i in damage if int(i) > 0]
        # gene_int = range(diff_start, end_codon_pos)
        # # Convert a2g to fwd coordinates
        a2g = [i for i in damage if int(i) < 0]

        # CDS coords are 0-based
        # Damaged is 1-based
        # do we have damage
        if a2g:
            a2g_conv = {i: (int(i) + read_len + 1) - 1 for i in a2g}
            a2g = [
                (a2g_conv[k] - cstart)
                for k in a2g_conv
                if cstart <= int(a2g_conv[k]) < end_codon_pos
            ]
        else:
            a2g = []
        if c2t:
            c2t = [(i - cstart) for i in c2t if (cstart <= int(i) < end_codon_pos)]
        else:
            c2t = []

        damage_pos = c2t + a2g

        # Which codons
        codons_table = get_codon_table()
        codons = []
        for i in damage_pos:
            codon = i // 3 * 3
            codons.append(
                "{}:{}>{}:{}>{}".format(
                    codon,
                    intersect_seq_inframe_nt[codon : codon + 3],
                    damaged_seq_inframe_nt[codon : codon + 3],
                    codons_table[intersect_seq_inframe_nt[codon : codon + 3]],
                    codons_table[damaged_seq_inframe_nt[codon : codon + 3]],
                )
            )
        x["damage_codon_diffs"] = ",".join(codons)
        x["damage_inframe_ag"] = ",".join(map(str, a2g))
        x["damage_inframe_ct"] = ",".join(map(str, c2t))

        if intersect_seq_inframe_aa == damaged_seq_inframe_aa:
            x["damage_aaseq_diffs"] = None
        else:
            positions = list(
                map(
                    str,
                    [
                        i
                        for i in range(len(intersect_seq_inframe_aa))
                        if intersect_seq_inframe_aa[i] != damaged_seq_inframe_aa[i]
                    ],
                )
            )
            diff = list(
                map(
                    str,
                    [
                        intersect_seq_inframe_aa[i] + ">" + damaged_seq_inframe_aa[i]
                        for i in range(len(intersect_seq_inframe_aa))
                        if intersect_seq_inframe_aa[i] != damaged_seq_inframe_aa[i]
                    ],
                )
            )
            x["damage_aaseq_diffs"] = ",".join(
                [x + ":" + y for x, y in zip(positions, diff)]
            )

        x["intersect_seq_inframe_nt"] = intersect_seq_inframe_nt
        x["damaged_seq_inframe_nt"] = damaged_seq_inframe_nt
        x["intersect_seq_inframe_aa"] = intersect_seq_inframe_aa
        x["damaged_seq_inframe_aa"] = damaged_seq_inframe_aa
        x["sequence_aa"] = sequence_aa

    return x


log = logging.getLogger("my_logger")


def analyze_proteins(x, files, gene_predictions, min_len, outdir, debug, nproc):
    pandarallel.initialize(nb_workers=nproc, progress_bar=False, verbose=False)

    comm = x[0]
    genome = x[1]

    deam_file = files["reads"][comm]["deamSim_ofile"]["SR"]
    fragsim_file = files["reads"][comm]["fragSim_ofile"]["SR"]
    fasta_file = [s for s in files["genomes"] if genome in s]

    # Genes
    faa = gene_predictions[genome]["faa"]
    fna = gene_predictions[genome]["fna"]
    translation_table = gene_predictions[genome]["translation_table"]
    if Path(deam_file).exists():
        recs = []
        encoding = guess_type(deam_file)[1]  # uses file extension
        _open = partial(gzip.open, mode="rt") if encoding == "gzip" else open

        pattern = re.compile(
            "(\S+)___(\S+)---(\d+):(\S+):([+\-]):(\d+):(\d+):(\d+):(.*)"
        )
        with _open(deam_file) as f:
            records = enumerate(SeqIO.parse(f, "fasta"))
            recs = get_headers(records=records, pattern=pattern)
            df = pd.DataFrame(recs)
        if debug:
            log.info("Reading reads...")

        # Get reads
        df_reads = pr.PyRanges(df)

        if df_reads.df["Chromosome"].str.contains(genome).any():
            if debug:
                log.info("Reading genes...")
            # Get genes info
            genes = get_gene_coordinates(fna=fna)
            # Get intersections
            # Will rop intersects that are too short for coding AAs
            if debug:
                log.info("Finding read-gene intersections...")

            r2g_intersections = get_intersections(
                reads=df_reads, genes=genes, genome=fasta_file, min_len=min_len
            )
            names = {
                "Start_b": "Start_gene",
                "End_b": "End_gene",
                "Strand": "Strand_gene",
            }
            r2g_intersections = r2g_intersections.join(
                genes.drop("type"), report_overlap=True
            )
            r2g_intersections = pr.PyRanges(
                r2g_intersections.df[
                    r2g_intersections.df["Overlap"]
                    == r2g_intersections.df["intersect_length"]
                ].rename(columns=names)
            )
            r2g_intersections_df = r2g_intersections.df
            gn = r2g_intersections_df["Name"]
            read_multi_span = r2g_intersections_df[
                gn.isin(gn[gn.duplicated()])
            ].sort_values("Name")
            # len(read_multi_span.index)
            if debug:
                log.info("Retrieving non-damaged sequences from intersections...")
            # Get sequence from interval
            r2g_intersections.nondamaged_seq = get_nondamaged_seqs(
                intersections=r2g_intersections, genome=fasta_file
            )
            if debug:
                log.info("Retrieving damaged sequences from intersections...")
            r2g_intersections = get_damaged_seqs(
                intersections=r2g_intersections, deam_file=deam_file
            )

            r2g_intersections.damage_intersection = r2g_intersections.df.parallel_apply(
                find_damage, axis=1
            )
            if debug:
                log.info("Retrieving damaged codons from genes...")
            aa_damage = r2g_intersections.df.merge(
                fasta_to_dataframe(fna)[["name", "sequence"]].rename(
                    columns={"name": "gene_name"}
                )
            ).parallel_apply(get_seqs_inframe, axis=1)
            aa_damage = aa_damage.reset_index()
            aa_damage["community"] = comm

            # aa_damage = aa_damage[~aa_damage.damage_codon_diffs.isnull()].reset_index()
            aa_damage = aa_damage.rename(
                columns={"Name": "read_name", "Chromosome": "chromosome"}
            )

            out_suffix = ".tsv.gz"

            fname = f"{comm}---{genome}_aa-damage{out_suffix}"
            outfile = Path(outdir, fname)
            columns = [
                "community",
                "read_name",
                "gene_name",
                "chromosome",
                "damage",
                "damage_intersection",
                "damage_inframe_ct",
                "damage_inframe_ag",
                "damage_aaseq_diffs",
                "damage_codon_diffs",
                "intersect_seq",
                "intersect_seq_inframe_aa",
                "intersect_seq_inframe_nt",
                "damaged_seq",
                "damaged_seq_inframe_aa",
                "damaged_seq_inframe_nt",
            ]
            aa_damage.to_csv(
                path_or_buf=outfile,
                sep="\t",
                header=True,
                index=False,
                compression="gzip",
            )
            return outfile
        else:
            if debug:
                print("Genome {} not found in sample".format(genome))
