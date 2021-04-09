#!/usr/bin/env python

"""
ancient-reads: Simulate ancient reads for each taxon in each synthetic
               community

Usage:
  ancient-reads [options] <genome_table> <config>
  ancient-reads -h | --help
  ancient-reads --version

Options:
  <genome_table>      Taxon genome info.
  <config>            Config parameters
  -h --help           Show this screen.
  -d --debug          Debug mode (no subprocesses; verbose output)
  --version           Show version.

Description:
  Simulating ancient reads for each taxon in each synthetic community

  genome_table
  ------------
  * tab-delimited
  * must contain 2 columns
    * "Taxon" = taxon name
    * "Fasta" = genome fasta file path
  * other columns are allowed

  config
  ------
  * YAML with config parameters 
    (Check https://github.com/genomewalker/aMGSIM for details)

  Output
  ------
  * A set of read files for each sample (fragSim, deamSim, ART)
    - read sequences are named by the taxon they originate from
    - directory structure: OUTPUT_DIR/COMMUNITY/ancient-read_files
  * A JSON file with the location of each file type (fragSim, deamSim, ART)
"""

# import
# batteries
from docopt import docopt
import sys
import os
import re
import logging
from functools import partial
from multiprocessing import Pool
from shutil import rmtree

# application
from aMGSIM.library import defaults as d
from aMGSIM.library import functions as f
import json
import pandas as pd
import pyfaidx
from pathlib import Path
import subprocess

# import tqdm
import tqdm
from dataclasses import dataclass
from Bio import SeqIO
import itertools
import gzip
from mimetypes import guess_type
from aMGSIM import __version__

debug = None


def exceptionHandler(
    exception_type, exception, traceback, debug_hook=sys.__excepthook__
):
    """Print user friendly error messages normally, full traceback if DEBUG on.
    Adapted from http://stackoverflow.com/questions/27674602/hide-traceback-unless-a-debug-flag-is-set
    """
    if debug:
        print("\n*** Error:")
        # raise
        debug_hook(exception_type, exception, traceback)
    else:
        print("{}: {}".format(exception_type.__name__, exception))


logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)


# logging
logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)


def _rename_name(x, output_dir):
    taxon = x[0]
    file = x[1]

    out_dir = os.path.join(output_dir, "genomes")
    out_file = os.path.join(out_dir, "{}.fasta".format(taxon))

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # regexes
    re0 = re.compile(r".+ complete genome. ")
    re1 = re.compile(r"\W")
    re2 = re.compile(r"^_*(.*?)_*$")
    re3 = re.compile(r"_*complete_genome")
    re4 = re.compile(r"(.{78}).+")
    ambig_chars = re.compile(r"[RYSWKMBVDH]")
    taxon = re0.sub("", taxon)
    taxon = re1.sub("_", taxon)
    taxon = re2.sub(r"\1", taxon)
    taxon = re3.sub("", taxon)
    taxon = re4.sub(r"\1", taxon)
    # iterating through sequence
    ambig_cnt = 0
    encoding = guess_type(file)[1]  # uses file extension
    _open = partial(gzip.open, mode="rt") if encoding == "gzip" else open

    with open(out_file, mode="w") as outFH:
        with _open(file) as f:
            for i, record in enumerate(SeqIO.parse(f, "fasta")):
                name = record.name
                name = re0.sub("", name)
                name = re1.sub("_", name)
                name = re2.sub(r"\1", name)
                name = re3.sub("", name)
                name = re4.sub(r"\1", name)
                name = "{}__seq-{}".format(taxon, i)
                # name = name.lstrip('>') + '__seq{}'.format(i)
                record.id = name
                record.description = name
                ambig_cnt += len(ambig_chars.findall(str(record.seq.upper())))
                SeqIO.write(record, outFH, "fasta")
    return x + [out_file]


def rename_genomes(genome_table, cpus, output_dir):

    func = partial(_rename_name, output_dir=output_dir)
    genomes = genome_table[["Taxon", "Fasta", "Genome_size"]].values.tolist()

    if debug is True:
        files = list(map(func, genomes))
    else:
        p = Pool(cpus)
        files = list(tqdm.tqdm(p.imap_unordered(func, genomes), total=len(genomes)))
    df = pd.DataFrame.from_records(
        files, columns=["Taxon", "Fasta", "Genome_size", "Fasta_normalized"]
    )
    return df


def obj_dict(obj):
    return obj.__dict__


def run_fragSim(exe, params, seq_depth, ofile, tmp_dir, frag_ofile, fasta, debug):

    parms = process_params(params)

    cmd = "{exe} {params} -n {seq_depth} -o {ofile} -tmp {tmp_dir} -f {frag_ofile} {fasta}"
    cmd = cmd.format(
        exe=exe,
        params=parms,
        seq_depth=seq_depth,
        ofile=ofile,
        fasta=fasta,
        tmp_dir=tmp_dir,
        frag_ofile=frag_ofile,
    )

    # system call
    if debug is True:
        sys.stderr.write("CMD: " + cmd + "\n")
    try:
        res = subprocess.run(
            cmd, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    except subprocess.CalledProcessError as e:
        raise e
    if debug is True:
        sys.stderr.write(res.stderr.decode() + "\n")
        sys.stderr.write(res.stdout.decode() + "\n")


def process_params(params):
    parms = {}
    for k in params:
        if params[k] is not False:
            if params[k] is True:
                parms[k] = ""
            elif type(params[k]) is list:
                parms[k] = ",".join(map(str, params[k]))
            else:
                parms[k] = params[k]

    parms = " ".join(["{} {}".format(k, v) for k, v in parms.items()])
    return parms


def run_art(exe, params, seqSys, fasta, ofile, read_len, library, debug):
    parms = process_params(params)

    if library == "pe":
        if params["--qprof1"] and params["--qprof2"]:
            cmd = "{exe} {params} -p -c 1 -l {read_len} -amp -na -o {ofile} -i {fasta}"
        else:
            cmd = "{exe} {params} -p -c 1 -ss {seqSys} -l {read_len} -amp -na -o {ofile} -i {fasta}"
    elif library == "se":
        if params["--qprof1"] and params["--qprof2"]:
            cmd = "{exe} {params} -c 1 -l {read_len} -amp -na -o {ofile} -i {fasta}"

        else:
            cmd = "{exe} {params} -c 1 -ss {seqSys} -l {read_len} -amp -na -o {ofile} -i {fasta}"

    cmd = cmd.format(
        exe=exe,
        params=parms,
        ofile=ofile,
        fasta=fasta,
        read_len=read_len,
        seqSys=seqSys,
    )

    # system call
    if debug is True:
        sys.stderr.write("CMD: " + cmd + "\n")
    try:
        res = subprocess.run(
            cmd, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    except subprocess.CalledProcessError as e:
        raise e
    if debug is True:
        sys.stderr.write(res.stderr.decode() + "\n")
        sys.stderr.write(res.stdout.decode() + "\n")


def run_deamSim(exe, params, ofile, fasta, libprep, debug):

    parms = process_params(params)

    if libprep == "double" and params["-mapdamage"]:
        cmd = "{exe} {params} double -o {ofile} -name {fasta}"
    elif libprep == "single" and params["-mapdamage"]:
        cmd = "{exe} {params} single -o {ofile} -name {fasta}"
    else:
        cmd = "{exe} {params} -o {ofile} -name {fasta}"

    cmd = cmd.format(exe=exe, params=parms, ofile=ofile, fasta=fasta)

    # system call
    if debug is True:
        sys.stderr.write("CMD: " + cmd + "\n")
    try:
        res = subprocess.run(
            cmd, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    except subprocess.CalledProcessError as e:
        raise e
    if debug is True:
        sys.stderr.write(res.stderr.decode() + "\n")
        sys.stderr.write(res.stdout.decode() + "\n")


def run_adptSim(exe, params, ofile, fasta, library, read_len, debug):

    parms = process_params(params)

    if library == "pe":
        cmd = "{exe} {params} -artp {ofile} -l {read_len} {fasta}"
    else:
        cmd = "{exe} {params} -arts {ofile} -l {read_len} {fasta}"

    cmd = cmd.format(exe=exe, params=parms, ofile=ofile, read_len=read_len, fasta=fasta)

    # system call
    if debug is True:
        sys.stderr.write("CMD: " + cmd + "\n")
    try:
        res = subprocess.run(
            cmd, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    except subprocess.CalledProcessError as e:
        raise e
    if debug is True:
        sys.stderr.write(res.stderr.decode() + "\n")
        sys.stderr.write(res.stdout.decode() + "\n")


def prepare_data_fragments(
    x,
    frag_type,
    tmp_dir,
):
    fa_suffix = ".fasta"
    fa_suffix_gz = ".fasta.gz"

    fragSim_frag_suffix = ".tsv"

    entry = "fragments_{}".format(frag_type)

    fragment_data = {}

    fragSim_fname = "{}---{}___fragSim-{}".format(x["comm"], x["taxon"], frag_type)

    fragments = pd.DataFrame(x[entry]["fragments"])

    fragment_data["seq_depth"] = x[entry]["seq_depth"]

    # Save fragmets
    fragment_data["fragSim_frag_ofile"] = Path(tmp_dir, fragSim_fname).with_suffix(
        fragSim_frag_suffix
    )

    fragments.to_csv(
        path_or_buf=fragment_data["fragSim_frag_ofile"],
        sep="\t",
        header=False,
        index=False,
    )

    fragment_data["fragSim_ofile"] = Path(tmp_dir, fragSim_fname).with_suffix(
        fa_suffix_gz
    )

    deamSim_fname = "{}---{}___deamSim-{}".format(x["comm"], x["taxon"], frag_type)

    fragment_data["deamSim_ofile"] = Path(tmp_dir, deamSim_fname).with_suffix(
        fa_suffix_gz
    )

    adptSim_fname = "{}---{}___adptSim-{}".format(x["comm"], x["taxon"], frag_type)

    fragment_data["adptSim_ofile"] = Path(tmp_dir, adptSim_fname).with_suffix(fa_suffix)

    art_fname = "{}---{}___art-{}.".format(x["comm"], x["taxon"], frag_type)

    fragment_data["art_ofile"] = Path(tmp_dir, art_fname)

    return fragment_data


@dataclass
class Returnvalue:
    ancient: dict()
    modern: dict()


def collect_file_names(
    x, art_ofile, fragSim_ofile, deamSim_ofile, adptSim_ofile, frag_type, library
):
    if library == "pe":
        r1 = "{}1.fq".format(art_ofile)
        r2 = "{}2.fq".format(art_ofile)
        sr = None
    else:
        r1 = None
        r2 = None
        sr = "{}fq".format(art_ofile)

        # Collect file names
    files = {
        "comm": x["comm"],
        "taxon": x["taxon"],
        "frag_type": frag_type,
        "fragSim_ofile": fragSim_ofile,
        "deamSim_ofile": deamSim_ofile,
        "adptSim_ofile": adptSim_ofile,
        "art_ofile_r1": r1,
        "art_ofile_r2": r2,
        "art_ofile_sr": sr,
    }
    return files


def generate_fragments(
    x,
    fragSim_exe,
    fragSim_params,
    deamSim_exe,
    deamSim_params,
    adptSim_exe,
    adptSim_params,
    art_exe,
    art_params,
    libprep,
    tmp_dir,
    exp_data,
    debug,
    genome_table,
):
    genome_data = genome_table[(genome_table["Taxon"] == x["taxon"])]
    fasta = genome_data["Fasta_normalized"].item()
    # Create index
    fasta_seq = pyfaidx.Faidx(fasta)
    read_len = exp_data["read_length"]
    library = exp_data["library"]
    seqSys = exp_data["seqSys"]
    seq_depth = exp_data["n_reads"]
    files_modern = {}
    files_ancient = {}
    # Case when onlyAncient is False
    # Here we will need to run for modern and ancient
    if x["onlyAncient"] is False:
        if x["fragments_ancient"] is not None:
            # Run fragSim
            frag_type = "ancient"
            frag_data = prepare_data_fragments(
                x=x, frag_type=frag_type, tmp_dir=tmp_dir
            )
            run_fragSim(
                exe=fragSim_exe,
                params=fragSim_params["ancient"],
                seq_depth=frag_data["seq_depth"],
                ofile=frag_data["fragSim_ofile"],
                tmp_dir=tmp_dir,
                frag_ofile=frag_data["fragSim_frag_ofile"],
                fasta=fasta,
                debug=debug,
            )
            # Run deamSim
            run_deamSim(
                exe=deamSim_exe,
                params=deamSim_params,
                ofile=frag_data["deamSim_ofile"],
                fasta=frag_data["fragSim_ofile"],
                libprep=libprep,
                debug=debug,
            )
            # Regex to parse deamSim headers: (\S+):([+-]):(\d+):(\d+):(\d+)(?:_DEAM:(.*))?
            run_adptSim(
                exe=adptSim_exe,
                params=adptSim_params,
                ofile=frag_data["adptSim_ofile"],
                fasta=frag_data["deamSim_ofile"],
                read_len=read_len,
                library=library,
                debug=debug,
            )

            # rename_sequences(fasta=frag_data["adptSim_ofile"])

            run_art(
                exe=art_exe,
                params=art_params,
                seqSys=seqSys,
                fasta=frag_data["adptSim_ofile"],
                ofile=frag_data["art_ofile"],
                read_len=read_len,
                library=library,
                debug=debug,
            )

            files_ancient = collect_file_names(
                fragSim_ofile=frag_data["fragSim_ofile"],
                deamSim_ofile=frag_data["deamSim_ofile"],
                adptSim_ofile=frag_data["adptSim_ofile"],
                library=library,
                art_ofile=frag_data["art_ofile"],
                frag_type=frag_type,
                x=x,
            )

        if x["fragments_modern"] is not None:
            frag_type = "modern"
            # Run fragSim
            frag_data = prepare_data_fragments(
                x=x, frag_type=frag_type, tmp_dir=tmp_dir
            )
            run_fragSim(
                exe=fragSim_exe,
                params=fragSim_params["modern"],
                seq_depth=frag_data["seq_depth"],
                ofile=frag_data["fragSim_ofile"],
                tmp_dir=tmp_dir,
                frag_ofile=frag_data["fragSim_frag_ofile"],
                fasta=fasta,
                debug=debug,
            )
            run_adptSim(
                exe=adptSim_exe,
                params=adptSim_params,
                ofile=frag_data["adptSim_ofile"],
                fasta=frag_data["fragSim_ofile"],
                read_len=read_len,
                library=library,
                debug=debug,
            )
            run_art(
                exe=art_exe,
                params=art_params,
                seqSys=seqSys,
                fasta=frag_data["adptSim_ofile"],
                ofile=frag_data["art_ofile"],
                read_len=read_len,
                library=library,
                debug=debug,
            )
            files_modern = collect_file_names(
                fragSim_ofile=frag_data["fragSim_ofile"],
                deamSim_ofile=frag_data["deamSim_ofile"],
                adptSim_ofile=frag_data["adptSim_ofile"],
                library=library,
                art_ofile=frag_data["art_ofile"],
                frag_type=frag_type,
                x=x,
            )

    if x["onlyAncient"] is True:
        frag_type = "ancient"
        # Run fragSim
        frag_data = prepare_data_fragments(x=x, frag_type=frag_type, tmp_dir=tmp_dir)
        run_fragSim(
            exe=fragSim_exe,
            params=fragSim_params["ancient"],
            seq_depth=frag_data["seq_depth"],
            ofile=frag_data["fragSim_ofile"],
            tmp_dir=tmp_dir,
            frag_ofile=frag_data["fragSim_frag_ofile"],
            fasta=fasta,
            debug=debug,
        )
        # Run deamSim
        run_deamSim(
            exe=deamSim_exe,
            params=deamSim_params,
            ofile=frag_data["deamSim_ofile"],
            fasta=frag_data["fragSim_ofile"],
            libprep=libprep,
            debug=debug,
        )
        run_adptSim(
            exe=adptSim_exe,
            params=adptSim_params,
            ofile=frag_data["adptSim_ofile"],
            fasta=frag_data["deamSim_ofile"],
            read_len=read_len,
            library=library,
            debug=debug,
        )
        run_art(
            exe=art_exe,
            params=art_params,
            seqSys=seqSys,
            fasta=frag_data["adptSim_ofile"],
            ofile=frag_data["art_ofile"],
            read_len=read_len,
            library=library,
            debug=debug,
        )
        files_ancient = collect_file_names(
            fragSim_ofile=frag_data["fragSim_ofile"],
            deamSim_ofile=frag_data["deamSim_ofile"],
            adptSim_ofile=frag_data["adptSim_ofile"],
            library=library,
            art_ofile=frag_data["art_ofile"],
            frag_type=frag_type,
            x=x,
        )
        # Regex to parse deamSim headers: (\S+):([+-]):(\d+):(\d+):(\d+)(?:_DEAM:(.*))?
    if x["onlyAncient"] is None:
        frag_type = "modern"
        # Run fragSim
        frag_data = prepare_data_fragments(x=x, frag_type=frag_type, tmp_dir=tmp_dir)
        run_fragSim(
            exe=fragSim_exe,
            params=fragSim_params["modern"],
            seq_depth=frag_data["seq_depth"],
            ofile=frag_data["fragSim_ofile"],
            tmp_dir=tmp_dir,
            frag_ofile=frag_data["fragSim_frag_ofile"],
            fasta=fasta,
            debug=debug,
        )
        run_adptSim(
            exe=adptSim_exe,
            params=adptSim_params,
            ofile=frag_data["adptSim_ofile"],
            fasta=frag_data["fragSim_ofile"],
            read_len=read_len,
            library=library,
            debug=debug,
        )
        run_art(
            exe=art_exe,
            params=art_params,
            seqSys=seqSys,
            fasta=frag_data["adptSim_ofile"],
            ofile=frag_data["art_ofile"],
            read_len=read_len,
            library=library,
            debug=debug,
        )
        files_modern = collect_file_names(
            fragSim_ofile=frag_data["fragSim_ofile"],
            deamSim_ofile=frag_data["deamSim_ofile"],
            adptSim_ofile=frag_data["adptSim_ofile"],
            library=library,
            art_ofile=frag_data["art_ofile"],
            frag_type=frag_type,
            x=x,
        )

    os.remove(frag_data["fragSim_frag_ofile"])
    return Returnvalue(files_ancient, files_modern)


def get_ancient_genomes_data(self):
    self = {
        k: self.get(k, None)
        for k in ("comm", "taxon", "onlyAncient", "library", "seqSys", "coverage")
    }
    return self


def validations(self):
    mapdamage = self["deamSim"]["-mapdamage"]
    damage = self["deamSim"]["-damage"]

    if (mapdamage is not False) and (damage is not False):
        raise IOError("deamSim: You cannot specify both mapdamage and damage")
    elif (mapdamage is False) and (damage is False):
        raise IOError("deamSim: Both mapdamage and damage cannot be False")
    else:
        pass


def rename_reads(records, pattern, comm, taxon, frag_type, fastx):
    for i, record in records:
        # renaming fastq read
        m1 = re.match(pattern, record.id)
        # Read name as:
        # 1: Sample
        # 2: genome or genome_contig
        # 3: read number
        # 4: type (ancient/modern)
        # 5: strand
        # 6: start
        # 7: end
        # 8: length
        # 9: damage positions in read
        if fastx == "fastq":
            name = "{}___{}---{}:{}:{}:{}:{}:{}:{}/{}".format(
                comm,
                m1.group(1),
                i,
                frag_type,
                m1.group(2),
                m1.group(3),
                m1.group(4),
                m1.group(5),
                m1.group(6),
                m1.group(7),
            )
        elif fastx == "fasta":
            name = "{}___{}---{}:{}:{}:{}:{}:{}:{}".format(
                comm,
                m1.group(1),
                i,
                frag_type,
                m1.group(2),
                m1.group(3),
                m1.group(4),
                m1.group(5),
                m1.group(6),
            )
        record.id = name
        record.description = name
        # record.description = 'type:{} strand:{} start:{} stop:{} len:{} deam:{}'.format(
        #     frag_type, m1.group(2), m1.group(3), m1.group(4), m1.group(5), m1.group(6), m1.group(7))
        yield record


def _combine_reads(read_files, output_dir, output_file, fastx):
    """Combine fastq read files into 1 read file.
    Parameters
    ----------
    read_files : list
        All read files to combine
    output_dir : str
        Output directory path
    output_file : str
        Output file path
    """
    output_file = os.path.join(output_dir, output_file)

    if fastx == "fastq":
        p0 = re.compile("(\S+)---(\S+)___art-(\S+).\d.fq")
        p1 = re.compile("(\S+):([+\-]):(\d+):(\d+):(\d+)(?:_DEAM:(.*))?\-\d\/(\d)$")
    elif fastx == "fasta":
        p0 = re.compile("(\S+)---(\S+)___\S+-(\S+).fasta.gz")
        p1 = re.compile("(\S+):([+\-]):(\d+):(\d+):(\d+)(?:_DEAM:(.*))?")

    with open(output_file, "w") as outFH:
        for in_file in read_files:
            m0 = re.match(p0, os.path.split(in_file)[1])
            comm = m0.group(1)
            taxon = m0.group(2)
            frag_type = m0.group(3)
            encoding = guess_type(in_file)[1]  # uses file extension
            _open = partial(gzip.open, mode="rt") if encoding == "gzip" else open
            with _open(in_file) as f:
                records = enumerate(SeqIO.parse(f, fastx))
                SeqIO.write(
                    rename_reads(
                        records,
                        pattern=p1,
                        comm=comm,
                        taxon=taxon,
                        frag_type=frag_type,
                        fastx=fastx,
                    ),
                    outFH,
                    fastx,
                )
            # delete temporary file
            # os.remove(in_file)
    return output_file


# TODO: Add single reads option


def _combine_fastq_types(
    file_type, data_modern, data_ancient, output_dir, comm, suffix
):

    output_dir = os.path.join(output_dir, str(comm))
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    if "r1" in file_type:
        r1 = list(
            set(
                [s[file_type] for s in data_modern]
                + [s[file_type] for s in data_ancient]
            )
        )
        r1 = sorted(r1)
        R1_files = _combine_reads(
            read_files=r1,
            output_dir=output_dir,
            output_file="comm-{}_{}.1.fq".format(comm, suffix),
            fastx="fastq",
        )
        files = {
            "comm": comm,
            "file_type": file_type.replace("_r1", ""),
            "pair": "R1",
            "file": R1_files,
        }

    elif "r2" in file_type:
        r2 = list(
            set(
                [s[file_type] for s in data_modern]
                + [s[file_type] for s in data_ancient]
            )
        )
        r2 = sorted(r2)
        R2_files = _combine_reads(
            read_files=r2,
            output_dir=output_dir,
            output_file="comm-{}_{}.2.fq".format(comm, suffix),
            fastx="fastq",
        )
        files = {
            "comm": comm,
            "file_type": file_type.replace("_r2", ""),
            "pair": "R2",
            "file": R2_files,
        }
    else:
        sr = list(
            set(
                [s[file_type] for s in data_modern]
                + [s[file_type] for s in data_ancient]
            )
        )

    return files


# TODO: CHeck if list of files are empty


def _combine_fasta_types(
    file_type, data_modern, data_ancient, output_dir, comm, suffix
):
    if suffix == "deamSim":
        sr = list([s[file_type] for s in data_ancient])
    else:
        sr = list(
            set(
                [s[file_type] for s in data_modern]
                + [s[file_type] for s in data_ancient]
            )
        )
    # if not sr:
    sr = sorted(sr)
    output_dir = os.path.join(output_dir, str(comm))

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    # read1
    SR_files = _combine_reads(
        read_files=sr,
        output_dir=output_dir,
        output_file="comm-{}_{}.fa".format(comm, suffix),
        fastx="fasta",
    )
    files = {"comm": comm, "file_type": file_type, "pair": "SR", "file": SR_files}
    return files


def combine_fastx_files(x, output_dir, ancient_files, modern_files, exp_data):
    comm = x[0]
    file_type = x[1]
    data_modern = list(filter(lambda d: d["comm"] == comm, modern_files))
    data_ancient = list(filter(lambda d: d["comm"] == comm, ancient_files))
    files = {
        "art_ofile_r1": "art",
        "art_ofile_r2": "art",
        "fragSim_ofile": "fragSim",
        "deamSim_ofile": "deamSim",
    }
    if any(x in file_type for x in ["r1", "r2"]):
        files = _combine_fastq_types(
            file_type=file_type,
            suffix=files[file_type],
            data_modern=data_modern,
            data_ancient=data_ancient,
            output_dir=output_dir,
            comm=comm,
        )
    else:
        files = _combine_fasta_types(
            file_type=file_type,
            suffix=files[file_type],
            data_modern=data_modern,
            data_ancient=data_ancient,
            output_dir=output_dir,
            comm=comm,
        )
    return files


def main(args):

    global debug
    if args["--debug"]:
        debug = True
    else:
        debug = None

    sys.excepthook = exceptionHandler

    # simulating reads
    args = f.validate_schema(args, d.schema_init_ar, debug)

    config = f.get_config(
        config=args["<config>"], schema=d.ar_schema_config, debug=debug
    )

    # art_params = config_params['art_config']
    config_params = config["global"]

    val = validations(config)

    fragSim_params = config["fragSim"]
    deamSim_params = config["deamSim"]
    adptSim_params = config["adptSim"]
    art_params = config["art"]
    adptRem_params = config["AdapterRemoval"]

    libprep = config_params["libprep"]

    # Create folders
    tmp_dir = config_params["tmp_dir"]
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir, exist_ok=True)
    else:
        logging.info("Folder {} already exists. Deleting...".format(tmp_dir))
        rmtree(tmp_dir)
        os.makedirs(tmp_dir, exist_ok=True)

    output_dir = config_params["output_dir"]
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    else:
        logging.info("Folder {} already exists. Deleting...".format(output_dir))
        rmtree(output_dir)
        os.makedirs(output_dir, exist_ok=True)

    logging.info("Loading genomes...")
    # load tables
    genome_table = f.load_genome_table(args["<genome_table>"])

    logging.info("Normalizing genome names...")
    genome_table = rename_genomes(
        genome_table=genome_table, output_dir=output_dir, cpus=config_params["cpus"]
    )

    logging.info("Simulating reads by sample...")
    with open(fragSim_params["ancient_genomes"], "r") as json_file:
        filename = json_file
        ancient_genomes = json.load(json_file)
    ancient_genomes_data = ancient_genomes["data"]
    ancient_genomes_exp = ancient_genomes["experiment"]
    # fragSim_params.pop('ancient_genomes', None)

    func = partial(
        generate_fragments,
        fragSim_exe=config_params["fragSim"],
        fragSim_params=fragSim_params,
        deamSim_exe=config_params["deamSim"],
        deamSim_params=deamSim_params,
        adptSim_exe=config_params["adptSim"],
        adptSim_params=adptSim_params,
        art_exe=config_params["art"],
        art_params=art_params,
        libprep=libprep,
        tmp_dir=tmp_dir,
        debug=debug,
        genome_table=genome_table,
        exp_data=ancient_genomes_exp,
    )

    if debug is True:
        files = list(map(func, ancient_genomes_data))
    else:
        p = Pool(config_params["cpus"])
        files = list(
            tqdm.tqdm(
                p.imap_unordered(func, ancient_genomes_data),
                total=len(ancient_genomes_data),
            )
        )
        # fragments_data = process_map(func, ancient_genomes_data, max_workers=config_params['cpus'])
    ancient_files = list(filter(None, map(lambda x: x.ancient, files)))
    modern_files = list(filter(None, map(lambda x: x.modern, files)))
    comms = list(
        set([s["comm"] for s in ancient_files] + [s["comm"] for s in modern_files])
    )

    func = partial(
        combine_fastx_files,
        output_dir=output_dir,
        ancient_files=ancient_files,
        modern_files=modern_files,
        exp_data=ancient_genomes_exp,
    )

    logging.info("Combining reads by sample...")

    files = ["art_ofile_r1", "art_ofile_r2", "fragSim_ofile", "deamSim_ofile"]

    comm_files = list(itertools.product(comms, files))
    if debug is True:
        files = list(map(func, comm_files))
    else:
        p = Pool(config_params["cpus"])
        files = list(
            tqdm.tqdm(p.imap_unordered(func, comm_files), total=len(comm_files))
        )

    df = pd.DataFrame(files).sort_values("comm")

    read_files = {}

    df_files = df.groupby(["comm", "file_type"])[["pair", "file"]].apply(
        lambda x: dict(x.values)
    )

    for comm in df["comm"].unique():
        comm = str(comm)
        read_files[comm] = {}
        for index, val in df_files.iteritems():
            if comm == str(index[0]):
                read_files[comm][index[1]] = val

    read_files = {
        "reads": read_files,
        "genomes": genome_table["Fasta_normalized"].values.tolist(),
    }
    file_name = Path(filename.name).stem
    dir_name = Path(filename.name).parent

    suffix = ".json"
    read_files_ofile = Path(dir_name, "{}_read-files".format(file_name)).with_suffix(
        suffix
    )

    read_files_json = json.dumps(
        read_files, default=obj_dict, ensure_ascii=False, indent=4
    )

    with open(read_files_ofile, "w", encoding="utf-8") as outfile:
        print(read_files_json, file=outfile)

    logging.info("Ancient synthetic reads generated.")


def opt_parse(args=None):
    version = "Version: " + __version__
    if args is None:
        args = docopt(__doc__, version=version)
    else:
        args = docopt(__doc__, version=version, argv=args)
    main(args)


# %%
