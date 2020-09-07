from __future__ import print_function
# import
# batteries
import os
import sys
import re
import logging
import subprocess
from functools import partial
from multiprocessing import Pool

# 3rd party
import pandas as pd
from Bio import SeqIO

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


# functions
def tidy_taxon_names(x):
    """Remove special characters from taxon names
    """
    x = re.sub(r'[()\/:;, ]+', '_', x)
    return x


def load_genome_table(in_file):
    """Loading genome table
    Parameters
    ----------
    in_file : str
        input file path
    """
    df = pd.read_csv(in_file, sep='\t')

    # check headers
    diff = set(['Taxon', 'Fasta']) - set(df.columns.values)
    if len(diff) > 0:
        diff = ','.join(diff)
        raise ValueError('Cannot find table columns: {}'.format(diff))

    # getting genome sizes
    df['Genome_size'] = [_genome_size(x) for i, x in df.iterrows()]

    # tidy taxon names
    df['Taxon'] = df['Taxon'].astype(str).apply(tidy_taxon_names)

    return df

# Load ancient genomes parameters
def load_ancient_genomes(json_file):
    pass


def _genome_size(x):
    """Get the total bp for the genome sequence of all genomes

    Parameters
    ----------
    x : pd.Series
       Series that includes 'Fasta' in the index
    """
    bp = 0
    for record in SeqIO.parse(x['Fasta'], 'fasta'):
        bp += len(record.seq)
    return bp


def load_abund_table(in_file):
    """Loading abundance table
    Parameters
    ----------
    in_file : str
        input file path
    """
    df = pd.read_csv(in_file, sep='\t')

    # check headers
    diff = set(['Community', 'Taxon', 'Perc_rel_abund']) - \
        set(df.columns.values)
    if len(diff) > 0:
        diff = ','.join(diff)
        raise ValueError('Cannot find table columns: {}'.format(diff))

    # tidy taxon names
    df['Taxon'] = df['Taxon'].astype(str).apply(tidy_taxon_names)

    return df


def create_fragSim_params(fragSim_params, type):
    # Create params for ancient and modern
    if type == "ancient":

        fragSim_params = ' '.join(['--{} {}'.format(k, v)
                                   for k, v in fragSim_params.items()])
    else:
        pass


def sample_taxon_list(genome_table, abund_table):
    """Creating [sample,taxon] list of lists
    Parameters
    ----------
    genome_table : pd.DataFrame
    abund_table : pd.DataFrame
    """
    # joining tables
    df = abund_table.merge(genome_table, on=['Taxon'])
    
    # convert to a list of lists
    sample_taxon = []
    cols = ['Community', 'Taxon', 'Genome_size',
            'Fasta', 'Perc_rel_abund']
    for i, x in df[cols].iterrows():
        sample_taxon.append(x.tolist())
    return sample_taxon

#def generate_random_ancient_taxon(self):
    
    


def sim_illumina_ancient(sample_taxon, output_dir, seq_depth,
                         fragSim_exe, fragSim_params,
                         art_params, art_exe,
                         temp_dir, nproc=1, debug=False):
    """Simulate illumina reads
    Parameters
    ----------
    sample_taxon : list
        [Community,Taxon,Genome_size,Fasta,Perc_rel_abund,Fold]
    output_dir : str
        Output director for all read files
    seq_depth : int
        Sequencing depth per sample
    fragSim_exe: str
        fragSim binary
    fragSim_params: dict
        Parameters provided to fragSim
    art_exe: str
        ART binary
    art_params : dict
        Parameters provided to art_illumina
    temp_dir : str
        Temporary file directory
    nproc : int
        Number of parallel processes
    debug : bool
        Debug mode
    """

    # calculating fold per genome
    for x in sample_taxon:
        genome_size = float(x[2])
        perc_rel_abund = float(x[4])
        try:
            _ = art_params['paired']
            paired = 2
        except KeyError:
            try:
                art_params['mflen']
                paired = 2
            except KeyError:
                paired = 1
        read_length = art_params['len'] * paired
        seq_depth_ancient = int((cov_ancient * genome_size)/read_length)
        seq_depth_modern = int(seq_depth - seq_depth_ancient)
        fold_modern = (perc_rel_abund/100) * \
            ((seq_depth_modern * read_length) / genome_size)
        fold_ancient = (perc_rel_abund/100) * \
            ((seq_depth_ancient * read_length) / genome_size)
        print('genome:{} fold_modern: {} fold_ancient:{} perc_rel_abund:{} seq_depth_modern:{} seq_depth_ancient: {} read_length:{} genome_size:{}'.format(
            x[1], fold_modern, fold_ancient, perc_rel_abund, seq_depth_modern, seq_depth_ancient, read_length, genome_size))
        x.append(seq_depth_modern)
        x.append(fold_modern)
        x.append(seq_depth_ancient)
        x.append(fold_ancient)

    # directories
    # temp dir
    if not os.path.isdir(temp_dir):
        os.makedirs(temp_dir)
    # output dir
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # simulate per sample
    logging.info('Simulating reads...')
    func = partial(sim_art_ancient,
                   art_exe=art_exe,
                   art_params=art_params,
                   temp_dir=temp_dir,
                   debug=debug)
    if debug is True:
        fq_files = map(func, sample_taxon)
    else:
        p = Pool(nproc)
        fq_files = p.map(func, sample_taxon)
    fq_files = list(fq_files)

    # combining all reads by sample
    logging.info('Combining simulated reads by sample...')
    comms = list(set([x[0] for x in sample_taxon]))
    func = partial(combine_reads_by_sample,
                   fq_files=fq_files,
                   temp_dir=temp_dir,
                   file_prefix='illumina',
                   output_dir=output_dir,
                   debug=debug)
    if debug is True:
        res = map(func, comms)
    else:
        p = Pool(nproc)
        res = p.map(func, comms)
    res = list(res)

    # removing temp dir
    logging.info('Removing temp directory...')
    # rmtree(temp_dir)

    # status
    for sample_list in res:
        for file_name in sample_list:
            if file_name is not None:
                logging.info('File written: {}'.format(file_name))


def sim_art_ancient(x, art_exe, art_params, temp_dir, debug=False):
    """Simulate illumina reads
    Parameters
    ----------
    x : list
        [Community,Taxon,Genome_size,Fasta,Perc_rel_abund,Fold]
    """
    community = str(x[0])
    taxon = str(x[1])
    fasta = str(x[3])
    perc_rel_abund = float(x[5])
    fold_modern = float(x[6])

    # output
    # temporary directories
    out_dir = os.path.join(temp_dir, str(community), taxon)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    # temporary files
    output_prefix = os.path.join(out_dir, 'illumina')

    # 1. Create fragments with fragSim

    # Modern sequences

    # art command
    # cmd = '{fragSim_exe} {fragSim_params} --noALN -f {fold_modern} -i {input} -o {output_prefix}'
    # cmd = cmd.format(art_exe=art_exe,
    #                  art_params=art_params,
    #                  fold_modern=fold_modern,
    #                  input=fasta,
    #                  output_prefix=output_prefix)
    # print(cmd)
    # system call
    if debug is True:
        sys.stderr.write('CMD: ' + cmd + '\n')
    try:
        res = subprocess.run(cmd, check=True, shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        raise e
    if debug is True:
        sys.stderr.write(res.stderr.decode() + '\n')
        sys.stderr.write(res.stdout.decode() + '\n')

    # check that files have been created
    R0_file = output_prefix + '.fq'
    R1_file = output_prefix + '1.fq'
    R2_file = output_prefix + '2.fq'
    if os.path.isfile(R1_file) and os.path.isfile(R2_file):
        return [community, R1_file, R2_file]
    elif os.path.isfile(R0_file):
        return [community, R0_file]
    else:
        msg = 'Cannot find art_illumina output files!'
        raise ValueError(msg)

    # 2. Add deamination to the ancient sequences

    # 3. Add adapters

    # 4. Add errors with ART

    art_params = ' '.join(['--{} {}'.format(k, v)
                           for k, v in art_params.items()])
    # art command
    cmd = '{art_exe} {art_params} --noALN -f {fold_modern} -i {input} -o {output_prefix}'
    cmd = cmd.format(art_exe=art_exe,
                     art_params=art_params,
                     fold_modern=fold_modern,
                     input=fasta,
                     output_prefix=output_prefix)
    print(cmd)
    # system call
    if debug is True:
        sys.stderr.write('CMD: ' + cmd + '\n')
    try:
        res = subprocess.run(cmd, check=True, shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        raise e
    if debug is True:
        sys.stderr.write(res.stderr.decode() + '\n')
        sys.stderr.write(res.stdout.decode() + '\n')

    # check that files have been created
    R0_file = output_prefix + '.fq'
    R1_file = output_prefix + '1.fq'
    R2_file = output_prefix + '2.fq'
    if os.path.isfile(R1_file) and os.path.isfile(R2_file):
        return [community, R1_file, R2_file]
    elif os.path.isfile(R0_file):
        return [community, R0_file]
    else:
        msg = 'Cannot find art_illumina output files!'
        raise ValueError(msg)


def combine_reads_by_sample(sample, fq_files, temp_dir, file_prefix, output_dir, debug=False):
    """ Concat all sample-taxon read files into per-sample read files
    Parameters
    ----------
    sample : str
        Sample ID
    temp_dir : str
        Temporary directory path
    file_prefix : str
        Output file prefix
    output_dir : str
        Output directory path
    debug : bool
        Debug mode
    """
    # sorting fastq files by read pair
    sample = str(sample)
    R1_files = [x[1] for x in fq_files if x[0] == sample]
    try:
        R2_files = [x[2] for x in fq_files if x[0] == sample]
    except IndexError:
        R2_files = None

    # writing files
    output_dir = os.path.join(output_dir, str(sample))
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    # read1
    R1_files = _combine_reads(R1_files, output_dir, 'R1.fq')
    # read2
    if R2_files is not None:
        R2_files = _combine_reads(R2_files, output_dir, 'R2.fq')

    return [R1_files, R2_files]


def _combine_reads(read_files, output_dir, output_file):
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
    with open(output_file, 'w') as outFH:
        for in_file in read_files:
            taxon = os.path.split(os.path.split(in_file)[0])[1]
            for i, record in enumerate(SeqIO.parse(in_file, 'fastq')):
                # renaming fastq read
                name = '{}__SEQ{}'.format(taxon, i)
                record.id = name
                record.description = name
                SeqIO.write(record, outFH, 'fastq')
            # delete temporary file
            os.remove(in_file)

    return output_file
