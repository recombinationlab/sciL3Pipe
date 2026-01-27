
from datetime import datetime
from pathlib import Path
import yaml
from shutil import which
from collections import defaultdict
import json
import glob
import pandas as pd
import os, sys, pwd, re
import socket, platform
import subprocess
import logging
from snakemake_interface_executor_plugins.settings import ExecMode
from snakemake.logging import logger, logger_manager

v = datetime.now()
run_date = v.strftime('%Y.%m.%d.')


snakemake.utils.min_version("9.11.6")

# TODO: config enforced by schema, don't need try except blocks
################################################################################
# Parse config
################################################################################

# The final output is tabular, we might need to indent subsequent lines correctly.
indent = 24

INFO = ""

try:
    out_dir = config['output_dir']
except:
    out_dir = ''
    INFO += 'Defaulting to working directory as output directory' + '\n' + (' ' * indent)

try:
    aligner = config['aligner']
except:
    sys.exit('ERROR: Need to specify aligner')

if aligner == 'bwa':
    try:
        bwa_index = config['bwa_index'][config['assembly']]
    except:
        sys.exit('ERROR: BWA index not specified in config.yaml')
elif aligner == 'bowtie2':
    try:
        bowtie2_index = config['bowtie2_index'][config['assembly']]
    except:
        sys.exit('ERROR: Bowtie2 index not specified in config.yaml')

try:
    assembly = config['assembly']
except:
    sys.exit('ERROR: Need to specify assembly e.g. "hg38" or "mm10"')


if config['generate_summary']:
    try:
        summary_for = config['summary_for']
        if summary_for not in ['wprefix', 'woprefix']:
            sys.exit('ERROR: summary_for can only be wpredix or woprefix')
    except:
        sys.exit('ERROR: Need to specify summary_for in config')
else:
    summary_for = ''

################################################################################
# Load samples 
################################################################################

# Fastq files
try:
    FILES = json.load(open(config['fastqs']))
except:
    FILES = None

# Samples should start with a general name e.g. yi292, followed by
# SSS file containing barcodes
ALL_SAMPLES = []
# POST_SSS_SAMPLES = defaultdict(list)
for samp in config['samples']:
    ALL_SAMPLES.append(samp)

################################################################################
# Load samples from trimmed
################################################################################

if config['from_trimmed'] != None:
    try:
        TRM_FILES = json.load(open(config['from_trimmed']))
        TRM_SAMPLES = [*TRM_FILES]
        trimmed_start = True
        INFO += 'Starting from trimmed fastq files. \n' + (' ' * indent)
    except:
        TRM_FILES = None
        trimmed_start = False
else:
    TRM_FILES = None
    trimmed_start = False
    
################################################################################
# Load SSS BAMs
################################################################################

if config['split_bam']:
    try:
        SSS_FILES = json.load(open(config['sss_bams']))
        SSS_SAMPLES = [*SSS_FILES]
        sss_bam_start = True
        # INFO += 'Starting from SSS BAM files. \n' + (' ' * indent)
    except:
        SSS_FILES = None
        sss_bam_start = False
else:
    SSS_FILES = None
    sss_bam_start = False


if config['split_bam'] and SSS_FILES != None:
    FILES = None

INFO += f'Starting from sss BAM: {sss_bam_start} \n' + (' ' * indent)


if SSS_FILES == None:
    samples = ALL_SAMPLES
    sss_all = list()
else:
    samples = SSS_SAMPLES
    sss_all = list()
    for s in samples:
        sss_all.extend([*SSS_FILES.get(s, [])])
    # INFO += f'SSS split samples to process {SSS_FILES} \n'

INFO += f'Samples to process: {samples} \n' + (' ' * indent)

#################################################################################
# Hybrid split and filter BAM
#################################################################################

if config["hybrid_filter"] is not None:
    if config["natsort_chromosomes"]:
        bam_filter="-f " + config["hybrid_filter"]
    else:
        bam_filter="-f " + config["hybrid_filter"] + "--natsort_off"
else:
    if config["natsort_chromosomes"]:
        bam_filter=""
    else:
        bam_filter="--natsort_off"

if config["additional_filter"] is not None:
    if config["natsort_chromosomes"]:
        bam_filter_2="-f " + config["additional_filter"]
    else:
        bam_filter_2="-f " + config["additional_filter"] + "--natsort_off"
else:
    if config["natsort_chromosomes"]:
        bam_filter_2=""
    else:
        bam_filter_2="--natsort_off"

    
if config["run_split"]:
    try:
        prefix_1 = config["prefix_1"]
        prefix_2 = config["prefix_2"]
        process_only = config["process_only"]
    except:
        sys.exit("ERROR: prefix or process only parameters not set in config")
else:
    try:
        process_only = config["process_only"]
        prefix_1 = "1"
        prefix_2 = "2"
    except:
        sys.exit("ERROR: process only parameter not set in config. Without split this can be considered a postfix.")




#################################################################################
# Log setup
#################################################################################

logdir = Path(out_dir) / "logs"
logdir.mkdir(parents=True, exist_ok=True)

extra_logfile = logdir / (datetime.now().isoformat().replace(":", "") + ".log")

file_handler = logging.FileHandler(extra_logfile)
file_handler.setLevel(logging.DEBUG)  # capture everything
file_handler.setFormatter(logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s"))
logger_manager.logger.addHandler(file_handler)


#################################################################################
# Basic Configuration
#################################################################################

# We want to report the sciL3Pipe version for the user, for reproducibility.
scil3pipe_version = '0.2.0'  

# After changing to our new scheme, we can verify the scheme to fit our expextation.
snakemake.utils.validate(config, schema='../schemas/config.schema.yaml')


#################################################################################
# Pipeline User Output
#################################################################################
# Adapted from https://github.com/moiexpositoalonsolab/grenepipe/blob/master/workflow/rules/initialize.smk

# Get a nicely formatted username and hostname
username = pwd.getpwuid(os.getuid())[0]
hostname = socket.gethostname()
hostname = hostname + ('; ' + platform.node() if platform.node() != socket.gethostname() else '')

# Get some info on the platform and OS
pltfrm = platform.platform() + '\n' + (' ' * indent) + platform.version()
try:
    # Not available in all versions, so we need to catch this
    ld = platform.linux_distribution()
    if len(ld):
        pltfrm += '\n' + (' ' * indent) + ld
    del ld
except:
    pass
try:
    # Mac OS version comes back as a nested tuple?!
    # Need to merge the tuples...
    def merge_tuple(x, bases=(tuple, list)):
        for e in x:
            if type(e) in bases:
                for e in merge_tuple(e, bases):
                    yield e
            else:
                yield e

    mv = ' '.join(merge_tuple(platform.mac_ver()))
    if not mv.isspace():
        pltfrm += '\n' + (' ' * indent) + mv
    del mv, merge_tuple
except:
    pass

# Get the git commit hash of sciL3Pipe, if available.
try:
    process = subprocess.Popen(
        ['git', 'rev-parse', '--short', 'HEAD'], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    out, err = process.communicate()
    out = out.decode('ascii')
    scil3pipe_git_hash = out.strip()
    if scil3pipe_git_hash:
        scil3pipe_version += '-' + scil3pipe_git_hash
    del process, out, err, scil3pipe_git_hash
except:
    pass

# Get the conda version, if available.
try:
    process = subprocess.Popen(['conda', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    out = out.decode('ascii')
    conda_ver = out[out.startswith('conda') and len('conda') :].strip()
    del process, out, err
    if not conda_ver:
        conda_ver = 'n/a'
except:
    conda_ver = 'n/a'

# Same for mamba. This somehow can also give a differing conda version.
# Who knows what that means. I'm sick of conda. Just reporting the version here,
# and have someone else deal with it.
try:
    process = subprocess.Popen(['mamba', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    out = out.decode('ascii')
    mamba_ver = re.findall('mamba *(.*) *', out)[0]
    conda_ver_mamba = re.findall('conda *(.*) *', out)[0]
    del process, out, err
    if not mamba_ver:
        mamba_ver = 'n/a'
        conda_ver_mamba = ''
except:
    mamba_ver = 'n/a'
    conda_ver_mamba = ''
if conda_ver_mamba and conda_ver_mamba != conda_ver:
    conda_ver += ' (conda), ' + conda_ver_mamba + ' (mamba)'

# Get the conda env name, if available.
# See https://stackoverflow.com/a/42660674/4184258
conda_env = os.environ['CONDA_DEFAULT_ENV'] + ' (' + os.environ['CONDA_PREFIX'] + ')'
if conda_env == ' ()':
    conda_env = 'n/a'

# Get nicely wrapped command line
cmdline = sys.argv[0]
for i in range(1, len(sys.argv)):
    if sys.argv[i].startswith('--'):
        cmdline += '\n' + (' ' * indent) + sys.argv[i]
    else:
        cmdline += ' ' + sys.argv[i]

# Get abs paths of all config files
cfgfiles = []
for cfg in workflow.configfiles:
    cfgfiles.append(os.path.abspath(cfg))
cfgfiles = '\n                        '.join(cfgfiles)
      


logger.info('=====================================================================================')
logger.info(r'           _ _     _____ ____  _             ')
logger.info(r'  ___  ___(_) |   |___ /|  _ \(_)_ __   ___  ')
logger.info(r' / __|/ __| | |     |_ \| |_) | |  _ \ / _ \ ')
logger.info(r' \__ \ (__| | |___ ___) |  __/| | |_) |  __/ ')
logger.info(r' |___/\___|_|_____|____/|_|   |_|  __/ \___| ')
logger.info(r'                                |_|          ')
logger.info('')
logger.info('    Date:               ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
logger.info('    Platform:           ' + pltfrm)
logger.info('    Host:               ' + hostname)
logger.info('    User:               ' + username)
logger.info('    Conda:              ' + str(conda_ver))
logger.info('    Mamba:              ' + str(mamba_ver))
logger.info('    Python:             ' + str(sys.version.split(' ')[0]))
logger.info('    Snakemake:          ' + str(snakemake.__version__))
logger.info('    sciL3Pipe:          ' + str(scil3pipe_version))
logger.info('    Conda env:          ' + str(conda_env))
logger.info('    Command:            ' + cmdline)
logger.info('')
logger.info('    Base directory:     ' + workflow.basedir)
logger.info('    Working directory:  ' + os.getcwd())
logger.info('    Config file(s):     ' + cfgfiles)
logger.info('    Out directory:      ' + out_dir)
logger.info('    Additional info:    ' + INFO)
logger.info('')
logger.info('=====================================================================================')
logger.info('')

# No need to have these output vars available in the rest of the snakefiles
del indent
del pltfrm, hostname, username
del conda_ver, conda_env
del cmdline, cfgfiles
