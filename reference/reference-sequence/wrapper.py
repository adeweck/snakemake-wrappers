__author__ = "Johannes Köster"
__copyright__ = "Copyright 2019, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import subprocess as sp
import os
import sys
import pathlib
import re
from itertools import product
from snakemake.shell import shell


outfile = snakemake.output[0]


species = re.sub(r"[^a-zA-Z]+", '', snakemake.params.species.lower())
release = int(re.sub(r"^[^0-9]*", "", snakemake.params.release))
source = snakemake.params.source.lower()
datatype = snakemake.params.get("datatype", "")


outfile_suffixes = pathlib.Path(outfile).suffixes
if outfile_suffixes[-1].strip(".") == "gz":
    unzip_cmd = ""
else:
    unzip_cmd = "| gzip -d"
    
    
if source == "gencode":
    
    species_dict = {"homosapiens":"human",
                    "musmusculus":"mouse"}
    species = species_dict.get(species,species)
    release = f"M{release}" if species == "mouse" else release
    
    url_root = f"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{species}/release_{release}"
    
    if datatype == "dna":
        filename = f"GRC{species[0]}{{build_version}}.primary_assembly.genome.fa.gz"
    elif datatype == "cdna":
        filename = f"gencode.v{release}.transcripts.fa.gz"
    elif datatype == "cds":
        filename = f"gencode.v{release}.pc_translations.fa.gz"
    elif datatype == "ncrna":
        filename = f"gencode.v{release}.lncRNA_transcripts.fa.gz"
    else:
        raise ValueError("invalid datatype, must be one of dna, cdna, cds, ncrna")
    
    url0 = os.path.join(url_root, filename)
    build_version_start = 38 if release > 19 else 37 # !! different version between human and mouse. Dodgy heuristics!
        
    
elif source == "ensembl":
    
    species_dict = {"homosapiens":"homo_sapiense",
                    "human":"homo_sapiens",
                    "musmusculus":"mus_musculus",
                    "mouse":"mus_musculus"}
    species = species_dict.get(species,species)
    species_cap=species.capitalize()
    
    url_root = f"ftp://ftp.ensembl.org/pub/release-{release}/fasta/{species}/{datatype}"
    filename_root = f"{species_cap}.GRC{species[0]}{{build_version}}"
    if datatype == "dna":
        filename = f"{filename_root}.dna.primary_assembly.fa.gz"
    elif datatype == "cdna":
        filename = f"{filename_root}.cdna.all.fa.gz"
    elif datatype == "cds":
        filename = f"{filename_root}.cds.all.fa.gz"
    elif datatype == "ncrna":
        filename = f"{filename_root}.ncrna.fa.gz"
    elif datatype == "pep":
        filename = f"{filename_root}.pep.all.fa.gz"
    else:
        raise ValueError("invalid datatype, must be one of dna, cdna, cds, ncrna, pep")
        
    url0 = os.path.join(url_root, filename)
    build_version_start = 38 if release > 75 else 37 # Dodgy heuristics!
        

log = snakemake.log_fmt_shell(stdout=False, stderr=True)        
        

success = False
for build_version in range(build_version_start,build_version_start+5):
    url = url0.format(build_version=build_version)

    try:
        shell("curl -sSf {url} > /dev/null 2> /dev/null")
    except sp.CalledProcessError:
        continue

    shell("(curl -L {url} {unzip_cmd} > {snakemake.output[0]}) {log}")
    success = True
    break

if not success:
    print(
        "Unable to download requested sequence data from Ensembl. "
        "Did you check that this combination of species, build, and release is actually provided?",
        file=sys.stderr,
    )
    exit(1)

