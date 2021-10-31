__author__ = "Johannes Köster"
__copyright__ = "Copyright 2019, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "BSD3"


import subprocess
import sys
import pathlib
import re
from snakemake.shell import shell


outfile = snakemake.output[0]


species = re.sub(r"[^a-zA-Z]+", '', snakemake.params.species.lower())
release = int(re.sub(r"^[^0-9]*", "", snakemake.params.release))
source = snakemake.params.source.lower()
fmt = snakemake.params.get("fmt", "")
build = snakemake.params.get("build", "")
flavor = snakemake.params.get("flavor", "")


suffixes = pathlib.Path(outfile).suffixes
if fmt == "":
    fmt = suffixes[0].strip(".") #gtf or gff3 expected
if suffixes[-1].strip(".") == "gz":
    unzip_cmd = ""
else:
    unzip_cmd = "| gzip -d"

    
if flavor:
    flavor += "."
    
    
if source == "gencode":
    
    species_dict = {"homosapiens":"human",
                    "musmusculus":"mouse"}
    species = species_dict.get(species,species)
    release = f"M{release}" if species == "mouse" else release
    url = f"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{species}/release_{release}/gencode.v{release}.{flavor}annotation.{fmt}.gz"
    
elif source == "ensembl":
    
    species_dict = {"homosapiens":"homo_sapiense",
                    "human":"homo_sapiens",
                    "musmusculus":"mus_musculus",
                    "mouse":"mus_musculus"}
    species = species_dict.get(species,species)
    species_cap=species.capitalize()
    
    if build == "":
        # check if file exists
        build_core = "GRC{}".format(species[0])
        url_tmp = f"ftp://ftp.ensembl.org/pub/release-{release}/{fmt}/{species}/{species_cap}.{build_core}{{build_version}}.{release}.{flavor}{fmt}.gz"
        build_version = 38 if release > 75 else 37
        nomatch = True
        iter = 0
        while nomatch and iter < 5:
            url = url_tmp.format(build_version=build_version+iter)
            iter+=1
            try:
                nomatch = False
                subprocess.check_output(["curl", "--head", url])
            except:
                nomatch = True
    else:
        url = f"ftp://ftp.ensembl.org/pub/release-{release}/{fmt}/{species}/{species_cap}.{build}.{release}.{flavor}{fmt}.gz"
    
    
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

    
print(f"Downloading annotation from: {url}")

try:
    shell("(curl -L {url} {unzip_cmd} > {outfile}) {log}")
except subprocess.CalledProcessError as e:
    if snakemake.log:
        sys.stderr = open(snakemake.log[0], "a")
    print(
        "Unable to download annotation data from Ensembl. "
        "Did you check that this combination of species, build, and release is actually provided?",
        file=sys.stderr,
    )
    exit(1)
