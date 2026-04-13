"""
1. Download Blast from NCBI
2. Download UniRef90 database AND Build Blast index by Blast
"""

import os


def download_blast(download_folder, my_logger, version="2.16.0"):
    """
    :param download_folder: where NCBI blast was downloaded
    :param my_logger:
    :param version:
    :return:
    """
    os.makedirs(download_folder, exist_ok=True)
    os.chdir(download_folder)

    # Download
    blast_web_link = f"https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/{version}/ncbi-blast-{version}+-x64-linux.tar.gz"
    print(f"Download \n{blast_web_link}\nto the folder\n{download_folder} ...")
    my_logger.info(f"Download \n{blast_web_link}\nto the folder\n{download_folder} ...")
    os.system(f"wget -b -c {blast_web_link}")

    # Check Download
    if os.path.exists(os.path.join(download_folder, f"ncbi-blast-{version}+-x64-linux.tar.gz")):
        my_logger.info(f"Downloaded to {download_folder}")
    else:
        print(f"Fail to download {blast_web_link} to {download_folder}")
        return None

    # decompress
    print(f"Decompress ncbi-blast-{version}+-x64-linux.tar.gz ...")
    my_logger.info(f"Decompress ncbi-blast-{version}+-x64-linux.tar.gz ...")
    os.system(f"tar -xvzf ncbi-blast-{version}+-x64-linux.tar.gz")

    # check bin folder of blast
    blast_bin_path = os.path.join(download_folder, f"ncbi-blast-{version}+", "bin")
    if os.path.exists(blast_bin_path):
        my_logger.info(f"Success! The bin path of blast is {blast_bin_path}")
        print(f"Success! The bin path of blast is {blast_bin_path}")
        return blast_bin_path
    else:
        print(f"Failed! The bin path of blast {blast_bin_path} not exists!")
        return None


def download_build_uniref90(blast_bin, uniref_download_folder, my_logger):
    """
    :param blast_bin = "/gpfs/data/lilab/home/zhoub03/software/blast/ncbi-blast-2.16.0+/bin"
    :param uniref_download_folder:
    :param my_logger:
    :return:
    """
    os.makedirs(uniref_download_folder, exist_ok=True)
    os.chdir(uniref_download_folder)

    # download current version of UniRef90
    uniref90_web_link = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta.gz"
    uniref90_release_link = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.release_note"
    print(f"Download \n{uniref90_web_link}\nto the folder\n{uniref_download_folder} ...")
    my_logger.info(f"Download \n{uniref90_web_link}\nto the folder\n{uniref_download_folder} ...")
    os.system(f"wget -b -c {uniref90_release_link}")
    os.system(f"wget -b -c {uniref90_web_link}")

    # Check Download
    if os.path.exists(os.path.join(uniref_download_folder, f"uniref90.fasta.gz")):
        my_logger.info(f"Downloaded uniref90.fasta.gz to {uniref_download_folder}")
    else:
        print(f"Fail to download uniref90.fasta.gz to {uniref_download_folder}")
        return None

    # unzip uniref90.fasta.gz
    print(f"Unzip uniref90.fasta.gz to {uniref_download_folder}")
    os.system(f"gunzip -c uniref90.fasta.gz > uniref90.fasta")

    # make blast db for UniRef90
    os.environ["PATH"] = os.environ["PATH"] + ":" + blast_bin
    makedb_command = "makeblastdb -in uniref90.fasta -dbtype prot -parse_seqids"
    print(makedb_command)
    my_logger.info(makedb_command)
    os.system(makedb_command)

    # check pdb of blast
    if os.path.exists(os.path.join(uniref_download_folder, "uniref90.fasta.pdb")):
        print(f"Blastdb of UniRef90 has been built. The pre-construct step is completed.")
        my_logger.info(f"Blastdb of UniRef90 has been built. The pre-construct step is completed.")
        return 1
    else:
        print(f"Failed to make blastdb for UniRef90")
        my_logger.info(f"Failed to make blastdb for UniRef90")
        return None
