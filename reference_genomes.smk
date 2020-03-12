# from Bio import Entrez
Entrez.email="knutdrand@gmail.com"
entrez_entries = ["NC_017519.1", "NC_014921.1", "NC_010163.1", "NC_013511.1"]

rule myco_genomes:
    input:
        expand("../../Data/Entrez/{entry}.fa.sa", entry=entrez_entries)

# rule get_entrez_fasta:
#     output:
#         "../../Data/Entrez/{entrez_id}.{version}.fa"
#     run:
#         handle = Entrez.efetch(db="nuccore", id=f"{wildcards.entrez_id}.{wildcards.version}", rettype="fasta", retmode="text")
#         with open(output[0], "w") as f:
#             f.write(handle.read())

rule bwa_index:
    input:
        "{name}.fa"
    output:
        multiext("{name}.fa", ".ann", ".amb", ".bwt", ".sa", ".pac")
    shell:
        "bwa index {input}"

rule get_univec:
    output:
        "../../Data/univec/Univec.fa"
    shell:
        "wget https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec -O {output}"

rule download_reference:
    output:
        "../../Data/{species}/{species}.fa.gz"
    shell:
        "wget http://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.species}/bigZips/{wildcards.species}.fa.gz -O {output}"
