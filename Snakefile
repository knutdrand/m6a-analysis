include: "reference_genomes.smk"
import pandas as pd
data_path = "../../Data/"
base_path = "/media/knut/KnutData/m6a-data/output/"
read_path = "all_reads/" # "/media/knut/KnutData/m6a-data/200224_NB501273.Project_Li-libs2-2020-02-17/"
workingdir = ""# "/media/knut/KnutData/m6a-data/"
# read_path = base_path + "200224_NB501273.Project_Li-libs2-2020-02-17/"
bam_path = "/media/knut/KnutData/m6a-data/output/"
# bam_path = base_path+"output/"
chromosome_grep = "grep -Ew -e 'chr[0-9]{{1,2}}' -e chrX -e chrY"

#samples = [s.strip().replace(".fastq.gz", "") for s in open("config.csv")]
samples = [s.replace(".fastq.gz", "") for s in pd.read_csv("sherif_nels.csv").set_index("filename").index]
# samples = [s for s in samples if int(s.split("-")[0])<=6]
def get_species(sample_name):
    return "mm10"
    if int(sample_name.split("-")[0])<=6:
        return "mm10"
    return "danRer11"

species_dict = {"mm10": [], "danRer11": []}

for sample_name in samples:
    species_dict[get_species(sample_name)].append(sample_name)

def pyplot(code):
    return """python -c "import numpy as np; import matplotlib.pyplot as plt; %s" """ % code

rule trackhub:
    input:
        lambda wildcards: expand("trackhub/{{species}}/{folder}_{name}.bw", name=species_dict[wildcards.species], folder=("dedup", "mapped_reads"))
    output:
        "trackhub/{species}/trackDb.txt"
    shell:
        'chiptools trackdb single {input} > {output}'

rule all_zebra:
    input:
        ["{species}/logvplots/{place}/{name}.png".format(species=s, place=p, name=n) for s in ["danRer11"] for p in ["first", "last"] for n in species_dict[s]],

rule zebra_gb:
    input:
        [f"{species}/dedup_coverage/{name}.bw" for species in ["danRer11"] for name in species_dict[species]],
        [f"{species}/mapped_reads_coverage/{name}.sortedb.bw" for species in ["danRer11"] for name in species_dict[species]]
        
rule all:
    input:
        [f"{species}/logvplots/{place}/{name}.png" for species in ["danRer11", "mm10"] for place in ["first", "last"] for name in species_dict[species]],
        [f"{species}/metagene/{name}.png" for species in ["danRer11", "mm10"] for place in ["first", "last"] for name in species_dict[species]],
        # report(expand("fastqscreen/{name}_screen.png", name=samples), category="contamination"),
        # expand("{species}/reads_fig.svg", species=["danRer11"])
        # expand("fastqscreen/{name}_screen.png", name=samples)

rule move_data:
    input:
        read_path + "{samplename}"
    output:
        workingdir + "reads/{samplename}"
        # lambda wildcards: workingdir + get_species(wildcards.samplename) + "/reads/" + wildcards.samplename
    shell:
        "mv {input} {output}"

rule get_mouse_annotation:
    output:
        data_path + "mm10/annotation.gtf.gz"
    shell:
        "wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.annotation.gtf.gz -O {output}"

rule get_zebrafish_annotation:
    output:
        data_path + "danRer11/annotation.gtf.gz"
    shell:
        "wget ftp://ftp.ensembl.org/pub/release-99/gtf/danio_rerio/Danio_rerio.GRCz11.99.gtf.gz -O {output}"

rule gunzip:
    input:
        "{path}.gz"
    output:
        "{path}"
    wildcard_constraints:
        path=".*(fa|bed|gdb|gtf)"
    shell:
        "gunzip {input} --keep"

rule download_genes:
    output:
        "{species}/data/refGene.txt.gz"
    shell:
        "wget https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.species}/database/refGene.txt.gz -O {output}"

rule get_genes_bed:
    input:
        "{species}/data/refGene.txt.gz"
    output:
        "{species}/data/genes.bed"
    shell:
        """z%s {input} | awk '{{OFS="\t"}}{{print $3, $5, $6, ".", ".", $4}}' | uniq > {output}""" % chromosome_grep

rule get_first_exon_bed:
    input:
        "{species}/data/refGene.txt.gz"
    output:
        "{species}/data/first_exons.bed"
    shell:
        """zcat {input} | awk '{{OFS="\t"}}{{split($10, a, ","); split($11, b, ",");print $3, a[1], b[1], ".", ".", $4}}' | sortBed | uniq > {output}"""

rule get_last_exon_bed:
    input:
        "{species}/data/refGene.txt.gz"
    output:
        "{species}/data/last_exons.bed"
    shell:
        """zcat {input} | awk '{{OFS="\t"}}{{n=split($10, a, ","); split($11, b, ",");print $3, a[n-1], b[n-1], ".", ".", $4}}' | sortBed |uniq > {output}"""

rule get_exons:
    input:
        data_path + "{species}/annotation.gtf"
    output:
        data_path + "{species}/hisat/exons"
    shell:
        "hisat2_extract_exons.py {input} > {output}"

rule get_splice_sites:
    input:
        data_path + "{species}/annotation.gtf"
    output:
        data_path + "{species}/hisat/splice_sites"
    shell:
        "hisat2_extract_splice_sites.py {input} > {output}"

rule index_hisat:
    input:
        data_path + "{species}/hisat/splice_sites",
        data_path + "{species}/hisat/exons",
        data_path + "{species}/{species}.fa",
    output:
        data_path + "{species}/hisat/genome.1.ht2"
    shell:
        """
        hisat2-build --ss {input[0]} --exon {input[1]} {input[2]} %s{wildcards.species}/hisat/genome
        """ % data_path

# rule download_hisat_index:
#     output:
#         data_path+"mm10/mm10.tar.gz"
#     shell:
#         "wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/mm10.tar.gz -O {output}"
#         
# rule unpack_hisat_index:
#     input:
#         data_path+"mm10/mm10.tar.gz"
#     output:
#         expand(data_path+"mm10/hisat/genome.{i}.ht2", i=range(1, 9))
#     shell:
#         """
#         tar -xzf {input}
#         mv %smm10/mm10/ %smm10/hisat/
#         """ % (data_path, data_path)

rule map_hisat:
    input:
        data_path + "{species}/hisat/genome.1.ht2",
        workingdir + "all_reads/{name}.fastq.gz"
    output:
        temp(workingdir + "{species}/mapped_reads/{name}.bam")
    threads: 16
    resources:
        mem_gb=10
    shell:
        """hisat2 -x %s{wildcards.species}/hisat/genome -U {input[1]} -p {threads} | samtools view -q 30 -S -b > {output}""" % data_path

rule flagstat:
    input:
        workingdir + "{species}/mapped_reads/{name}.bam"
    output:
        report("{species}/flagstat/{name}.txt")
    shell:
        "samtools flagstat {input} > {output}"

rule get_summaries:
    input:
        lambda wildcards: expand(workingdir+ "{{species}}/{{folder}}/{name}.bam", name=species_dict[wildcards.species])
    output:
        "{species}/{folder}_count.txt"
    wildcard_constraints:
        name="[^/]+",
        folder="[^/]+",
        species="[^/]+"
    shell:
        """
        touch {output}
        for f in {input}
        do
            samtools view $f | wc -l >> {output}
        done
        """

# rule get_unique_summaries:
#     input:
#         expand(workingdir + "{{species}}/mapped_reads/{name}.unique.bed", name=samples)
#     output:
#         "{species}/read_count.unique.txt"
#     shell:
#         "wc -l {input} > {output}"

rule get_read_summary:
    input:
        lambda wildcards: expand(workingdir+"reads/{name}.fastq.gz", name=species_dict[wildcards.species])
    output:
        "{species}/read_count.txt"
    shell:
        """
        touch {output}
        for f in {input}
        do
            zgrep ^@ $f | wc -l >> {output}
        done
        """

rule combine_summaries:
    input:
        expand("{{species}}/{filename}.txt", filename=["read_count", "mapped_reads_count", "dedup_count"])
    output:
        report("{species}/mapping_summary.csv")
    shell:
        """
        paste {input} > {output}
        """
        # paste {input} | awk '{{OFS=","}}{{n=split($3, a, "/");print a[n],$1,$2,$4}}'| head -n -1 >> {output}

rule bamtobed:
    input:
        "{name}.bam"
    output:
        "{name}.bed.gz"
    shell:
        "samtools view {input} | chiptools rnabam2bed | gzip > {output}"

rule sort_bed:
    input:
        "{name}.bed.gz"
    output:
        "{name}.sorted.bed"
    resources:
        mem_gb=10
    shell:
        "bedtools sort -i {input} > {output}"

rule filter_dup:
    input:
        "{name}.sorted.bed"
    output:
        "{name}.unique.bed"
    shell:
        "chiptools filterdup {input} > {output}"

rule repeat_coverage_hack:
    input:
        "{species}/{folder}/{name}.bed.gz",
        "{species}/data/chrom.sizes.txt"
    output:
        "{species}/{folder}_coverage/{name}.bdg"
    wildcard_constraints:
        name="[^/]+",
        folder="[^/]+",
        species="[^/]+"
    shell:
        "z" +chromosome_grep + " {input[0]} | bedtools slop -i - -g {input[1]} -b 0 | bedtools genomecov -bga -i - -g {input[1]} > {output}"

rule download_chrom_sizes:
    output:
        "{species}/data/chromInfo.txt.gz"
    shell:
        "wget https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.species}/database/chromInfo.txt.gz -O {output}"

rule clean_chrom_sizes:
    input:
        "{species}/data/chromInfo.txt.gz"
    output:
        "{species}/data/chrom.sizes.txt"
    shell:
        "z"+chromosome_grep + " {input} > {output}"

rule exon_vplots:
    input:
        workingdir+"{species}/dedup_coverage/{name}.bdg",
        "{species}/data/{place}_exons.bed"
    output:
        "{species}/vplots/{place}/{name}.npy",
        "{species}/vplots/{place}/{name}.png",
    shell:
        "cat {input[0]} | chiptools vplot {input[1]} {output}"

rule metagene:
    input:
        workingdir+"{species}/dedup_coverage/{name}.bdg",
        "{species}/data/refGene.txt.gz"
    output:
        report("{species}/metagene/{name}.png", category="Metagene")
    shell:
        "cat {input[0]} | chiptools metagene {input[1]} {output}"

rule exon_logvplots:
    input:
        "{species}/vplots/{place}/{name}.npy",
    output:
        report("{species}/logvplots/{place}/{name}.png", category="VPlots")
    shell:
        pyplot("plt.imshow(np.log(np.load('{input}')+1));plt.savefig('{output}')")

rule sort_bam:
    input:
        "{name}.bam"
    output:
        "{name}.sortedb.bam"
    shell:
        "samtools sort {input} > {output}"

rule mark_duplicates:
    input:
        workingdir+"{species}/mapped_reads/{sample}.sortedb.bam"
    output:
        bam=workingdir+"{species}/dedup/{sample}.bam",
        metrics="{species}/dedup/{sample}.metrics.txt"
    log:
        "logs/{species}/picard/dedup/{sample}.log"
    params:
        "REMOVE_DUPLICATES=true"
    shell:
        """
        picard MarkDuplicates {params} INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics} &> {log}
        """

rule get_bdg2bw:
    output:
        "src/bdg2bw"
    shell:
        """
        wget https://gist.githubusercontent.com/taoliu/2469050/raw/34476e91ebd3ba9d26345da88fd2a4e7b893deea/bdg2bw -O {output}
        chmod a+x {output}
        """

rule create_bw_track:
    input:
        "src/bdg2bw",
        "{species}/{name}.bdg",
        "{species}/data/chrom.sizes.txt"
    output:
        "{species}/{name}.bw"
    wildcard_constraints:
        species="[^/]+"
    shell:
        "{input}"

rule move_to_trackhub_dedup:
    input:
        "{species}/{folder}_coverage/{name}.bw"
    output:
        "trackhub/{species}/{folder}_{name}.bw"
    wildcard_constraints:
        folder="[^(/,_)]+",
        species="[^/]+"
    shell:
        "mv {input} {output}"

rule move_to_trackhub:
    input:
        "{species}/{folder}_coverage/{name}.sortedb.bw",
    output:
        "trackhub/{species}/{folder}_{name}.bw"
    wildcard_constraints:
        folder="mapped_reads",
        species="[^/]+"
    shell:
        "mv {input} {output}"

rule do_fastq_screen:
    input:
        expand("/media/knut/KnutData/200224_NB501273.Project_Li-libs2-2020-02-17/{name}.fastq.gz", name=samples)
    output:
        report(expand("fastqscreen/{name}_screen.png", name=samples), category="contamination")
    shell:
        "fastq_screen --conf fastq_screen_config.txt --aligner bwa --outdir fastqscreen/ {input}"

rule get_rrnas_bed_file_zebra:
    input:
        "../../Data/danRer11/annotation.gtf"
    output:
        "danRer11/data/rrnas.bed"
    shell:
        """grep -i "rrna" {input} | grep "^[0-9]" | awk '{{OFS="\\t"}}{{print "chr"$1, $4, $5}}' | sortBed -i - | mergeBed -i - > {output}"""

rule get_rrnas_bed_file_mouse:
    input:
        "../../Data/mm10/annotation.gtf"
    output:
        "mm10/data/rrnas.bed"
    shell:
        """grep -i "rrna" {input} | """ + chromosome_grep + """ | awk '{{OFS="\\t"}}{{print $1, $4, $5}}' | sortBed -i - | mergeBed -i - > {output}"""

rule count_rrnas_reads:
    input:
        "{species}/{folder}/{name}.bam",
        "{species}/data/rrnas.bed"
    output:
        "{species}/{folder}_rrna_counts/{name}.txt"
    shell:
        "intersectBed -abam {input[0]} -b {input[1]} -u -bed | wc -l > {output}"

rule summarize_folder:
    input:
        lambda wildcards: expand("{{species}}/dedup_rrna_counts/{name}.txt", name=species_dict[wildcards.species])
    output:
        "{species}/dedup_rrna_counts/ALL.txt"
    shell:
        "cat {input} > {output}"

rule summarize_folder2:
    input:
        lambda wildcards: expand("{{species}}/mapped_reads_rrna_counts/{name}.sortedb.txt", name=species_dict[wildcards.species])
    output:
        "{species}/mapped_reads_rrna_counts/ALL.txt"
    shell:
        "cat {input} > {output}"

rule mapping_stats:
    input:
        "{species}/names.txt",
        "{species}/mapped_reads_count.txt",
        "{species}/dedup_count.txt",
        "{species}/mapped_reads_rrna_counts/ALL.txt",
        "{species}/dedup_rrna_counts/ALL.txt",
    output:
        report("{species}/reads_fig.svg", category="read_counts")
    shell:
        "python src/barchart.py all,dedup,rrna,rrna_dedup {input} {output}"

#rule gzip:
#    input:
#        "{filename}.bed"
#    output:
#        "{filename.bed.gz"
#    shell:
#        "gzip {input}"
