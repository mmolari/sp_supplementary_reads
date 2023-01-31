times = [1, 3, 5, 7, 11, 13]


rule minimap:
    input:
        fastq="data/2022-09-21_amoxicillin_run/vial_11/time_{t}/reads.fastq.gz",
        genome="data/2022-09-21_amoxicillin_run/vial_11/assembled_genome/genome.fa",
    output:
        sam="results/2022-09-21_amoxicillin_run/vial_11/sam_files/time_{t}.sam",
    shell:
        """
        minimap2 -a -x map-ont -t 2 {input.genome} {input.fastq} > {output.sam}
        """

rule extract_suppl_info:
    input:
        sam=rules.minimap.output.sam,
    output:
        csv="results/2022-09-21_amoxicillin_run/vial_11/suppl_reads_info/time_{t}.csv",
        png="results/2022-09-21_amoxicillin_run/vial_11/figures/suppl_reads_info/time_{t}.png"
    shell:
        """
        python3 scripts/extract_suppl_reads_info.py --sam {input.sam} --csv {output.csv} --png {output.png}
        """

rule plot_read_positions:
    input:
        sam=rules.minimap.output.sam,
    output:
        png="results/2022-09-21_amoxicillin_run/vial_11/figures/read_positions/time_{t}.png",
    shell:
        """
        python3 scripts/plot_read_positions.py --sam {input.sam} --png {output.png}
        """

rule all:
    input:
        expand(
            rules.minimap.output.sam,
            t=times,
        ),
        expand(
            rules.extract_suppl_info.output.csv,
            t=times,
        ),
        expand(
            rules.extract_suppl_info.output.png,
            t=times,
        ),
        expand(
            rules.plot_read_positions.output.png,
            t=times,
        ),
