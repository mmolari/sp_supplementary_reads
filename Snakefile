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


rule all:
    input:
        expand(
            "results/2022-09-21_amoxicillin_run/vial_11/sam_files/time_{t}.sam",
            t=times,
        ),
