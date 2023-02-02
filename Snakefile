times = {
    "vial_2" : [2, 3, 4, 13, 14],
    "vial_7" : [1, 10, 11, 12, 13, 14],
    "vial_8" : [1, 10, 11, 12, 13, 14],
    "vial_10" : [1, 2, 3, 4, 13, 14],
    "vial_11" : [1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14]
}
vials = list(times.keys())

rule minimap:
    input:
        fastq="data/2022-09-21_amoxicillin_run/{vial_n}/time_{t}/reads.fastq.gz",
        genome="data/2022-09-21_amoxicillin_run/assembled_genome/genome.fa",
    output:
        sam="results/2022-09-21_amoxicillin_run/{vial_n}/sam_files/time_{t}.sam",
    shell:
        """
        minimap2 -a -x map-ont -t 2 {input.genome} {input.fastq} > {output.sam}
        """

rule extract_suppl_info:
    input:
        sam=rules.minimap.output.sam,
    output:
        csv="results/2022-09-21_amoxicillin_run/{vial_n}/suppl_reads_info/time_{t}.csv",
        png="results/2022-09-21_amoxicillin_run/{vial_n}/figures/suppl_reads_info/time_{t}.png"
    shell:
        """
        python3 scripts/extract_suppl_reads_info.py --sam {input.sam} --csv {output.csv} --png {output.png}
        """

rule plot_read_positions:
    input:
        sam=rules.minimap.output.sam,
    output:
        png="results/2022-09-21_amoxicillin_run/{vial_n}/figures/read_positions/time_{t}.png",
    shell:
        """
        python3 scripts/plot_read_positions.py --sam {input.sam} --png {output.png}
        """

rule extract_suppl_freq_by_time:
    input:
        sam=rules.minimap.output.sam,
    output:
        csv="results/2022-09-21_amoxicillin_run/{vial_n}/suppl_freq/time_{t}.csv"
    shell:
        """
        python3 scripts/extract_supp_by_pos.py --sam {input.sam} --csv {output.csv} --header time_{wildcards.t}
        """

rule concat_suppl_freq:
    input:
        lambda w: expand(
            rules.extract_suppl_freq_by_time.output.csv,
            t=times[w.vial_n],
            vial_n=w.vial_n
        ),
    output:
        csv="results/2022-09-21_amoxicillin_run/{vial_n}/suppl_freq/trajectory_{vial_n}.csv"
    run:
        import pandas
        df = []
        for f in input:
            df.append(pandas.read_csv(f,index_col=0))
        c = pandas.concat(df, axis=1)
        c.to_csv(output.csv)

rule plot_traj:
    input:
        csv=rules.concat_suppl_freq.output.csv
    output:
        png="results/2022-09-21_amoxicillin_run/{vial_n}/figures/suppl_freq_trajectory_{vial_n}.png"
    shell:
        """
        python3 scripts/plot_trajectories.py --csv {input.csv} --png {output.png} --n 7
        """

rule all:
    input:
        [expand(
            rules.extract_suppl_info.output,
            t=times[v],
            vial_n=v
        ) for v in vials],
        [expand(
            rules.plot_read_positions.output.png,
            t=times[v],
            vial_n=v
        ) for v in vials],
        expand(
            rules.plot_traj.output.png,
            vial_n=vials
        )
