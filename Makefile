.PHONY: all

USR := molari0000
SRV := transfer.scicore.unibas.ch
FLD := /scicore/home/neher/GROUP/data/2022_morbidostat_shared_data/2022-09-21_neher_amoxicillin_run_september

SRC := $(USR)@$(SRV):$(FLD)

DT  := data/2022-09-21_amoxicillin_run
VL  := vial_11
DST := $(DT)/$(VL)
TMS := 1 3 5 7 11 13

# available times for vial 11: 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14


$(DST)/time_%/reads.fastq.gz:
	mkdir -p $(DST)/time_$*
	scp $(SRC)/$(VL)/time_$*/reads.fastq.gz $@

$(DST)/assembled_genome/genome.fa:
	mkdir -p $(DST)/assembled_genome
	scp $(SRC)/$(VL)/time_1/assembled_genome/onculture.fna $@
	
$(DST)/assembled_genome/genome.gbk:
	mkdir -p $(DST)/assembled_genome
	scp $(SRC)/$(VL)/time_1/assembled_genome/onculture.gbk $@

genomes-tgt := $(DST)/assembled_genome/genome.fa $(DST)/assembled_genome/genome.gbk

time-tgt := $(addprefix $(DST)/time_, $(TMS))
time-tgt := $(addsuffix /reads.fastq.gz, $(time-tgt))

all: $(genomes-tgt)  $(time-tgt)
