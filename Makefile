.PHONY: all

USR := molari0000
SRV := transfer.scicore.unibas.ch
FLD := /scicore/home/neher/GROUP/data/2022_morbidostat_shared_data/2022-09-21_neher_amoxicillin_run_september

SRC := $(USR)@$(SRV):$(FLD)

DT  := data/2022-09-21_amoxicillin_run
VL  := vial_8
DST := $(DT)/$(VL)
TMS := 1 10 11 12 13 14
GN := barcode01_genome


$(DST)/time_%/reads.fastq.gz:
	mkdir -p $(DST)/time_$*
	scp $(SRC)/$(VL)/time_$*/reads.fastq.gz $@

$(DT)/assembled_genome/genome.fa:
	mkdir -p $(DT)/assembled_genome
	scp $(SRC)/onculture/assembled_genome/$(GN).fna $@
	
$(DT)/assembled_genome/genome.gbk:
	mkdir -p $(DT)/assembled_genome
	scp $(SRC)/onculture/assembled_genome/$(GN).gbk $@

genomes-tgt := $(DT)/assembled_genome/genome.fa $(DT)/assembled_genome/genome.gbk

time-tgt := $(addprefix $(DST)/time_, $(TMS))
time-tgt := $(addsuffix /reads.fastq.gz, $(time-tgt))

all: $(genomes-tgt)  $(time-tgt)
