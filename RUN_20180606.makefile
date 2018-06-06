all:

SRC := BBJ_SUMMARY

LD := LD/fourier_ls-all.bed
nLD := $(shell cat $(LD) | tail -n+2 | wc -l)

GWAS := $(wildcard $(SRC)/*.autosome.txt.gz)

DATA := gwas_data/

TRAITS := $(shell ls -1 $(SRC)/*.autosome.txt.gz | xargs -I file basename file .autosome.txt.gz | sed 's/BBJ.//g')

TRAITS_STR := $(shell ls -1 $(SRC)/*.autosome.txt.gz | xargs -I file basename file .autosome.txt.gz | sed 's/BBJ.//g' | awk 'NR == 1 { out = $$1 } NR > 1 { out = out ":" $$1 } END { print out }')

################################################################

queue_data: $(foreach tr, $(TRAITS), jobs/20180606/data-$(tr).txt.gz)

jobs/20180606/data-%.txt.gz: $(SRC)/BBJ.%.autosome.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	echo ./make.data-ldblock.R $(LD) $< $(DATA) | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=8g -l h_rt=2:00:00 -b y -j y -N BBJ_DATA_$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

################################################################

queue_facto: $(foreach k, 20 50, jobs/20180606/facto-$(k).txt.gz)

jobs/20180606/facto-%.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(nLD) | awk '{ print "./make.factorization.R" FS $$1 FS "$(LD)" FS "$(DATA)" FS "$(TRAITS_STR)" FS $* FS ("result/20180606/factorization/$*/" $$1) }' | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=8g -l h_rt=2:00:00 -b y -j y -N BBJ_FACT_$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

################################################################

queue_conf: $(foreach k, 20 50, jobs/20180606/conf-$(k).txt.gz)

jobs/20180606/conf-%.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(nLD) | awk '{ print "./make.confounder.R" FS $$1 FS "$(LD)" FS "result/20180606/factorization/$*/" FS ("result/20180606/confounder/$*/" $$1 ".txt.gz") }' | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=4g -l h_rt=2:00:00 -b y -j y -N BBJ_CONF_$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

################################################################

queue_fgwas: $(foreach k, 20 50, jobs/20180606/fgwas-$(k).txt.gz)

jobs/20180606/fgwas-%.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(nLD) | awk '{ print "./make.fgwas.R" FS $$1 FS "$(LD)" FS "$(DATA)" FS "$(TRAITS_STR)" FS ("result/20180606/confounder/$*/" $$1 ".txt.gz") FS $* FS ("result/20180606/fgwas/$*/" $$1) }' | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=4g -l h_rt=2:00:00 -b y -j y -N BBJ_GWAS_$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@
