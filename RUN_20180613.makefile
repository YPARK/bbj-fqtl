all:

SRC := BBJ_SUMMARY

LD := LD/fourier_ls-all.bed
nLD := $(shell cat $(LD) | tail -n+2 | wc -l)

GWAS := $(wildcard $(SRC)/*.autosome.txt.gz)

DATA := gwas_data/

TRAITS := $(shell ls -1 $(SRC)/*.autosome.txt.gz | xargs -I file basename file .autosome.txt.gz | sed 's/BBJ.//g')

TRAITS_STR := $(shell ls -1 $(SRC)/*.autosome.txt.gz | xargs -I file basename file .autosome.txt.gz | sed 's/BBJ.//g' | awk 'NR == 1 { out = $$1 } NR > 1 { out = out ":" $$1 } END { print out }')

################################################################

queue_fgwas: $(foreach k, 50, jobs/20180613/fgwas-$(k).txt.gz)

queue_fgwas_long: $(foreach k, 50, jobs/20180613/fgwas-$(k)-long.gz)

jobs/20180613/fgwas-%.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(nLD) | awk '{ print "./make.fgwas.R" FS $$1 FS "$(LD)" FS "$(DATA)" FS "$(TRAITS_STR)" FS ("result/20180606/confounder/$*/" $$1 ".txt.gz") FS $* FS ("result/20180613/fgwas/$*/" $$1) FS "FALSE" }' | gzip > $@
	seq 1 $(nLD) | awk '{ print "./make.fgwas.R" FS $$1 FS "$(LD)" FS "$(DATA)" FS "$(TRAITS_STR)" FS ("result/20180606/confounder/$*/" $$1 ".txt.gz") FS $* FS ("result/20180613/fgwas_nn/$*/" $$1) FS "TRUE" }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=4g -l h_rt=2:00:00 -b y -j y -N BBJ_GWAS_$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

jobs/20180613/fgwas-%-long.gz: jobs/20180613/fgwas-%.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@zcat $< | awk 'system("! [ -f " $$8 ".zscore.gz ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -gt 0 ] && qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=4g -l h_rt=16:00:00 -b y -j y -N BBJ_LONG_FGWAS -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@


################################################################
queue_summary: $(foreach k, 50, jobs/20180613/summary-$(k).txt.gz)

jobs/20180613/summary-%.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	echo "./make.combine-fgwas.R result/20180613/fgwas/$*/ $(LD) 0.5 result/20180613/fgwas-$*-05.txt.gz" | gzip > $@
	echo "./make.combine-fgwas.R result/20180613/fgwas/$*/ $(LD) 0.9 result/20180613/fgwas-$*-09.txt.gz" | gzip >> $@
	echo "./make.combine-fgwas.R result/20180613/fgwas_nn/$*/ $(LD) 0.5 result/20180613/fgwas_nn-$*-05.txt.gz" | gzip >> $@
	echo "./make.combine-fgwas.R result/20180613/fgwas_nn/$*/ $(LD) 0.9 result/20180613/fgwas_nn-$*-09.txt.gz" | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=8g -l h_rt=2:00:00 -b y -j y -N BBJ_SUM_$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@


