all:

SRC := BBJ_SUMMARY

LD := LD/fourier_ls-all.bed
nLD := $(shell cat $(LD) | tail -n+2 | wc -l)
LODD := 0

################################################################
queue_qc: $(foreach k, 50, jobs/20180615/qc-$(k).txt.gz)

jobs/20180615/qc-%.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(nLD) | awk '{ print "./make.fgwas-qc.R" FS $$1 FS "$(LD)" FS "result/20180613/fgwas/$*/" FS "$(LODD)" FS ("result/20180615/fgwas-qc/$*/" $$1 ".txt.gz") }' | gzip > $@
	seq 1 $(nLD) | awk '{ print "./make.fgwas-qc.R" FS $$1 FS "$(LD)" FS "result/20180613/fgwas_nn/$*/" FS "$(LODD)" FS ("result/20180615/fgwas_nn-qc/$*/" $$1 ".txt.gz") }' | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=4g -l h_rt=2:00:00 -b y -j y -N BBJ_QC_$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

################################################################
pool: $(foreach k, 50, result/20180615/fgwas-qc-$(k).txt.gz result/20180615/fgwas_nn-qc-$(k).txt.gz)

result/20180615/fgwas-qc-%.txt.gz:
	ls -1 result/20180615/fgwas-qc/$*/*.txt.gz | xargs zcat | awk 'NR == 1 || $$1 != "chr"' | gzip > $@

result/20180615/fgwas_nn-qc-%.txt.gz:
	ls -1 result/20180615/fgwas_nn-qc/$*/*.txt.gz | xargs zcat | awk 'NR == 1 || $$1 != "chr"' | gzip > $@
