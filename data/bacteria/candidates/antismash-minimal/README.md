# antiSMASH Rule-Based (minimal)

Run on HPC using:

```bash
cd antismash-minimal
for file in ../../genome/*.fa; do
    qsub -V -b y -cwd -l mem_reserve=10G -l mem_free=10G -l h_vmem=10G -pe threads 4 \
        antismash $file --minlength 1 --minimal --disable-genbank --disable-embl --disable-html \
        --disable-svg --disable-xls --logfile $(basename "$file").log;
done;
```

Generate candidates.tsv using:

```bash
cd antismash-minimal
echo "Contig	BGC ID	BGC type	detection rules used	BGC_range	genes	subclusters	NRPSs/PKSs	signature_genes	RiPPs	predicted structure	monomers" > candidates.tsv
for file in *.genome; do
    name=$(basename "$file");
    cat $file/txt/*_BGC.txt | grep -v "detection rules used" | awk '{ print "'${name/.genome/}'\t"$0}' >> candidates.tsv;  
done
```