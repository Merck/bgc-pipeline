AntiSMASH command:

```
mkdir classification && cd classification && for file in ../../genome-contigs/*; do
  qsub -V -b y -cwd antismash $file --minlength 1 --disable-genbank
  --disable-embl --disable-html --disable-svg --disable-xls --logfile $(basename "$file").log
  --debug; done && ../../../../bgc_detection/hpc/wait_all.sh
```

Merging:

```
printf "Contig\t" > classification.tsv; cat classification/*/txt/*_BGC.txt |
  head -1 >> classification.tsv; for file in classification/*/; do name=$(basename
  "$file"); cat $file/txt/*_BGC.txt | grep -v "detection rules used" | awk '{ print
  "'$name'\t"$0 }' >> classification.tsv; done
```
