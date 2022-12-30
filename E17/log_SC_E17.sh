## Check md5sums

## Note: Cellranger counts only counts exonic reads by default: https://github.com/satijalab/seurat/issues/2472

## Aligning resequenced samples seperately:

## Using prebuilt reference: https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#mm10_2020A

cd /mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/E17/aligned_reads/separate_runs

ref="/home/shared/hg_align_db/GRCm38_Cellranger_gencode_vM23_ensembl_98"

cellranger=/opt/cellranger-6.1.2/bin/cellranger

for run_date in ("190322" "190501" "210604" "210614"); do 
  samples=$(cat $run_date/SampleSheet.csv | grep .*-GA | sed 's/,.*//')
  for ea in $samples; do
    sample_id=$(ls $run_date | grep $ea)
    fastq_path="$run_date/$sample_id"
     $cellranger count \
    --id=$sample_id \
    --fastqs=$fastq_path \
    --sample=$ea \
    --transcriptome=$ref \
    --localcores=13
  done;
done;

## QC:

cd /mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/E17/aligned_reads/separate_runs

samples=$(ls | grep .*Run.* | sed "s/_Run1//" | sed "s/_Run2//" | uniq)
for ea in $samples; do
  $cellranger aggr \
  --id="${ea}_aggregated" \
  --csv="${ea}_aggr.csv" \
  --normalize=none
  done;
done;

## Note: Ncre_MADM_EGFRF_Td runs show some batch effects, i.e. substantially different depth between runs.

## Combining resequenced samples:

cd /mnt/bdata/rebecca/collabs/ghashghaei/SC_expr/E17/aligned_reads/combined_runs

cellranger=/opt/cellranger-6.1.2/bin/cellranger

ref="/home/shared/hg_align_db/GRCm38_Cellranger_gencode_vM23_ensembl_98"

samples=$(cat ../../raw_reads/190322/SampleSheet.csv | grep .*-GA | sed 's/,.*//')
echo $samples
for ea in $samples; do
  run1=$(ls ../../raw_reads/190322/ | grep $ea)
  run2=$(ls ../../raw_reads/190501/ | grep $ea)
  $cellranger count \
    --id=$ea \
    --fastqs=../../raw_reads/190322/$run1,../../raw_reads/190501/$run2 \
    --sample=$ea \
    --force-cells 50000 \
    --transcriptome=$ref \
    --localcores=13 \
    --no-bam
done;

samples=$(cat ../../raw_reads/210604/SampleSheet.csv | grep .*-GA | sed 's/,.*//')
echo $samples
for ea in $samples; do
  run1=$(ls ../../raw_reads/210604/ | grep $ea)
  run2=$(ls ../../raw_reads/210614/ | grep $ea)
  $cellranger count \
   --id=$ea \
   --fastqs=../../raw_reads/210604/$run1,../../raw_reads/210614/$run2 \
   --sample=$ea \
   --force-cells 50000 \
   --transcriptome=$ref \
   --localcores=13 \
   --no-bam
done;

## Samples that weren't resequenced:

samples=$(cat ../../raw_reads/201210/SampleSheet.csv | grep .*-GA | sed 's/,.*//')
echo $samples
for ea in $samples; do
  $cellranger count \
    --id=$ea \
    --fastqs=../../raw_reads/201210/$ea \
    --sample=$ea \
    --force-cells 50000 \
    --transcriptome=$ref \
    --localcores=13 \
    --no-bam
done;

