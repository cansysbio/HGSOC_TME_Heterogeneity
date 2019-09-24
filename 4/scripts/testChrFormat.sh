# Script for determining whether BAM files are NCBI ('1') or UCSF ('chr1') aligned
# Peeks at first number lines for occurances of 'chr'.
# run on mmlab server

source activate titan

# Generate list of bam files with absolute paths
bam_files='/data/koplev01/projects/OVCT/TCGA/samples/sample_paths.txt'

cd /data/memon01/tcga/TCGA-OV-LEGACY
find "$PWD" -name "*.bam" > $bam_files

cd /data/memon01/tcga/TCGA-OV-LEGACY-ADD
find "$PWD" -name "*.bam" >> $bam_files

cd /data/memon01/tcga/TCGA-OV-LEGACY-WITH-WGS
find "$PWD" -name "*.bam" >> $bam_files


cd /data/koplev01/projects/OVCT/TCGA
output_file='bamfile_info_add_wgs.csv'

# Write heade
echo "path,chr_count,header"  > $output_file  # overwrites

while read path; do
	echo $path
	printf $path >> $output_file  # no newline

	printf "," >> $output_file

	samtools view $path | head -n 10000 | grep -c "chr" | tr -d '\n' >> $output_file

	printf "," >> $output_file

	# Get column 4, second line
	samtools view -H $path | sed -n '2p' | tr -d '\n' >> $output_file  # second line of bam header

	printf "\n" >> $output_file

done < $bam_files