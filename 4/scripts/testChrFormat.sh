# Script for determining whether BAM files are NCBI ('1') or UCSF ('chr1') aligned
# Peeks at first number lines for occurances of 'chr'.
# run on mmlab server

cd /data/koplev01/projects/OVCT/TCGA

bam_files='/data/koplev01/projects/OVCT/TCGA/samples/sample_paths.txt'
output_file='bamfile_info.csv'

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