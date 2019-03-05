#!/bin/sh

if [ -x "$(command -v wget)" ]; then
	echo "====================== Downloading test files =====================" 1>&2
	wget https://zenodo.org/record/2582855/files/lambda.phage_test_sample.bam
	wget https://zenodo.org/record/2582855/files/lambda.phage_test_sample.bam.bai
	wget https://zenodo.org/record/2582855/files/lambda.phage_test_sample_1.fq.gz
	wget https://zenodo.org/record/2582855/files/lambda.phage_test_sample_2.fq.gz
	wget https://zenodo.org/record/2582855/files/lambda_phage.fa
	wget https://zenodo.org/record/2582855/files/lambda_phage.fa.fai
	echo "======================== Download complete ========================" 1>&2
elif [ -x "$(command -v curl)" ]; then
	echo "====================== Downloading test files =====================" 1>&2
	curl -O https://zenodo.org/record/2582855/files/lambda.phage_test_sample.bam
	curl -O https://zenodo.org/record/2582855/files/lambda.phage_test_sample.bam.bai
	curl -O https://zenodo.org/record/2582855/files/lambda.phage_test_sample_1.fq.gz
	curl -O https://zenodo.org/record/2582855/files/lambda.phage_test_sample_2.fq.gz
	curl -O https://zenodo.org/record/2582855/files/lambda_phage.fa
	curl -O https://zenodo.org/record/2582855/files/lambda_phage.fa.fai
	echo "======================== Download complete ========================" 1>&2
else
	cat 1>&2 <<EOF
ERROR: neither curl nor wget are available, cannot download test files automatically.

Please manually download test files from
	
	https://zenodo.org/record/2582855#.XH1hgZyny00

All files should be in the tests/ directory.

EOF
fi