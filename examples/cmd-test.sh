# command to prepare input from SNP positions
../bin/SNP_Pos2polymarker_input.py snp_pos_example.txt polymarker_input_example.csv ./blastdb/test_reference.fa

# to run all at once
../run_getkasp.py polymarker_input_example.csv 200 1 1 63 25 0 ./blastdb/test_reference.fa

## to run step by step
../bin/parse_polymarker_input.py polymarker_input_example.csv
blastn -task blastn -db ./blastdb/test_reference.fa -query for_blast.fa -outfmt "6 std qseq sseq slen" -word_size 11 -num_threads 3 -out blast_out.txt
../bin/getflanking.py polymarker_input_example.csv blast_out.txt temp_range.txt
cat temp_range.txt 
gawk '{ print $2,$3,$4 > "temp_marker_"$1".txt" }' temp_range.txt
for i in temp_marker*; do blastdbcmd -entry_batch $i -db ./blastdb/test_reference.fa > flanking_$i.fa; done

../bin/getkasp3.py 63 25  0
../bin/getCAPS.py 200 63 25 0

# command to make a blastable reference
makeblastdb -in test_reference.fa -parse_seqids -dbtype nucl


