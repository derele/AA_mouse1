for i in ../data/seqApi/seq*fasta; 
do echo `basename $i` ;
blastn -query $i -db /SAN/db/blastdb/Eimeria_falciformis/Eimeria_contigs_final_AM.fa -outfmt "10 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid sseq qseq qframe sframe" -evalue 1e-5 -out Efal/`basename $i`.blastn; done

for i in ../data/seqNucl/seq*fasta; 
do echo `basename $i` ;
blastn -query $i -db /SAN/db/blastdb/Eimeria_falciformis/Eimeria_contigs_final_AM.fa -outfmt "10 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid sseq qseq qframe sframe" -evalue 1e-5 -out Efal/`basename $i`.blastn; done
