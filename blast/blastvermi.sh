for i in ../../seqSamples/seq*fasta; do blastn -query $i -db /SAN/Eimeria/vermiformis/vermiformis_final.fa -outfmt "10 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid sseq qseq qframe sframe" -evalue 1e-5 -out `basename $i`.blastn; done

for i in ../data/seqApi/seq*fasta; 
do echo `basename $i` ;
blastn -query $i -db /SAN/Eimeria/vermiformis/vermiformis_final.fa -outfmt "10 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid sseq qseq qframe sframe" -evalue 1e-5 -out Ever/`basename $i`.blastn; done

for i in ../data/seqNucl/seq*fasta; 
do echo `basename $i` ;
blastn -query $i -db /SAN/Eimeria/vermiformis/vermiformis_final.fa -outfmt "10 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid sseq qseq qframe sframe" -evalue 1e-5 -out Ever/`basename $i`.blastn; done

