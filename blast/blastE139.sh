#!/bin/bash
MYDB="/SAN/Eimeria_Wild_Genomes/Assemblies_Alice_final/E139_final.genome.scf.fasta"
OUTFOLDER="/SAN/Alices_sandpit/AA_mouse1/git/AA_mouse1/blast/E139/"

# do it once and for all
makeblastdb -in $MYDB -dbtype nucl

for i in /SAN/Alices_sandpit/AA_mouse1/data/seqApi/seq*fasta; 
do echo `basename $i` ;
blastn -query $i -db $MYDB -outfmt "10 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid sseq qseq qframe sframe" -evalue 1e-5 -out $OUTFOLDER/`basename $i`.blastn; done

for i in /SAN/Alices_sandpit/AA_mouse1/data/seqNucl/seq*fasta; 
do echo `basename $i` ;
blastn -query $i -db $MYDB -outfmt "10 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid sseq qseq qframe sframe" -evalue 1e-5 -out $OUTFOLDER/`basename $i`.blastn; done
