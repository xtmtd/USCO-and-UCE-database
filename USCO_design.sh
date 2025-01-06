#USCO_design

mkdir design && cd design
cp -rf ../proteins raw
cd raw/proteins/
for file in * ; do echo ${file%.*} >> ../species.list; done

#assessment protein completeness
for species in $(cat species.list); do run_BUSCO.py -i proteins/"$species".faa -m prot -c 28 -o "$species" -l ~/install/busco-3.1.0/lineage/arthropoda_odb10; mv run*/short* .; rm -rf run* tmp; done
cp raw/species.list .


#revise the sequence names for each species
mkdir proteins-raw
for species in $(cat species.list); do  cat proteins/$species.faa | awk '/^>/{print ">`$species`_" ++i; next}{print}' | seqkit replace -p .+ -r "`$species`_{nr}" --nr-width 5 | seqkit seq -u -w 0 | sed "s/*//g" > proteins-raw/$species.fa; sed -i "s/>/>"$species"/g" proteins-raw/$species.fa; done
mkdir 0-proteins
for species in $(cat species.list); do cat proteins-raw/"$species".fa | seqkit fx2tab | grep -v "\." | seqkit tab2fx -w 0 > ../0-proteins/"$species".fa; done


#search gene families
~/install/OrthoFinder-2.5.4/orthofinder.py -f 0-proteins/ -t 30 -S diamond_ultra_sens 
#check the species tree constructed by the defaults, if incorrect, provide the correct ones and re-start the analyses
#orthofinder -t 16 -s species_tree.tre -ft 0-proteins/OrthoFinder/Results_*/
mkdir 1-orthofinder && cd 1-orthofinder
cp -rf ../0-proteins/OrthoFinder/Results_* ./


#get HOGs of N0 node
python ~/install/OrthoFinder-2.5.4/tools/create_files_for_hogs.py Results_* Results_* N0


#get gene count of N0
cd Results_*/N0
cp ../Phylogenetic_Hierarchical_Orthogroups/N0.tsv .
python ~/install/OrthoFinder-2.5.4/tools/orthogroup_gene_count.py N0.tsv


#summarize Orthofinder HOG results
cat N0.GeneCount.csv | dos2unix | sed "s/^N0.//g" | awk -v FS="\t" -v OFS="\t" '{if(NR==1){print $0, "Total"} else{for(i=4;i<=NF;i++) sum+=$i;print $0,sum}sum=0}' > N0.GeneCount1.csv
awk 'NR!=1{for(i=4;i<=(NF-1);i++) if($i==$NF){print $1}}' N0.GeneCount1.csv > list.orthogroups.species-specific #species-specific orthogroups
grep -v -P "\t"0"\t" N0.GeneCount1.csv | awk 'NR!=1{print $1}' > list.universal_orthogroups #orthogroups with all species present
no_species=$(awk -v FS="\t" '/^HOG/{print NF-3}' N0.GeneCount.csv) #Number of species
grep -f list.universal_orthogroups N0.GeneCount1.csv | awk '$NF=='$no_species'{print $1}' > list.single-copy #single-copy orthogroups

no_genes=$(cat ../raw/proteins/*  | grep -c "^>") #Number of genes
no_genes_orthogroup=$(awk 'NR!=1{sum+=$NF} END{print sum}' N0.GeneCount1.csv) #Number of genes in orthogroups
no_orthogroup=$(grep -c "HOG0" N0.GeneCount.csv)
no_genes_species_specific=$(grep -f list.orthogroups.species-specific N0.GeneCount1.csv | awk '{sum+=$NF} END{print sum}')
no_genes_universal=$(cat N0.GeneCount1.csv | grep -f list.universal_orthogroups | awk 'NR!=1{sum+=$NF} END{print sum}')

echo -e "Number of species\t${no_species}" >> Statistics_Overall.tsv
echo -e "Number of genes\t${no_genes}" >> Statistics_Overall.tsv
echo -e "Number of genes in orthogroups\t${no_genes_orthogroup}" >> Statistics_Overall.tsv
echo "Number of unassigned genes" | awk '{print $0,"\t",'$no_genes'-'$no_genes_orthogroup'}' >> Statistics_Overall.tsv
echo "Percentage of genes in orthogroups" | awk '{print $0,"\t",'$no_genes_orthogroup'/'$no_genes'*100"%"}' >> Statistics_Overall.tsv
echo "Percentage of unassigned genes" | awk '{print $0,"\t",100-'$no_genes_orthogroup'/'$no_genes'*100"%"}' >> Statistics_Overall.tsv
echo -e "Number of orthogroups\t${no_orthogroup}" >> Statistics_Overall.tsv
echo "Mean orthogroup size" | awk '{print $0,"\t",'$no_genes_orthogroup'/'$no_orthogroup'}' >> Statistics_Overall.tsv
echo -e "Number of species-specific orthogroups\t$(cat list.orthogroups.species-specific | wc -l)" >> Statistics_Overall.tsv
echo -e "Number of genes in species-specific orthogroups\t${no_genes_species_specific}" >> Statistics_Overall.tsv
echo "Percentage of genes in species-specific orthogroups" | awk '{print $0,"\t",'$no_genes_species_specific'/'$no_genes'*100"%"}' >> Statistics_Overall.tsv
echo -e "Number of orthogroups with all species present\t$(cat list.universal_orthogroups | wc -l)" >> Statistics_Overall.tsv
echo -e "Number of genes in orthogroups with all species present\t${no_genes_universal}" >> Statistics_Overall.tsv
echo "Percentage of genes in orthogroups with all species present" | awk '{print $0,"\t",'$no_genes_universal'/'$no_genes'*100"%"}' >> Statistics_Overall.tsv
echo -e "Number of single-copy orthogroups\t$(cat list.single-copy | wc -l)" >> Statistics_Overall.tsv


#summary of each species
mkdir temp
echo -e """\n""Number of genes""\n""Number of genes in orthogroups""\n""Percentage of genes in orthogroups (%)""\n""Number of unassigned genes""\n""Percentage of unassigned genes (%)""\n""Number of orthogroups containing species""\n""Percentage of orthogroups containing species (%)""\n""Number of species-specific orthogroups""\n""Number of genes in species-specific orthogroups""\n""Percentage of genes in species-specific orthogroups (%)" > temp/head

no_species=$(awk -v FS="\t" '/^HOG/{print NF-3}' N0.GeneCount.csv) #Number of species
no_orthogroup=$(grep -c "HOG0" N0.GeneCount.csv) #number of orthogroups

for column in $(seq ${no_species})
  do
    species=$(head -n 1 N0.GeneCount.csv | cut -f 4- | cut -f $column | col -b)
    no_genes=$(grep -c ">" ../raw/proteins/${species}*)
    no_genes_orthogroup=$(cut -f 4- N0.GeneCount.csv | cut -f $column | awk 'NR!=1{sum+=$1} END{print sum}')
    pct_genes_orthogroup=$(echo | awk '{print '$no_genes_orthogroup'/'$no_genes'*100}')
    no_genes_unassign=$(echo | awk '{print '$no_genes'-'$no_genes_orthogroup'}')
    pct_genes_unassign=$(echo | awk '{print '$no_genes_unassign'/'$no_genes'*100}')
    no_orthogroup_species=$(sed 1d N0.tsv | grep -c $species)
    pct_orthogroup_species=$(echo | awk '{print '$no_orthogroup_species'/'$no_orthogroup*100'}')
    grep -f list.orthogroups.species-specific N0.tsv | awk '/'$species'/{print $1}' | sed "s/^N0.//g" > list.t1
    no_orthogroups_species_specific=$(cat list.t1 | wc -l)
    no_genes_species_specific=$(grep -f list.t1 N0.GeneCount.csv | cut -f 4- | cut -f $column | awk '{sum+=$1} END{print sum}')
    pct_genes_species_specific=$(echo | awk '{print '$no_genes_species_specific'/'$no_genes'*100}')
    echo -e $species"\n"$no_genes"\n"$no_genes_orthogroup"\n"$pct_genes_orthogroup"\n"$no_genes_unassign"\n"$pct_genes_unassign"\n"$no_orthogroup_species"\n"$pct_orthogroup_species"\n"$no_orthogroups_species_specific"\n"$no_genes_species_specific"\n"$pct_genes_species_specific > temp/${species}.txt
done

paste temp/head temp/*txt > Statistics_PerSpecies.tsv
rm -rf temp list.t1 


##at least 13 species with single-copy orthologs
cat Results_Sep14/N0/N0.tsv | sed 1d | cut -f1 > list.orthogroups

for id in $(cat list.orthogroups);   do     num=$(cat N0.GeneCount.csv | grep -P ^"$id" | cut -f1-18 | csvtk -t transpose | grep -c "^1$");     test "$num" -ge 13 && echo "$id" >> list.candidates; done

#extract HOGs in candidates.list
mkdir -p candidates_seqs
#revise the sequence name for each orthogroup sequence
for id in $(cat list.candidates); do cp Results_*/N0/HOG_Sequences/"$id".fa candidates_seqs/; done
sed -i 'N;s/_[0-9][0-9][0-9][0-9][0-9]//;' candidates_seqs/*fa


#delete the duplicated sequences for each HOG
for loci in $(cat list.candidates)
  do 
	  for species in $(cat ../species.list)
	    do
		    num=$(cat candidates_seqs/"$loci".fa | grep -c "$species")
		    if [ "$num" -ge 2 ]; then
			echo "$species" >> temp.list
			cat candidates_seqs/"$loci".fa | seqkit grep -v -f temp.list > $loci.fa
			mv $loci.fa candidates_seqs/
			rm temp.list
		    fi
	    done
  done	  


#Rename files in the candidates_seqs folder.
rename 's/^N0\.//' *
sed -i 's/^N0\.//' list.candidates


#check the possible errors
seqkit stat -T *fa > ../stat.candidates0
#delete genes of length shorter than 100 (candidates.list1)
cat stat.candidates0 | sort -k7n > stat.candidates1   #mean length of all loci >=100

grep -E '^\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+[0-9]+\s*$' stat.candidates1 | awk '$7<100{print $1}' > t1
for filename in $(cat t1); do rm "candidates_seqs/$filename"; done
rm t1
cd candidates_seqs/
for file in * ; do echo ${file%.*} >> ../list.candidates1; done

#delete the sequence of length < 1/2 or 2 times of mean length for each OG
mkdir temp && cd temp

for id in $(cat ../list.candidates1)
  do
    mean=$(cat ../stat.candidates0 | grep "$id" | cut -f7)
    cat ../candidates_seqs/"$id".fa | seqkit fx2tab -n -l | sed "s/\t\t\t/\t/g" > "$id".length
    min=$(echo "scale=0;$mean/2"|bc)
    max=$(echo "scale=0;$mean*2/1"|bc)
    line=$(cat "$id".length | wc -l)
    for LINE in $(seq $line)
      do
        num=$(cat "$id".length | sed -n "$LINE"p | cut -f2)
	touch "$id".filter.list
        test "$num" -le "$min" -o "$num" -ge "$max" && cat "$id".length | sed -n "$LINE"p | cut -f1 >> "$id".filter.list
      done
    test -s "$id".filter.list && cat ../candidates_seqs/"$id".fa | seqkit grep -v -f "$id".filter.list > "$id".fa
    rm "$id".length "$id".filter.list
  done

mv *fa ../candidates_seqs/

seqkit stat -T *.fa > ../stat.candidates2

#filter loci with sequence number <8
cat stat.candidates2 | sort -k4n | sed -n 2,8p | cut -d "." -f1 > list.filter
for id in $(cat list.filter); do rm candidates_seqs/"$id".fa; done
cat list.candidates1 | grep -v -f list.filter > list.candidates2
# 4360 loci remaning

#align sequences for each gene with custum script using MAGUS
cd 2-align
bash align_MAFFT.sh
#filter loci with sequence number <13
#seqkit stat -T *.fa > ../stat.candidates3
#cat stat.candidates3 | sort -k4n | sed -n 2,288p | cut -d "." -f1 > list.filter2
#for id in $(cat list.filter2); do rm ../2-align/"$id".fa; done
#cat list.candidates2 | grep -v -f list.filter2 > list.candidates3


#Using taper, treeshrink detects possible errors in the sequence
cd 3-Error_detection
bash loci_filtering_alignment-based.sh
bash gene_trees.sh  #Select LG model
bash treeshrink.sh   #The fault tolerance rates were selected as 0.05 and 0.1

##build amino acid-level hidden Markov model (HMM) profiles
cd 2-align
for file in * ; do echo ${file%.*} >> ../align.list; done 

mkdir 4-hmms && cd 4-hmms
for gene in $(cat ../align.list); do cp ../2-align/$gene.fa .; done
cat ../align.list | parallel -j 28 hmmbuild {}.hmm {}.fa
#for gene in $(sed -n '1,1000p' ../1-orthofinder/candidates.list1); do hmmbuild $gene.hmm $gene.fa; done

#check errors
#for gene in $(cat ../align.list); do test -s $gene.hmm && echo || echo "$gene" >> ../hmm.fail.list; done

#search BUSCO sequences against the complete library of HMM profiles
mkdir 5-hmmsearch && cd 5-hmmsearch/
mkdir output
temp_fun() {   hmmsearch ../4-hmms/"$1".hmm ../1-orthofinder/candidates_seqs/"$1".fa > output/"$1".hmmsearch.out;   }
export -f temp_fun
cat ../align.list | parallel -I% -j 28 --max-args 1 temp_fun %

#Check if the HMM search output file is empty, and if it is empty, append the gene name to "../hmmsearch.fail.list"  file
for gene in $(cat ../1-orthofinder/candidates.list1); do test -s output/$gene.hmmsearch.out && echo || echo "$gene" >> ../hmmsearch.fail.list; done

mkdir statistics
for gene in $(cat ../align.list); do cat output/$gene.hmmsearch.out | sed -n '15,29p' > statistics/$gene.hmmsearch.stats; done

#extract the bitscore for each HOG
mkdir bit_score
for gene in $(cat ../align.list); do cat statistics/$gene.hmmsearch.stats | awk '{print $2}' | egrep '^[0-9]' > bit_score/$gene.bit_score; done
for gene in $(cat ../align.list); do num1=$(cat bit_score/"$gene".bit_score | sort -k1n | head -n 1); num2=$(echo "scale=2;$num1*0.9"|bc); echo $num2 > lowest_bit_score/$gene.lowest_bit_score; done
for gene in $(cat ../align.list); do sed -i "s/^/"$gene"\t/g" lowest_bit_score/$gene.lowest_bit_score; done
cat lowest_bit_score/* > scores_cutoff


#set lengths_cutoff
mkdir 6-lengths && cd 6-lengths
for id in $(cat ../align.list)
  do
    cat ../1-orthofinder/candidates_seqs/"$id".fa | seqkit fx2tab -n -l | cut -f4 > length."$id"
    mean=$(cat length."$id" | awk '{sum+=$1} END {print sum/NR}' | awk '{printf "%.0f", $1}')
    sigma=$(awk '{x[NR]=$0; s+=$0; n++} END{a=s/n; for (i in x){ss += (x[i]-a)^2} sd = sqrt(ss/n); print sd}' length."$id")
    diff=$(echo "scale=4;($sigma-5)"|bc)
    num1=`echo "$diff >= 0" |bc`
    test "$num1" = 0 && sigma=5  #when sigma<5, restrict it as 5
    echo -e "$id""\t"0"\t""$sigma""\t""$mean" >> lengths_cutoff
    sigma2=$(echo "scale=15;($sigma*2)"|bc)
    echo -e "$id""\t"0"\t""$sigma2""\t""$mean" >> lengths_cutoff_2sigma #calculate 2σ
    rm length."$id"
  done
#keep the 15 decimal places for sigma and integer for mean length 


#generate consensus sequences
mkdir 7-consensus && cd 7-consensus/
for gene in $(cat ../align.list); do hmmemit ../4-hmms/$gene.hmm > $gene.ancestral; done
cat * | sed "s/-sample1//g" > ancestral

#generate variant sequences
for gene in $(cat ../align.list); do hmmemit -N 10 ../4-hmms/$gene.hmm > $gene.ancestral_variants; done
cat *variants | sed "s/-sample/_/g" > ancestral_variants
rm *.ancestral *.ancestral_variants


#generate block profile for each OG
mkdir 8-prfl && cd 8-prfl
for gene in $(cat ../align.list); do msa2prfl.pl ../2-align/$gene.fa > $gene.prfl; done


#annotate the buscos
mkdir 9-ipr && cd 9-ipr/
cp ../7-consensus/ancestral .
 ~/install/interproscan-5.53-87.0/interproscan.sh -dp -b iprscan -f TSV -goterms -iprlookup -pa -t p -cpu 30 -i ancestral -appl Pfam,Smart,Superfamily,CDD


#generate Refseq.faa file required for busco version 5
mkdir 10-refseq && cd 10-refseq
&&cat ../collembola_odb2/info/species.info | sed "s/ /_/g;s/^/_/g;s/\t/_0\t/g" > map.txt

for gene in $(cat ../align.list)
  do
    cat ../1-orthofinder/candidates_seqs/"$gene".fa | seqkit fx2tab | cut -f1 > "$gene".t1
    cat ../1-orthofinder/candidates_seqs/"$gene".fa | seqkit fx2tab | cut -f2 > "$gene".t2
    for id in $(cat "$gene".t1)
      do
	taxid=$(cat map.txt | grep "$id" | cut -f1)
	echo -e "$gene""$taxid" >> "$gene".t3
	paste "$gene".t3 "$gene".t2 | seqkit tab2fx -w 0 > loci/"$gene".fa
      done
    rm "$gene".*
  done

cat loci/*fa > refseq_db.faa


awk 'match($0, /[1-9][0-9]*/) {print $1 "\t" substr($0, RSTART, RLENGTH) "at30001\t" $0}' ../align.list | awk '{print $1 "\t" $2}' > ogID.map.txt

#Rename the sequence in the "ancestral" file
while read line; do   key=$(echo $line | awk '{print $1}');   value=$(echo $line | awk '{print $2}');   sed -i "s/^>${key}\$/>${value}/" ancestral; done < info/ogID.map.

#Rename the sequence in the "ancestral_variants" file
awk '{ for (i=1; i<=10; i++) print $2 }' info/ogID.map.txt > t1
seqkit seq -n ancestral_variants | sed 's/.*_//' > t2
paste -d "_" t1 t2 > t3
cat ancestral_variants | seqkit seq -s -w 0 | sed "s/\.//g" > t4
paste t3 t4 | seqkit tab2fx | seqkit seq -w 0 > t5   ##(t5 is the final result)
rm t1 t2 t3 t4

for file in *; do new_prefix=$(echo "$file" | awk 'match($0, /[1-9][0-9]*/) {print substr($0, RSTART, RLENGTH)}'); mv "$file" "${new_prefix}at30001.prfl"; done
for file in *; do new_prefix=$(echo "$file" | awk 'match($0, /[1-9][0-9]*/) {print substr($0, RSTART, RLENGTH)}'); mv "$file" "${new_prefix}at30001.hmm"; done
for file in *; do new_prefix=$(echo "$file" | awk -F. '{print $1}'); sed -i "2s/\S\+/${new_prefix}/2" "$file"; done


#Rename the sequence in the "refseq_db.faa" file
seqkit seq -n refseq_db.faa | sed 's/_.*//' > t1
sed -E 's/^[^1-9]*([1-9][0-9]*)\b.*$/\1/' t1 | sed 's/$/at30001/' > t2
seqkit seq -n refseq_db.faa | sed 's/^[^_]*_//' > t3
paste -d "_" t2 t3 > t4
awk '{printf "%s:%06d\n", $0, NR}' t4 > t5
cat refseq_db.faa | seqkit seq -s -w 0 | sed "s/\.//g" > t6
paste t5 t6 | seqkit tab2fx | seqkit seq -w 0 > t7   ##(t7 is the final result)
 rm t1 t2 t3 t4 t5 t6

#Rename the contents of the "scores_cutoff, lengths_cutoff_2sigma, lengths_cutoff" file
paste <(cut -f 2 info/ogID.map.txt) <(cut -f 2- lengths_cutoff) > lengths_cutoff-1 #
paste <(cut -f 2 info/ogID.map.txt) <(cut -f 2- lengths_cutoff_2sigma) > lengths_cutoff_2sigma-1 #
paste <(cut -f 2 info/ogID.map.txt) <(cut -f 2- scores_cutoff) > scores_cutoff-1 #


awk '{print $1}' interproscan.annotations.info > t1
awk '{$1=""; print $0}' interproscan.annotations.info > t2
sed -E 's/^[^1-9]*([1-9][0-9]*)\b.*$/\1/' t1 | sed 's/$/at30001/' > t3
paste t3 t2 > t4


#Integrate "ancestral, ancestral_variants, dataset.cfg, hmms, info, lengths_cutoff, prfl, refseq_db.faa, scores_cutoff" file to "collembola_odb The 10" folder generates the data set.


#"Collembola_odb10" Data set BUSCO test
cd species/
ls > ../species.list
sed -i -E 's/\..*//' ../species.list
cd ..
for species in $(cat species.list); do busco -m genome -i "$species".fa -c $"THREADS" -o "$species" -l ~/install/busco-5.6.0/lineage/collembola_odb10 --offline -f; done