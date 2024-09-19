# prepare COMBINED E coli PPI network data and preprocess it

#download this e. coli PPI datasets from https://data.mendeley.com/datasets/3rj8rsbmhd/1
#which were used in the paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9249336/


#use only Nature Biotech. paper's PPI network from the 3 networks available in that zip file
cat e-coli-data/Nat\ Biotechnol.\ 2014\ Table\ S5-geneName.txt | tr '\t' ' ' | grep -v protein > ecoli-all.sif

#prepared e-coli-ppi-cong-2019.edgelist from https://www.science.org/doi/suppl/10.1126/science.aaw6718/suppl_file/aaw6718_tables_s1-16.xlsx

# append edges from e-coli-ppi-cong-2019.edgelist to ecoli-all.sif
cat e-coli-data/e-coli-ppi-cong-2019.edgelist >> ecoli-all.sif

#prepare a graph without self-loop and multi-edges
python ../Graph_check.py --i ecoli-all.sif --o ecoli-all.pel

# prepare an edgelist file (ecoli-all-transnum.pel) where protein names are renamed by their indexes 
# and the mapping from protein name to their indexes are saved in ecoli-all.map
../transnum.pl ecoli-all.map < ecoli-all.pel > ecoli-all-transnum.pel

# get an adjacency list representation of the above graph in ecoli-all.grh
../transgrh.pl < ecoli-all-transnum.pel > ecoli-all.grh

#downloaded protein E coli protein complexes from ECOCYC: https://ecocyc.org/group?id=:ALL-PROTEINS-2&orgid=ECOLI
#and saved that as All-protein-complexes-of-E.-coli-K-12-substr.-MG1655.txt

# save list of protein complexes (each line is a space-separated list of proteins in a complex) in ecocyc-complexes.txt
cut -f2 e-coli-data/All-protein-complexes-of-E.-coli-K-12-substr.-MG1655.txt | sed 's| //||g' | grep -v Gene | sed '/^$/d' > ecocyc-complexes.txt

perl -lane 'print "C$.\t$_"' ecocyc-complexes.txt > ecocyc-complexes-w-dummy-names.txt

# select all protein complexes with 5 or more proteins and save them in ecoli-protein-complexes-gte-5-w-names-densities.txt
# along with each of their sizes, #of nodes+edges in its induced subgraph, and that subgraph's density
k=5
python subgraph_densities.py --g ecoli-all.pel --s ecocyc-complexes-w-dummy-names.txt --m $k > ecoli-protein-complexes-gte-"$k"-w-names-densities.txt

# for each density cutoff, prepare a file containing protein complexes having at least that much density in its induced subgraph
for i in $(seq 0.9 -0.1 0.6); do perl -slanF'\t' -e 'print "$F[0]\t$F[1]" if $F[-1] >= "$i"' -- -i="$i" ecoli-protein-complexes-gte-"$k"-w-names-densities.txt > ecoli-complexes-w-names-den-gte-"$i".txt; done


## ============================= COMPUTING pseudo-cliques via ClusterONE and FPCE for each density cutoff in {0.6,...0.9} ====================================
for i in $(seq 0.9 -0.1 0.6); do echo "Cluster1 Run for d = $i"; /usr/bin/time -o rss-time-cl1-"$k"-"$i".txt -f "RSS=%M TIME=%S+%U" java -jar ../code/cluster_one-1.0.jar -d "$i" -s "$k" --no-merge ecoli-all.pel | sed 's/\t/ /g' > cl1-ecoli-"$k"-"$i".out; done > cl1-"$k"-"$i".log

for i in $(seq 0.9 -0.1 0.6); do echo "FPCE Run for d = $i"; /usr/bin/time -o rss-time-fpce-"$k"-"$i".txt -f "RSS=%M TIME=%S+%U" ../code/FPCE/fpce M -l "$k" ecoli-all.grh "$i" fpce-ecoli-"$k"-"$i".out; ../untransnum.pl ecoli-all.map < fpce-ecoli-"$k"-"$i".out > fpce-ecoli-"$k"-"$i"-genenames.out; done > fpce-"$k"-"$i".log



########################## FastQC result comparison with FPCE #######################
../transfastqc.pl < ecoli-all-transnum.pel > ecoli-all.fastqc
for i in $(seq 0.9 -0.1 0.6); do echo "FastQC Run for d = $i"; /usr/bin/time -o rss-time-fastqc-"$k"-"$i".txt -f "RSS=%M TIME=%S+%U" ../code/SIGMOD24-MQCE-main/code/FastQC/build/FastQC -f ecoli-all.fastqc -u 5 -g $i -q 1; sed -i 's/[^ ]* //' output; ../untransnum.pl ecoli-all.map < output > fastqc-ecoli-5-$i-genenames.out; done > fastqc-"$k"-"$i".log

################# enrichment analysis ######################
Download ecocyc.gaf from https://current.geneontology.org/products/pages/downloads.html
grep -v ^! ecocyc.gaf | cut -f3,5 | sort -k1 > ecocyc-gene2go.txt
perl -lanF'\t' -e '$h{$F[0]} .= "$F[1];"; END {print "$_\t$h{$_}" for keys %h}' ecocyc-gene2go.txt > ecocyc-gene2golist.txt; sed -i 's/;$//g' ecocyc-gene2golist.txt

cd goatools-main

#to get enrichment results of FPCE-clusters with 0.8 density cutoff
i=0.8; ln=1; cat ../fpce-ecoli-5-"$i"-genenames.out | while read line; do echo $line | tr ' ' '\n' > tmp.study; python scripts/find_enrichment.py tmp.study ecoli-all.pop ecocyc-gene2golist.txt --outfile=goea-fpce-ecoli-5-"$i"-comm-"$ln".tsv --pval_field=fdr_bh --ns=CC --ev_exc=IEA --pval=0.05; ((ln+=1)); done

#to get enrichment results of FPCE-clusters with 0.8 density cutoff
i=0.8; ln=1; cat ../fpce-ecoli-5-"$i"-genenames.out | while read line; do echo $line | tr ' ' '\n' > tmp.study; python scripts/find_enrichment.py tmp.study ecoli-all.pop ecocyc-gene2golist.txt --outfile=goea-fpce-ecoli-5-"$i"-comm-"$ln".tsv --pval_field=fdr_bh --ns=CC --ev_exc=IEA --pval=0.05; ((ln+=1)); done


#compare enrichment results
wc -l ../cl1-ecoli-5-"$i".out
ls -1 goea-cl1-ecoli-5-"$i"-comm-*.tsv | wc -l

wc -l ../fpce-ecoli-5-"$i".out
ls -1 goea-fpce-ecoli-5-"$i"-comm-*.tsv | wc -l

