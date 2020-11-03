# silver-carp-genome
Set of scripts for processing primary analysis including the data filter,genome assemlby, genome annotation (repeat and gene annoation) and evolution analysis.

### S0.WGS.filter.pl
Perl script of filtering reads in fastq files of Illumina HiSeq4000 sequencing.
<pre>
Usage:perl S0.WGS.filter.pl <lane.lst> <lib.lst> [maxjob, default 10]
        -t <int>    thread number for each filter_data job,default 8
	-m <int>    the reads pair number in buffer,default 2000000
	-Q <int>    quality system: 64 or 33 (default: 33)
	-p <str>    set the qsub -P,default=rdtest
        -q <str>    set the qsub -q,default=bc.q
	-h     output help information to screen
</pre>

### S0.WGS.filter.sh
</br>Shell script for processing WGS.filter.pl. 
<pre>
perl ./S0.WGS.filter.pl lane.list lib.list
</pre>
The lane.list file was formatted like this:
<pre>
xx_1.fq 10 10 40
xx_2.fq 10 10 40
</br>read_x_file read_front_cutoff_bases_number read_tail_cutoff_bases_number low_quality_bases_cutoff
</pre>

</br>The lib.lst file was formatted like this:
<pre>
DFAADTAAPEI-55 500
library_ID  insert_size
</pre>

### S1.Assembly.Soapdenovo.sh
<pre>
./SOAPdenovo-63mer pregraph -s lib.cfg -K 33 -p 40 -d 1 -o bighead > pregraph.log
./SOAPdenovo-63mer contig -g bighead -R -p 40 > contig.log
./SOAPdenovo-63mer map -s lib.cfg -p 30 -g bighead -k 33 > map.log
./SOAPdenovo-63mer scaff -g bighead -F -p 40 > scaff.log
</pre>
</br>The script includes four steps and main parameters of Soapdenovo assembly.


### S2.RepeatMasker.find.sh and S2.RepeatProteinMask.find.sh
<pre>
./RepeatMasker -nolow -no_is -norna -engine ncbi  -parallel 1 -lib ./RepBase21.01/RepeatMaskerLib.embl.lib bighead.genome.fa
./RepeatProteinMask -engine ncbi -noLowSimple -pvalue 0.0001 bighead.genome.fa
</pre>
</br>We used the known repeat database to find repeats using corresponding software RepeatMasker and RepeatProteinMask, identifying TEs at the DNA and protein level, respectively.

### S3.Augustus.gene.sh
<pre>
/augustus-3.2.1/bin/augustus --species=zebrafish --AUGUSTUS_CONFIG_PATH=./config/ --uniqueGeneId=true --noInFrameStop=true --gff3=on --strand=both bighead.genome.rmask.fa > bighead.genome.rmask.fa.augustus
</pre>
</br>With the help of HMM model, de novo prediction was performed based on the repeat-masked genome. Programs we applied were AUGUSTUS.

### S3.Genewise.homo.gene.sh
<pre>
./genewise/wise2.4.1/src/bin/genewise -trev -sum  -genesf  -gff ref.gene.gff bighead.genome.fa >  bighead.genome.fa.genewise
</pre>
</br>Homologous proteins of other species (e.g. from Ensembl) were mapped to the genome using TblastN with an E-value cutoff 1e-5, the aligned sequences as well as their corresponding query proteins were then filtered and passed to GeneWise for searching accurate spliced alignments.

### S4.phylo_tree.sh
<pre>
./S4.phylo_tree.pl single-copy.pep.phy --type ml -d aa -outdir result/
</pre>
</br>Single-copy gene families were used to reconstruct the phylogentic tree. 4-fold degenerate sites or amino acid sites were extracted from each single-copy family and concatenated to one super gene for each species. Then, PhyML or Mrbayes is used to construct the phylogenetic tree.

### S4.phylo_tree.pl
<pre>

Name
     phylo_tree.pl - Using TreeBeST or PhyML or MrBayes to 
     reconstruct phylogenetic tree
Usage
     perl phylo_tree.pl <infile> [options]
      <infile>         input phylip format(sequential) or mfa format sequence file;
      -format <str>    file formats: phylip(default); mfa;
      -type <str>      method: nj(default); ml; bayes;
      -outdir <str>    save results to this directory, default ./;
      -clean           delete temporary files, default not;
      -help            show help information;
 
     TreeBeST's relevant parameters(nj):
      -d <str>         data type: nt, nucleotide(default); aa, amino acid;
      -t <str>         model for codon & nt: ntmm; dn; ds; dm(default); 
                       model for aa: mm(default); jtt; kimura;
      -b <num>         bootstrap times, default 100; 

     PhyML's relevant parameters(ml):
      -d <str>         data type: nt, nucleotide(default); aa, amino acid;
      -m <str>         models for nt: HKY85(default); JC69; K80; F81; F84; TN93; GTR;
                       models for aa: WAG(default); JTT; MtREV; Dayhoff;
      -b <num>         set bootstrap(int>0) or aLRT method(-1,-2,-4), default -4;
      -rates <str>     rates across sites: gamma(default); invgamma; equal; propinv;

     MrBayes's relevant parameters(bayes):
      -data <str>      data type: dna(default); protein;
      -nst <num>       nucleotide models: 6, GTR(default); 2, HKY; 1, JC;
      -aamodel <str>   amino-acid models: Poisson; Jones(default); Dayhoff; WAG; BLOSUM;
      -rates <str>     rates across sites: gamma(default); invgamma; equal; propinv;
      -ngen <num>      generation number, default 100000;
      -samplefreq <num>    sample one time for every 100 (default) generations;
      -mpirun <num>    run MrBayes with mpirun using multi CPUs or not, default 1;

Example
     perl phylo_tree.pl <in.phy> -type nj -d nt
     perl phylo_tree.pl <in.phy> -type ml -d nt
     perl phylo_tree.pl <in.phy> -type bayes -data dna
</pre>
