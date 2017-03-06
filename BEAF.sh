#!/usr/bin/env bash

ans="yes" #If you're sure you want to keep databases files from one loop to the next, change this to 'Y'. If you're sure you're not reusing .udb, change to 'N'

make_kp ()
{
if [[ "$CFLR" == "Y" ]]; then echo "******Reusing previous settings and files"; else
	cd $address
	rm -rf *.kp; rm -rf config.tmp; rm -rf config.file1
	cat config.file | awk NF > config.file1
	echo "# Checking if any buckets must be stored..."
	echo "X" > LastR1.kp
	sort -k3,3 config.file1 > doconfig.kp
	while read T1 T2 R1 R2 Ref SubRef Out; do
		LastR1=`cat LastR1.kp`	
		if [[ "$R1" == "$LastR1" ]]
		then
			echo "Y" > LastKeep.kp
		else
			echo "N" > LastKeep.kp
		fi
		cat Keep_config.kp LastKeep.kp > Keep_config2.kp
		mv Keep_config2.kp Keep_config.kp
		echo "$R1" > LastR1.kp
	done < doconfig.kp
	echo "N" > LastKeep.kp
	cat Keep_config.kp LastKeep.kp > Keep_config2.kp
	mv Keep_config2.kp Keep_config.kp
	sed -i -e 1,1d Keep_config.kp
	paste doconfig.kp Keep_config.kp > config.tmp
	rm -rf *.kp; rm -rf config.file1
	mv config.tmp config.kp
	cd $address; echo "1" > CR.step; CFLR="N"
fi
}

Check ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "1" ]]; then echo "******Jumping autocheking system"; else
	cd $address
	rm -rf *.check
	echo "# Checking your config file..."
	sort -k7,7 config.kp > config.check
	while read T1 T2 R1 R2 Ref SubRef Out Keep; do
		cd $address
		case $T1 in
			G|g|P|p|N|n)
			;;
			*)
				echo "# You're using a wrong config.file format. On the first column (T1), where you're currently using '$T1', use only G (for genome analysis), P (for protein analysis) or N (for protein nucleotide sequences analysis)"
				rm -rf config.kp config.check
				exit
			;;
		esac
		case $T2 in
			R|r)
				case $R2 in
					*.gz)
						if ! [ -s $R2 ];
						then
							echo "# Check your R2 file in '$R2'. The program either couldn't find the file or the file is empty."
							rm -rf config.kp config.check
							exit
						fi
					;;
					*)
						echo "# You're using a wrong config format. On the fourth column (R2), where you're currently using '$R2', you must use a gzipped file instead."
						rm -rf config.kp config.check
						exit
					;;
				esac
			;;
			I|i|F|f)
			;;
			*)
				echo "# You're using a wrong config.file On the second column (T2), where you're currently using '$T2', use only R (for paired end fastq files), I (for interleaved fastq file) or F (for interleaved fasta file)"
				rm -rf config.kp config.check
				exit
			;;
		esac
		case $R1 in
			*.gz)
				if ! [ -s $R1 ];
				then
					echo "# Check your R1 file in '$R1'. The program either couldn't find the file or the file is empty."
					rm -rf config.kp config.check
					exit
				fi
			;;
			*)
				echo "# You're using a wrong config.file On the third column (R1), where you're currently using '$R1', you must use a gzipped file instead."
				rm -rf config.kp config.check
				exit
			;;
		esac
		LastOut=`cat LastOut.check`
		if [[ "$Out" == "$LastOut" ]]
			then
				echo "# You're using a wrong config.file On your seventh column (Out), you've used the same name for your output folder more than once, repeating '$Out'"
				rm -rf config.kp config.check
				exit

		fi
		echo "$Out" > LastOut.check
	done < config.check
	rm -rf *.check
	cd $address; echo "2" > CR.step; CFLR="N"
fi
}

TimeHeader ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "2" ]]; then echo "******In this module, time is not measured"; else
	cd $address
	rm -rf header; rm -rf time.log2
	echo "_________________________________________________________________________________________________________
	Output|Sequence|Type1|Type2|Reference|Subref|Time|Reads|Buckets|ppm1|contigs|ORFs" > header
	cat time.log header > time.log2
	rm -rf header time.log
	mv time.log2 time.log
	cd $address; echo "3" > CR.step; CFLR="N"
fi
}

SoftTimeHeader ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "2" ]]; then echo "******In this module, time is not measured"; else
	cd $address
	rm -rf header; rm -rf time.log2
	echo "_________________________________________________________________________________________________________
	Output|Sequence|Type1|Type2|Reference|Subref|Time|Reads|Buckets|ppm1|contigs" > header
	cat time.log header > time.log2
	rm -rf header time.log
	mv time.log2 time.log
	cd $address; echo "3" > CR.step; CFLR="N"
fi
}

Trim ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "3" && `cat CR.step` != "24" ]]; then echo "******Skipping reads trimming opperations"; else
	cd $address
	rm -rf Buckets
	mkdir $address/Buckets
	cd $address/Buckets
	case $T2 in 
		R|r)
			echo "# Trimming your files..."
			echo "# Starting file 1..."
			cutadapt --minimum-length 80 --quality-base 24 --trim-n -o R1.trimmed.fastq.gz $R1 # Parameters of reads trimming should be specified here
			echo "# Starting file 2..."
			cutadapt --minimum-length 80 --quality-base 24 --trim-n -o R2.trimmed.fastq.gz $R2 # Parameters of reads trimming should be specified here
			echo "# Now we are merging files, it could take several minutes away..."
			zcat *.gz | gzip -c > FastaQ-zcat.gz
		;;
		I|i)
			echo "# Now we are trimming your files..."
			cutadapt --minimum-length 80 --quality-base 24 --trim-n -o FastaQ-zcat.gz $R1 # Parameters of reads trimming should be specified here
		;;
		F|f)
			gunzip -c <`ls *.gz`> FastaQ-zcat.fa
			rm -rf *.gz
		;;
	esac
	cd $address; echo "4" > CR.step; CFLR="N"
fi
}

QAnConversion ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "4" ]]; then echo "******Skipping FASTQ handling / FASTQ conversion"; else
	cd $address/Buckets
	case $T2 in
		R|r|I|i)
			cd $address/Buckets
			echo "# Starting quality assessment of trimming..."
			rm -rf FASTQCresults; rm -rf FastaQ-zcat.fa
			mkdir FASTQCresults
			fastqc -f fastq -o FASTQCresults FastaQ-zcat.gz
			cd FASTQCresults
			rm -rf *_fastqc
			cd $address/Buckets
			echo "# We will convert merged file to fasta format."
			gunzip -c FastaQ-zcat.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > FastaQ-zcat.fa
			echo "# Removing unwanted files..."
			rm -rf *.gz
		;;
		*)
			rm -rf FASTQCresults
			echo "# FASTA file type identified. Since FASTA does not have PHRED values, It wont be accessed by FASTQC algorithm."
		;;
	esac
	cd $address; echo "5" > CR.step; CFLR="N"
fi
}

CopyFile ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "3" && `cat CR.step` != "24" ]]; then echo "******Skipping files copying process"; else
	cd $address
	rm -rf Buckets
	mkdir Buckets
	cd $address/Buckets
	case $T2 in 
		R|r)
			echo "# Copying file 1 from storage..."
			cp -r $R1 $address/Buckets
			echo "# Copying file 2 from storage..."
			cp -r $R2 $address/Buckets
		;;
		I|i|F|f)
			echo "# Copying file from storage..."
			cp -r $R1 $address/Buckets
		;;
	esac
	cd $address; echo "4" > CR.step; CFLR="N"
fi
}

SoftMergeRename ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "4" ]]; then echo "******Skipping files merging / conversion process"; else
	cd $address/Buckets
	case $T2 in 
		R|r)
			rm -rf FastaQ-zcat.gz FastaQ-zcat.fa
			echo "Merging files"
			zcat *.gz | gzip -c > FastaQ-zcat.gz
			echo "# We will convert merged file to fasta format."
			gunzip -c FastaQ-zcat.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > FastaQ-zcat.fa
			rm -rf *.gz
		;;
		I|i)
			rm -rf FastaQ-zcat.gz FastaQ-zcat.fa
			gunzip -c `ls *.gz` | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > FastaQ-zcat.fa
		;;
		F|f)
			rm -rf FastaQ-zcat.fa
			gunzip -c <`ls *.gz`> FastaQ-zcat.fa
			rm -rf *.gz
		;;
	esac
	cd $address; echo "5" > CR.step; CFLR="N"
fi
}

BucketEngine ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "5" && `cat CR.step` != "24" ]]; then echo "******Skipping Buckets System"; else
	cd $address/Buckets
	if [ -s bk_list.txt ]
	then
		reads=`cat reads.nmb`
		for n in `cat bk_list.txt`; do
			head -9000000 FastaQ-zcat.fa > $n.bk
			sed -i -e 1,1d bk_list.txt
			sed -i '1,9000000d' FastaQ-zcat.fa		
		done
		if [ -s FastaQ-zcat.fa ];
		then
			mv FastaQ-zcat.fa last.bk
		fi
		rm -rf bk_list.txt
	else
		if [ -s buckets_list.txt ];
		then
			buckets=`ls *.bk | wc -l`
			reads=`cat reads.nmb`
		else
			if [ -s reads.nmb ]
			then
				touch reads.nmb
			else
				grep ">" FastaQ-zcat.fa | wc -l > reads.nmb
				reads=`cat reads.nmb`
				buckets=`expr $reads / 4500000`
			fi
			if [ "$buckets" -eq "0" ];
			then
				echo "# We identified $reads reads. It would take no buckets, avoiding this step."
				mv FastaQ-zcat.fa unique.bk
			else
				echo "# We identified $reads reads. It would take $buckets bucket steps."
				echo "### Starting operation of cutting and readapting..."
				echo "## Generating buckets..."
				shuf -i 0-1000000 -n $buckets > bk_list.txt
				echo "# Sorting buckets and setting up configurations..."
				for n in `cat bk_list.txt`; do
					head -9000000 FastaQ-zcat.fa > $n.bk
					sed -i -e 1,1d bk_list.txt		
					sed -i '1,9000000d' FastaQ-zcat.fa		
				done
				if [ -s FastaQ-zcat.fa ];
				then
					mv FastaQ-zcat.fa last.bk
				fi
				rm -rf bk_list.txt
			fi
		fi
	fi
	echo "### Removing temporary files (stage 1)..."
	rm -rf *.txt; rm -rf *.fa; rm -rf *.gz
	ls *.bk > buckets_list.txt
	cp -r buckets_list.txt buckets_search.txt
	rm -rf bk_list.txt
	buckets=`ls *.bk | wc -l`
	cd $address; echo "6" > CR.step; CFLR="N"
fi
}

Filter1 ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "6" ]]; then echo "******Skipping Filtering System 1"; else
	cd $address
	case $T1 in
		G|g)	
			echo "# Starting searches..."
			echo "### Making database from reference genome... "
			cd $address/Reference_seqs
			cp -r $Ref none
			cat none | awk NF > none1
			sed -i '/>/d' none1
			cat none1 | tr -d '\n' | sed 's/.\{100\}/&\n>\n/g' | sed '1s/.*/>\n&/' | awk -vRS=">" '{$0=n$0;ORS=RT}++n' > md8
			rm -rf none none1
			usearch -makeudb_usearch md8 -output $Ref.udb
			rm -rf md8
			mv $Ref.udb $address/Buckets
			cd $address/Buckets
			rm -rf *.m7
			case $Keep in
				Y|y)
					echo "# Searching against reference $Ref and keeping buckets..."
					for buck in `cat buckets_search.txt`; do
						if [ -s $buck.m8 ]
						then
							touch $buck.m8
						else
							usearch -usearch_global $buck -db $Ref.udb -strand both -id 0.95 -evalue 1e-20 -matched $buck.m7 # Parameters of reads search by Usearch algorithm should be specified here
							mv $buck.m7 $buck.m8
							sed -i -e 1,1d buckets_search.txt
						fi
					done
				;;
				*)
					echo "# Searching against reference $Ref and removing buckets..."
					for buck in `cat buckets_search.txt`; do
						usearch -usearch_global $buck -db $Ref.udb -strand both -id 0.95 -evalue 1e-20 -matched $buck.m7 # Parameters of reads search by Usearch algorithm should be specified here
						mv $buck.m7 $buck.m8
						rm -rf $buck
						touch $buck
						sed -i -e 1,1d buckets_search.txt
					done
					rm -rf *.bk
				;;
			esac
			rm -rf *.udb; rm -rf *.m7; rm -rf buckets_search.txt
		;;
		P|p|N|n)
			cd $address/Buckets
			if [ -s udblist ]
			then
				touch udblist
			else
				rm -rf *.udb
				cd $address/Reference_seqs
				rm -rf TestExtension
				mkdir TestExtension
				cp -r $Ref TestExtension
				cd $address/Reference_seqs/TestExtension
				ls *.fa > list.fa; ls *.fasta > list.fasta; ls *.fas > list.fas; ls *.faa > list.faa; ls *.fna > list.fna; ls *.fsa > list.fsa
				cat list.f* > list; rm -rf list.f*
				if [ -s list ]
				then
					echo "# Recognized $Ref file as fasta format. Making udb..."
					cat $Ref | awk NF > none
					sed -i '/>/d' none
					cat none | tr -d '\n' | sed 's/.\{100\}/&\n>\n/g' | sed '1s/.*/>\n&/' | awk -vRS=">" '{$0=n$0;ORS=RT}++n' > md8
					usearch -makeudb_usearch md8 -output $Ref.udb
					rm -rf md8 $Ref
					mv $Ref.udb $address/Buckets
					cd $address/Reference_seqs
					rm -rf TestExtension
					cd $address/Buckets
					ls *.udb > udblist
				else
					ls *.udb > listudb
					if [ -s listudb ]
					then
						echo "# Recognized $Ref file as .udb format."
						mv $Ref $address/Buckets
						cd $address/Reference_seqs
						rm -rf TestExtension
						cd $address/Buckets
						ls *.udb > udblist
					else
						echo "Couldn't recognize $Ref file as neither fasta nor udb format. Will try to use it as udb regardless."
						mv $Ref $address/Buckets
						cd $address/Reference_seqs
						rm -rf TestExtension
						cd $address/Buckets
						ls $Ref > udblist
					fi
				fi
			fi
			cd $address/Buckets
			dbinuse=`cat udblist`
			rm -rf *.m7
			case $Keep in
				Y|y)
					echo "# Searching against reference $Ref and keeping buckets..."
					for buck in `cat buckets_search.txt`; do
						usearch -usearch_local $buck -db $dbinuse -strand both -id 0.25 -evalue 1e-5 -matched $buck.m7
						mv $buck.m7 $buck.m8
						sed -i -e 1,1d buckets_search.txt
					done
				;;
				*)
					echo "# Searching against reference $Ref and removing buckets..."
					for buck in `cat buckets_search.txt`; do
						usearch -usearch_local $buck -db $dbinuse -strand both -id 0.25 -evalue 1e-5 -matched $buck.m7
						mv $buck.m7 $buck.m8
						rm -rf $buck
						touch $buck
						sed -i -e 1,1d buckets_search.txt
					done
					rm -rf *.bk
				;;
			esac
			rm -rf *.udb; rm -rf udblist; rm -rf hits; rm -rf *.m7; rm -rf buckets_search.txt
		;;
	esac
	cd $address; echo "7" > CR.step; CFLR="N"
fi
}

PreLogGen ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "7" ]]; then echo "******Skipping generation of pre-log"; else
	cd $address/Buckets
	cat *.m8 > hits
	echo "## Calculating statistics..."
	hits=`grep ">" hits | wc -l`
	if [[ -s $address/Buckets/reads.nmb ]]
	then
		reads=`cat $address/Buckets/reads.nmb`
	else
		if [[ -s $address/OUTPUT/$Out/reads.nmb ]]
		then
			reads=`cat $address/OUTPUT/$Out/reads.nmb`
		fi
	fi
	ppm1=`expr 1000000 \* $hits / $reads`
	echo -e "RESULTS:
	File from: $R1\t$R2
	Results: $Out
	Reference: $Ref
	Reads: $reads
	Buckets: $buckets
	Hits: $hits
	Portion in ppm: $ppm1" > Log.txt
	echo "$ppm1" > $address/OUTPUT/$Out/ppm1.nmb
	cp -r reads.nmb $address/OUTPUT/$Out
	mv Log.txt $address/OUTPUT/$Out
	if [[ -d FASTQCresults ]]
	then
		cp -r FASTQCresults $address/OUTPUT/$Out
	fi
	cntg="0"
	echo "$cntg" > $address/OUTPUT/$Out/cntg.nmb
	rm -rf *.txt; rm -rf *.m8
	cd $address; echo "8" > CR.step; CFLR="N"
fi
}

G_Prepare_SPADES ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "8" ]]; then echo "******Skipping Assembly Preparation Module"; else
	cd $address/Buckets
	echo "# Making contigs for $Out..."
	cp hits hits.fasta
	cp -r hits.fasta $address/spades/bin
	mv hits.fasta $address/OUTPUT/$Out
	cd $address/spades/bin
	rm -rf assembly_*
	cd $address; echo "9" > CR.step; CFLR="N"
fi
}

G_SPADES1 ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "9" ]]; then echo "******Skipping Assembly Module with High Kmers"; else
	cd $address/spades/bin
	python spades.py -k 21,31,41,51,61,71,81,91,101,111,121 --only-assembler -s hits.fasta -o assembly_$Out
	rm -rf hits.fasta
	cd $address; echo "10" > CR.step; CFLR="N"
fi
}

G_SPADES2 ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "9" && `cat CR.step` != "10" ]]; then echo "******Skipping Assembly Mode with Lower Kmers"; else
	if [[ -d $address/spades/bin/assembly_$Out ]]
	then
		cd $address/spades/bin/assembly_$Out
	else
		cd $address/spades/bin
		rm -rf scaffolds.fasta
	fi
	if [ -s scaffolds.fasta ]
	then
		echo "SPADES ran properly with high kmers"
	else
		cd $address/OUTPUT/$Out
		echo "# Trying for $Out with lower kmers"
		rm -rf $address/spades/bin/assembly_*
		cp -r hits.fasta $address/spades/bin
		cd $address/spades/bin
		rm -rf assembly_*
		python spades.py -k 11,15,21,25,31,35,41,45,51 --only-assembler -s hits.fasta -o assembly_$Out
		rm -rf hits.fasta
	fi
	cd $address; echo "11" > CR.step; CFLR="N"
fi
}

GA ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "11" ]]; then echo "******Skipping Genome Analysis"; else
	cd $address/OUTPUT/$Out
	if [[ -s assessment.tar.gz ]]
	then
		echo "******Skipping QUAST service"
	else	
		cd $address/spades/bin/assembly_$Out
		if [ -s scaffolds.fasta ]
		then
			echo "# Analyzing draft putative genome..."
			cp -r scaffolds.fasta $address/quast
			cp -r scaffolds.fasta $address/OUTPUT/$Out
			cd $address/quast
			python metaquast.py -R $address/Reference_seqs/$Ref -o assessment scaffolds.fasta
			echo "### Compressing results..."
			tar -zcvf assessment.tar.gz assessment --remove-files
			mv assessment.tar.gz $address/OUTPUT/$Out/assessment.tar.gz
			echo "# QUAST service is finished for the file going to OUTPUT/$Out"
		else
			echo "# The proposed analysis of $Out could not continue due to problems in SPADES assembly."
		fi
	fi
	cd $address/OUTPUT/$Out
	if [ -s scaffolds.fasta ]
	then
		cd $address/OUTPUT/$Out
		echo "# Initiating ORF finding process for file going to $Out"
		rm -rf ORFs.$Out.fna
		perl $address/bb.orffinder.pl --infile=scaffolds.fasta --outfile=ORFs.$Out.fna --minlen=200 --fasta # If user wants to find orfs bigger or smaller just change parameter "minsize"
		cntg=`grep ">" scaffolds.fasta | wc -l`
		echo "$cntg" > $address/OUTPUT/$Out/cntg.nmb
		if [ -s ORFs.$Out.fna ]
		then
			ORFs=`grep ">" ORFs.$Out.fna | wc -l`
			echo "$ORFs" > $address/OUTPUT/$Out/ORFs.nmb
			gzip ORFs.$Out.fna
		else
			rm -rf ORFs.$Out.fna
			ORFs="0"
			echo "$ORFs" > $address/OUTPUT/$Out/ORFs.nmb
		fi
		echo "A total of $ORFs ORFs were found for reference $Ref, from $cntg contigs"
		gzip scaffolds.fasta
	else
		echo "# The proposed analysis of $Out could not continue due to problems in SPADES assembly."
	fi
	cd $address/spades/bin
	rm -rf assembly_$Out
	cd $address; echo "12" > CR.step; CFLR="N"
fi
}

SoftGA ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "11" ]]; then echo "******Skipping calculation of contigs and output compression"; else
	cd $address/spades/bin/assembly_$Out
	gzip scaffolds.fasta
	mv scaffolds.fasta.gz $address/OUTPUT/$Out
	cd $address/spades/bin
	rm -rf assembly_$Out
	cd $address/OUTPUT/$Out
	gzip hits.fasta
	cntg=`grep ">" scaffolds.fasta | wc -l`
	echo "$cntg" > $address/OUTPUT/$Out/cntg.nmb
	ORFs="[Warning: ORFs are not calculated in Soft version]"
	echo "$ORFs" > $address/OUTPUT/$Out/ORFs.nmb
	gzip scaffolds.fasta
	cd $address; echo "12" > CR.step; CFLR="N"
fi
}

BlastDBGen ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "8" ]]; then echo "******Skipping Preparation and generation of Blast Databases"; else
	cd $address/Buckets
	cp hits $address/Reference_seqs/$SubRef
	cd $address/Reference_seqs/$SubRef
	rm -rf list.f*; rm -rf glist; rm -rf BlastDBlist
	ls *.fa > list.fa; ls *.fasta > list.fasta; ls *.fas > list.fas; ls *.faa > list.faa; ls *.fna > list.fna; ls *.fsa > list.fsa
	cat list.f* > glist; rm -rf list.f*
	if [ -s glist ]
	then
		echo "# Recognized $SubRef files in fasta format. Making blast databases..."
		for sub in `cat glist`; do
			case $T1 in
				P|p)
					makeblastdb -in $sub -dbtype prot -out $sub.db
				;;
				N|n)
					makeblastdb -in $sub -dbtype nucl -out $sub.db
				;;
			esac
		done
	fi
	case $T1 in
		P|p)
			ls *.psq | sed 's/.psq//' | sort -k1,1 > BlastDBlist
		;;
		N|n)
			ls *.nsq | sed 's/.nsq//' | sort -k1,1 > BlastDBlist
		;;
	esac
	cd $address; echo "13" > CR.step; CFLR="N"
fi
}

Filter2 ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "13" ]]; then echo "******Skipping Second Filtering System"; else
	cd $address/Reference_seqs/$SubRef
	for sub in `cat BlastDBlist`; do
		echo "# Searching against $sub..."
		case $T1 in
			P|p)
				dsp1=`date -u "+%s"`
				blastx -db $sub -query hits -out $sub.tmp -evalue 1e-5 -strand both -max_target_seqs 1 -num_threads 4 -outfmt 6 # Parameters of reads search by blast should be specified here
				dsp2=`date -u "+%s"`
				echo $dsp2 - $dsp1 |bc -l > $sub.ft.time2
			;;
			N|n)
				dsp1=`date -u "+%s"`
				blastn -db $sub -query hits -out $sub.tmp -evalue 1e-5 -strand both -max_target_seqs 1 -num_threads 4 -outfmt 6 # Parameters of reads search by blast should be specified here
				dsp2=`date -u "+%s"`
				echo $dsp2 - $dsp1 |bc -l > $sub.ft.time2
			;;
		esac
		cat $sub.tmp | sort -k3,3 -k4,4 -n -r | awk '$3 > 90 && $4 > 25' | uniq > $sub.ft
		touch $sub.ft
		rm -rf $sub.tmp
		sed -i -e 1,1d BlastDBlist
	done
	rm -rf BlastDBlist; rm -rf *.tmp
	cd $address; echo "14" > CR.step; CFLR="N"
fi
}

SaveDBs ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "14" ]]; then echo "******Skipping the process of saving/removing databases"; else
	cd $address/Reference_seqs/$SubRef
	case $ans in
		Y|y|YES|Yes|yes)
			echo "# Databases of subreference $SubRef now saved to References_seqs folder in blastdb format. Fasta files used to make the databases have been realocated to Reference_seqs/$SubRef/Fasta_files."
			mkdir Fasta_files
			for file in `cat glist`; do
				mv $file Fasta_files
			done
			rm -rf glist
		;;
		*)
			echo "# Removing blastdb databases generated using fasta files in the subreference folder (Reference_seqs/$SubRef)"
			for file in `cat glist`; do
				case $T1 in
					P|p)
						rm -rf $file.db.p*
					;;
					N|n)
						rm -rf $file.db.n*
					;;
				esac
			done
			rm -rf glist
		;;
	esac
	cd $address; echo "15" > CR.step; CFLR="N"
fi
}

Extraction ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "15" ]]; then echo "******Skipping Extraction Module"; else
	cd $address/Reference_seqs/$SubRef
	echo "# Arranging data..."
	mv *.ft $address/OUTPUT/$Out
	mv *.time2 $address/OUTPUT/$Out
	mv hits $address/OUTPUT/$Out
	cd $address/OUTPUT/$Out
	rm -rf ext.py
	echo "#exclamationmark/usr/bin/python

import string
import sys
ListOfIds = sys.argv[1]
fastafile = sys.argv[2]

try:
    ids=open(ListOfIds, 'r')
except IOError, e:
    print \"File error: \",ListOfIds
    pass


lignes = ids.readlines()
req=[]
for ligne in lignes:
    req.append(ligne.strip())

handle = open(fastafile)

bank={}
seqIDmap={}
seq_id = handle.next()
while (seq_id[0]!=\">\"):
    seq_id = handle.next()
while True:
    try:
        seq = handle.next()
        line = handle.next()
        while (line[0]!=\">\"):
            seq = seq+line
            line = handle.next()
        bank[seq_id]=seq
        IDclean=string.split(seq_id, \" \")[0][1:].strip()
        seqIDmap[IDclean]=seq_id
        seq_id = line # for the next
    except StopIteration:
        break
# last line
bank[seq_id]=seq
seqIDmap[string.split(seq_id, \" \")[0][1:].strip()]=seq_id

handle.close()

faName=fastafile.split(\"/\")[-1]
listName=ListOfIds.split(\"/\")[-1]
subsetName=listName+\"-\"+faName
subset = open(subsetName,\"w\")
nbNF=0
for i in req:
    try:
        subset.write(seqIDmap[i].strip()+\"\substn\")
        subset.write(bank[seqIDmap[i]].strip()+\"\substn\")
    except KeyError:
        print i, \"not found in fasta\"
        nbNF+=1

subset.close()

print
print nbNF, \"IDs (listed above) from\",listName, \"have not been found in\", faName
print
print \"the Subset fasta file\", subsetName, \"is now created\"" > ext.py
	sed -i 's/exclamationmark/!/' ext.py
	sed -i 's/substn/n/' ext.py
	mkdir read_hits
	mkdir contigs
	mkdir blast_hits
	mkdir ORFs
	cd $address; echo "16" > CR.step; CFLR="N"
fi
}

PN_Prepare_SPADES ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "17" ]]; then echo "******Skipping Assembly Preparation"; else
	cd $address/OUTPUT/$Out
	rm -rf *.rev; rm -rf *.hits; rm -rf *.rev-hits
	cut -f 1 $File > $File.rev
	python ext.py $File.rev hits > $File.hits
	rm -rf $File.hits
	sq=`grep ">" $File.rev-hits | wc -l`
	echo "$sq" > $address/OUTPUT/$Out/sq.nmb
	ppm2=`expr 1000000 \* $sq / $reads`
	echo "$ppm2" > $address/OUTPUT/$Out/ppm2.nmb
	cp -r $File.rev-hits $address/spades/bin/$File.fasta
	mv $File.rev-hits read_hits/$File.fasta
	mv $File $address/OUTPUT/$Out/blast_hits
	cd $address/spades/bin
	rm -rf $address/spades/bin/assembly_*
	cd $address; echo "18" > CR.step; CFLR="N"
fi
}

PN_SPADES1 ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "18" ]]; then  echo "Skipping Assembly with High Kmers"; else
	cd $address/spades/bin
	rm -rf assembly_*
	python spades.py -k 21,31,41,51,61,71,81,91,101,111,121 --only-assembler -s $File.fasta -o assembly_$Out_$File
	cd $address; echo "19" > CR.step; CFLR="N"
fi
}

PN_SPADES2 ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "18" && `cat CR.step` != "19" ]]; then echo "******Skipping Assembly with Low Kmers"; else
	cd $address/spades/bin/assembly_$Out_$File
	if [ -s scaffolds.fasta ]
	then
		echo "# No need to try with lower kmers "
	else
		echo "# Trying for $File for $Out with lower kmers"
		cd $address/spades/bin
		sq=`grep ">" $File.fasta | wc -l`
		echo "$sq" > $address/OUTPUT/$Out/sq.nmb
		ppm2=`expr 1000000 \* $sq / $reads`
		echo "$ppm2" > $address/OUTPUT/$Out/ppm2.nmb
		cd $address/spades/bin
		rm -rf assembly_*
		python spades.py -k 9,11,13,15,17,19,21,31 --only-assembler -s $File.fasta -o assembly_$Out_$File
	fi
	cd $address; echo "20" > CR.step; CFLR="N"
fi
}

PNA ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "20" ]]; then echo "******Skipping calculation of contigs"; else
	cd $address/spades/bin/assembly_$Out_$File
	if [ -s scaffolds.fasta ]
	then
		cntg=`grep ">" scaffolds.fasta | wc -l`
		echo "$cntg" > $address/OUTPUT/$Out/cntg.nmb
		cd $address/spades/bin/assembly_$Out_$File
		echo "# SPADES worked on $File for $Out, finding $cntg contigs"
		mv scaffolds.fasta contigs.$File.fasta
		mv contigs.$File.fasta $address/OUTPUT/$Out/contigs
	else
		echo "# The proposed analysis could not continue due to problems in SPADES assembly."

		Warnings="WARNING: Did not run SPADES properly"
	fi
	cd $address/spades/bin
	rm -rf assembly_$Out_$File; rm -rf $File.fasta
	cd $address/OUTPUT/$Out 
	cd $address; echo "21" > CR.step; CFLR="N"
fi
}

PNORFs ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "21" && `cat CR.step` != "12" ]]; then echo "******Skipping ORF finding process"; else
	cd $address/OUTPUT/$Out/contigs
	if [ -s contigs.$File.fasta ]
	then
		echo "# Initiating ORF finding process"
		rm -rf ORFs.$File.fna
		perl $address/bb.orffinder.pl --infile=contigs.$File.fasta --outfile=ORFs.$File.fna --minlen=300 --fasta # If user wants to find orfs bigger or smaller just change parameter "minsize"
		if [ -s ORFs.$File.fna ]
		then
			ORFs=`grep ">" ORFs.$File.fna | wc -l`
			echo "$ORFs" > $address/OUTPUT/$Out/ORFs.nmb
			mv ORFs.$File.fna $address/OUTPUT/$Out/ORFs
		else
			ORFs="0"
			echo "$ORFs" > $address/OUTPUT/$Out/ORFs.nmb
			rm -rf ORFs.$File.fna
		fi
		if [[ -s cntg.nmb ]]
		then 
			cntg=`cat cntg.nmb`
		fi
		echo "A total of $ORFs ORFs were found for file $File, from $cntg contigs"
	fi
	echo "$ORFs" > $address/OUTPUT/$Out/ORFs.nmb
	cd $address; echo "22" > CR.step; CFLR="N"
fi
}

CleaningTheMess ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "21" && `cat CR.step` != "22"  && `cat CR.step` != "12" ]]; then echo "******Skipping Self Organizing Module"; else
	cd $address/Buckets
	rm -rf hits
	cd $address
	case $T1 in
		P|p|N|n)
			cd $address/OUTPUT/$Out
			if [[ -s log.tmps ]]
			then
				touch log.tmps
			else
				cat tmp1.* > log.tmps
				rm -rf tmp1.*
			fi
			cd $address/OUTPUT/$Out
			mv *.ft blast_hits
			cd $address/OUTPUT/$Out/blast_hits
			ls *.ft > list
			for file in `cat list`; do
				mv $file ${file/.f*/.blast.tsv}
			done
			rm -rf list
		
			cd $address/OUTPUT/$Out/contigs	
			ls *.fasta > contigs_list
			for file in `cat contigs_list`; do
				mv $file ${file/.f*/.fasta}
			done
			rm -rf contigs_list

			cd $address/OUTPUT/$Out/ORFs
			ls *.fna > orfs_list
			for file in `cat orfs_list`; do
				mv $file ${file/.f*/.fna}
			done
			rm -rf orfs_list

			cd $address/OUTPUT/$Out/read_hits
			ls *.fasta > hits_list
			for file in `cat hits_list`; do
				mv $file ${file/.f*/.hits.fasta}
			done
			rm -rf hits_list

			cd $address/OUTPUT/$Out
			rm -rf log.header log.tmp1; rm -rf cont_log*; rm -rf ORF_log*
			echo "-----------------------------------------------------------------------------
Subref_database|Hits_seq.|Ppm|Contigs|ORFs|Blast_Time|SPADES_Time|Total_SubRef_Time|Status" > log.header
			cat log.header log.tmps > log.tmp1
			sed -i 's/|/\t/g' log.tmp1
			cut -f 4 log.tmp1 > cont_log$(( n %= 100001))
			cut -f 5 log.tmp1 > ORF_log$(( n %= 100001))
			mv cont_log* $address; mv ORF_log* $address
			rm -rf log.header
		;;
	esac
	rm -rf $address/spades/bin/assembly_*
	cd $address/OUTPUT/$Out
	d2=`date -u "+%s"`
	if [[ $TimeLoss == "Y" ]]
	then
		d3="at least $(echo "$d2 - $d1" |bc -l)"
	else
		d3=$(echo "$d2 - $d1" |bc -l)
	fi
	echo "Total time for Reference $Ref with Subreference $SubRef: $d3 seconds" > fulltime.tmp
	cd $address/OUTPUT/$Out
	cat log.tmp1 fulltime.tmp > log.tmp
	cat Log.txt log.tmp > FL
	rm -rf Log.txt
	mv FL Log.txt
	mv hits hits.fasta
	gzip hits.fasta
	rm -rf hits hits.fasta
	if [[ -s scaffolds.fasta ]]
	then 
		gzip scaffolds.fasta
	fi
	rm -rf scaffolds.fasta log.tmp log.tmps log.tmp1 fulltime.tmp ext.py list; rm -rf *.time2; rm -rf *.rev; rm -rf *.hits
	cd $address; echo "23" > CR.step; CFLR="N"
fi
}

PrintResults ()
{
cd $address; if [[ "$CFLR" == "Y" && `cat CR.step` != "23" ]]; then echo "******Skipping Print Results"; else
	cd $address/OUTPUT/$Out
	rm -rvf tmp1.*
	cd $address
	d2=`date -u "+%s"`
	if [[ $TimeLoss == "Y" ]]
	then
		d3="at least $(echo "$d2 - $d1" |bc -l)"
	else
		d3=$(echo "$d2 - $d1" |bc -l)
	fi
	if [ -s cont_log* ]
	then
		cat cont_log* > c_int
		sed -i '/[a-z]/d' c_int
		cntg=`awk '{s+=$1} END {print s}' c_int`
		echo "$cntg" > $address/OUTPUT/$Out/cntg.nmb
	else
		cd $address/OUTPUT/$Out
		if [[ -s $address/OUTPUT/$Out/cntg.nmb ]]
		then 
			cntg=`cat $address/OUTPUT/$Out/cntg.nmb`
		fi
		cd $address
	fi
	if [[ `cat CR.mode` == "Soft" ]]; 
	then
		ORFs="[Warning: ORFs are not calculated in Soft version]"
	else
		if [ -s ORF_log* ]
		then
			cat ORF_log* > o_int
			sed -i '/[a-z]/d' o_int
			ORFs=`awk '{s+=$1} END {print s}' o_int`
			echo "$ORFs" > $address/OUTPUT/$Out/ORFs.nmb
		else
			cd $address/OUTPUT/$Out
			if [[ -s $address/OUTPUT/$Out/ORFs.nmb ]]
			then 
				ORFs=`cat $address/OUTPUT/$Out/ORFs.nmb`
			fi
			cd $address
		fi
	fi
	cd $address/OUTPUT/$Out
	if [[ -s $address/OUTPUT/$Out/ppm1.nmb ]]
	then
		ppm1=`cat $address/OUTPUT/$Out/ppm1.nmb`
	fi
	if [[ -s $address/Buckets/reads.nmb ]]
	then
		reads=`cat $address/Buckets/reads.nmb`
	else
		if [[ -s $address/OUTPUT/$Out/reads.nmb ]]
		then
			reads=`cat $address/OUTPUT/$Out/reads.nmb`
		fi
	fi
	cd $address
	echo "$Out|$R1|$T1|$T2|$Ref|$SubRef|$d3|$reads|$buckets|$ppm1|$cntg|$ORFs" > par.time
	cat time.log par.time > time.log2; rm -rf time.log par.time; mv time.log2 time.log
	sed -i 's/|/\t/g' time.log
	rm -rf c_int o_int; rm -rf ORF_log*; rm -rf cont_log*
	echo -e "\n [BEAF12.01.17 worked in $R1 with reference as $Ref (output as $Out) for $d3 seconds] \n"
	sed -i -e 1,1d config.kp
	case $Keep in
		Y|y)
			touch Buckets
		;;
		*)
			rm -rf Buckets
		;;
	esac
	cd $address; echo "24" > CR.step; CFLR="N"
fi
}

ErrorRevision ()
{
	cd $address; echo "25" > CR.step; CFLR="N"
	echo "Starting Error Revision"
	while read T1 T2 R1 R2 Ref SubRef Out; do
		if [[ -d $address/OUTPUT ]]
		then
			touch OUTPUT
		else
			mkdir $address/OUTPUT
		fi
		cd $address/OUTPUT
		if [[ -d $Out ]]
		then
			touch $Out
		else
			mkdir $Out
			cd $Out
			echo "This line was skipped in the config.file" > Error.msg
		fi
	done < config.file
	
	cd $address/OUTPUT
	rm -rf list
	ls > list
	sed -i '/list/d' list
	sed -i '/Errors/d' list
	if [ -d Errors ]
	then
		mv Errors Errors.old
		mkdir Errors
		mv Errors.old Errors
	else
		mkdir Errors
	fi
	for folder in `cat list`; do        
	        cd $address/OUTPUT/$folder
	        if [ -s scaffolds.fasta.gz ]
		then
			touch scaffolds.fasta.gz
	        else
	                if [ -d contigs ]
	                then
	                        find "contigs" -type f -exec echo Found file {} \; > test.txt
	                        if [ -s test.txt ]
				then
					rm -rf test.txt
	                        else
					cd $address/OUTPUT
					mv $folder Errors
					rm -rf test.txt
	                        fi
	                else
				cd $address/OUTPUT
				mv $folder Errors
	                fi
	        fi
	        cd $address/OUTPUT
	done
	rm -rf list
	cd $address/OUTPUT/Errors
	ls > redolist
	sed -i '/redolist/d' redolist
	sed -i '/Errors.old/d' redolist
	mv redolist $address
	cd $address
	for out in `cat redolist`; do
		grep "$out" config.file > $out.tmpredo
	done
	cat *.tmpredo > config.redo
	rm -rf *.tmpredo
	mv config.redo $address/OUTPUT/Errors
	rm -rf redolist
}

	# ======================================================================================================================================================================================== #
	# =====================================================================================PIPELINE BEAF====================================================================================== #
	# ======================================================================================================================================================================================== #


BEAF ()
{
	echo "##### Running BEAF, full version #####"
	d0=`date -u "+%s"`
	address=$(dirname $(readlink -f $0))
	
	make_kp
	Check
	TimeHeader
	while read T1 T2 R1 R2 Ref SubRef Out Keep; do
		d1=`date -u "+%s"`
		cd $address
		echo -e "\n# Starting work in file $R1 with $Ref as reference, going to $Out\n"
		if [ -d Buckets ]
		then
			cd $address/Buckets
			ls *.bk > buckets_list.txt
			if [ -s buckets_list.txt ];
			then
				echo "# Using previous buckets"
			else
				cd $address
				Trim
				QAnConversion
			fi
		else
			mkdir $address/Buckets
			Trim
			QAnConversion
		fi
		if [[ -d $address/OUTPUT/$Out ]]
		then
			echo "Folder $Out already exists, continuing work..."
		else
			if [[ -d $address/OUTPUT ]]
			then
				touch OUTPUT
			else
				mkdir $address/OUTPUT
			fi
			mkdir $address/OUTPUT/$Out
			echo "Folder $Out created in OUTPUT"
		fi
		BucketEngine
		cd $address/Buckets
		buckets=`ls *.bk | wc -l`
		Filter1
		PreLogGen
		cd $address
		cd $address/Buckets
		case $T1 in
			G|g)
				if [ -s hits ]
				then
					G_Prepare_SPADES
					G_SPADES1
					G_SPADES2
					GA
					cd $address/spades/bin
					rm -rf assembly_$Out
					cd $address/OUTPUT/$Out
					gzip hits.fasta
				else
					echo "# The proposed analysis could not continue due to its lacking of homology between provided sequences and reference genome."
					rm -rf hits
					echo "12" > CR.step
				fi
			;;
			P|p|N|n)
				if [ -s hits ]
				then
					echo "# Submitting to Blast per subreference family..."
					BlastDBGen
					Filter2
					SaveDBs
					Extraction
					if [[ -s $address/Buckets/reads.nmb ]]
					then
						reads=`cat $address/Buckets/reads.nmb`
					else
						if [[ -s $address/OUTPUT/$Out/reads.nmb ]]
						then
							reads=`cat $address/OUTPUT/$Out/reads.nmb`
						fi
					fi	
					cd $address/OUTPUT/$Out
					if [[ "$CFLR" == "Y" ]]
					then
						if [[ -s list ]]
						then
							touch list
						else
							ls *.ft > list
						fi
					else
						ls *.ft > list
					fi
					for File in `cat list`; do
						echo "# Working on file $File for $Out..."
						cd $address/OUTPUT/$Out
						BTime="0"
						STime="0"
						TTime="0"
						d4=`date -u "+%s"`
						ppm2="0"
						cntg="0"
						sq="0"
						Warnings="[OK]"
						ORFs="0"
						if [ -s $File ]
						then
							PN_Prepare_SPADES
							PN_SPADES1
							PN_SPADES2
							PNA
						else
							echo "# $File in $Out did not reach the minimum criteria to be considered homologus"
							cd $address/OUTPUT/$Out
							Warnings="WARNING: Did not reach minimum criteria to be considered homologus"
							rm -rf $File
						fi
						PNORFs
						d5=`date -u "+%s"`
						cd $address/OUTPUT/$Out
						if [ -s $File.time2 ]
						then
							BTime=`cat $File.time2`
						fi
						if [[ $TimeLoss == "Y" && $CFLR == "Y" ]]
						then
							STime=$(echo $d5 - $d4 |bc -l)
							TTime="at least $(echo $BTime + $STime |bc -l)"
							STime="at least $(echo $d5 - $d4 |bc -l)"
						else
							STime=$(echo $d5 - $d4 |bc -l)
							TTime=$(echo $BTime + $STime |bc -l)
						fi
						rm -rf *.time2
						cd $address/OUTPUT/$Out
						if [[ -s $address/OUTPUT/$Out/sq.nmb ]]
						then
							sq=`cat $address/OUTPUT/$Out/sq.nmb`
						fi
						if [[ -s $address/OUTPUT/$Out/ppm2.nmb ]]
						then
							ppm2=`cat $address/OUTPUT/$Out/ppm2.nmb`
						fi
						if [[ -s $address/OUTPUT/$Out/cntg.nmb ]]
						then
							cntg=`cat $address/OUTPUT/$Out/cntg.nmb`
						fi
						if [[ -s $address/OUTPUT/$Out/ORFs.nmb ]]
						then
							ORFs=`cat $address/OUTPUT/$Out/ORFs.nmb`
						fi
						cd $address/OUTPUT/$Out
						echo "${File/.f*/}|$sq|$ppm2|$cntg|$ORFs|$BTime|$STime|$TTime|$Warnings" > tmp1.$File
						sed -i -e 1,1d list
					done
					cd $address/OUTPUT/$Out; rm -rf list
				else
					echo "# The proposed analysis could not continue due to its lacking of homology between provided sequences and reference sequence."
					ppm2="0"
					echo "$ppm2" > $address/OUTPUT/$Out/ppm2.nmb
					cntg="0"
					echo "$cntg" > $address/OUTPUT/$Out/cntg.nmb
					ORFs="0"
					echo "$ORFs" > $address/OUTPUT/$Out/ORFs.nmb
					echo "22" > CR.step
				fi
				cd $address
			;;
		esac
		CleaningTheMess
		PrintResults
		cd $address/OUTPUT/$Out; rm -rf list; rm -rf *.nmb
		CFLR="N"; TimeLoss="N"
	done < config.kp
	rm -rf *.kp *.nmb
	rm -rf par.time
	
	ErrorRevision
	
	cd $address; 
	CFLR="N"
	rm -rf CR.step CR.mode

	d99=`date -u "+%s"`
	dtotal=$(echo "$d99 - $d0" |bc -l)
	echo "BEAF1011.65 worked for $dtotal seconds, ending at $(date "+%X")."
}

SoftBEAF ()
{
	echo "##### Running Soft BEAF, a faster but simpler version of the BEAF software #####"
	d0=`date -u "+%s"`
	address=$(dirname $(readlink -f $0))
	
	make_kp
	Check
	SoftTimeHeader
	while read T1 T2 R1 R2 Ref SubRef Out Keep; do
		d1=`date -u "+%s"`
		cd $address
		echo -e "\n# Starting work in file $R1 with $Ref as reference, going to $Out\n"
		if [ -d Buckets ]
		then
			cd $address/Buckets
			ls *.bk > buckets_list.txt
			if [ -s buckets_list.txt ];
			then
				echo "# Using previous buckets"
			else
				cd $address
				CopyFile
				SoftMergeRename
			fi
		else
			mkdir $address/Buckets
			CopyFile
			SoftMergeRename
		fi
		if [[ -d $address/OUTPUT/$Out ]]
		then
			echo "Folder $Out already exists, continuing work..."
		else
			if [[ -d $address/OUTPUT ]]
			then
				touch OUTPUT
			else
				mkdir $address/OUTPUT
			fi
			mkdir $address/OUTPUT/$Out
			echo "Folder $Out created in OUTPUT"
		fi
		BucketEngine
		cd $address/Buckets
		buckets=`ls *.bk | wc -l`
		Filter1
		PreLogGen
		cd $address/Buckets
		case $T1 in
			G|g)
				if [ -s hits ]
				then
					G_Prepare_SPADES
					G_SPADES2
					SoftGA
				else
					echo "# The proposed analysis could not continue due to its lacking of homology between provided sequences and reference genome."
					rm -rf hits
					echo "12" > CR.step
				fi
			;;
			P|p|N|n)
				if [ -s hits ]
				then
					echo "# Submitting to Blast per subreference family..."
					BlastDBGen
					Filter2
					SaveDBs
					Extraction
					if [[ -s $address/Buckets/reads.nmb ]]
					then
						reads=`cat $address/Buckets/reads.nmb`
					else
						if [[ -s $address/OUTPUT/$Out/reads.nmb ]]
						then
							reads=`cat $address/OUTPUT/$Out/reads.nmb`
						fi
					fi					
					cd $address/OUTPUT/$Out
					if [[ -s list ]]
					then
						touch list
					else
						ls *.ft > list
					fi
					for File in `cat list`; do
						cd $address/OUTPUT/$Out
						echo "# Working on file $File for $Out..."
						BTime="0"
						STime="0"
						TTime="0"
						d4=`date -u "+%s"`
						ppm2="0"
						cntg="0"
						sq="0"
						Warnings="[OK]"
						if [ -s $File ]
						then
							PN_Prepare_SPADES
							PN_SPADES2
							PNA
						else
							echo "# $File in $Out did not reach the minimum criteria to be considered homologus"
							cd $address/OUTPUT/$Out
							Warnings="WARNING: Did not reach minimum criteria to be considered homologus"
							rm -rf $File
						fi
						d5=`date -u "+%s"`
						cd $address/OUTPUT/$Out
						if [ -s $File.time2 ]
						then
							BTime=`cat $File.time2`
						fi
						if [[ $TimeLoss == "Y" && $CFLR == "Y" ]]
						then
							STime=$(echo $d5 - $d4 |bc -l)
							TTime="at least $(echo $BTime + $STime |bc -l)"
							STime="at least $(echo $d5 - $d4 |bc -l)"
						else
							STime=$(echo $d5 - $d4 |bc -l)
							TTime=$(echo $BTime + $STime |bc -l)
						fi
						rm -rf *.time2
						cd $address/OUTPUT/$Out
						if [[ -s $address/OUTPUT/$Out/sq.nmb ]]
						then
							sq=`cat $address/OUTPUT/$Out/sq.nmb`
						fi
						if [[ -s $address/OUTPUT/$Out/ppm2.nmb ]]
						then
							ppm2=`cat $address/OUTPUT/$Out/ppm2.nmb`
						fi
						if [[ -s $address/OUTPUT/$Out/cntg.nmb ]]
						then
							cntg=`cat $address/OUTPUT/$Out/cntg.nmb`
						fi
						echo "${File/.f*/}|$sq|$ppm2|$cntg|NotTested|$BTime|$STime|$TTime|$Warnings" > tmp1.$File
						sed -i -e 1,1d list
					done
					cd $address/OUTPUT/$Out; rm -rf ORFs; rm -rf list; rm -rf *.nmb
				else
					echo "# The proposed analysis could not continue due to its lacking of homology between provided sequences and reference sequence."
					ppm2="0"
					echo "$ppm2" > $address/OUTPUT/$Out/ppm2.nmb
					cntg="0"
					echo "$cntg" > $address/OUTPUT/$Out/cntg.nmb
					echo "21" > CR.step
				fi
				cd $address
			;;
		esac
		CleaningTheMess
		PrintResults
		CFLR="N"; TimeLoss="N"
		cd $address/OUTPUT/$Out; rm -rf list; rm -rf *.nmb
	done < config.kp
	rm -rf *.kp *.nmb 
	rm -rf par.time
	
	ErrorRevision

	cd $address; 
	CFLR="N"
	rm -rf CR.step CR.mode

	d99=`date -u "+%s"`
	dtotal=$(echo "$d99 - $d0" |bc -l)
	echo "BEAF1011.65 worked for $dtotal seconds, ending at $(date "+%X")."
}

	# ======================================================================================================================================================================================== #
	# ========================================================================================MAIN============================================================================================ #
	# ======================================================================================================================================================================================== #

Main ()
{
	case $1 in
		S|s|-s|-S|Soft|soft|SOFT|-Soft|-soft)
			echo "Soft" > CR.mode
			SoftBEAF
		;;
		*)
			echo "FullVersion" > CR.mode
			BEAF
		;;
	esac
}

	# ======================================================================================================================================================================================== #
	# =====================================================================================PROGRAM START====================================================================================== #
	# ======================================================================================================================================================================================== #


address=$(dirname $(readlink -f $0))
CFLR="N"
if [ -s CR.step ]
then
	echo "A file from your previous run was detected, indicating that your last run of BEAF was interrupted before the program finished. BEAF is capable of detecting at which step the program was interrupted and continuing interrupted runs, but note that if you removed files from the BEAF folder after the program was abruptly interrupted, the new run can be disrupted. Do you want to continue from where the program stopped? (Y/N)"
	read CFLR # continue from last run
	case $CFLR in
		Y|y|Yes|yes|continue|Continue)
			CFLR="Y"
			TimeLoss="Y"
			echo "BEAF will continue from last run."
			ver=`cat CR.mode`
			Main $ver
		;;
		*)
			CFLR="N"
			rm -rf CR.step CR.mode; rm -rf *.kp; rm -rf Buckets
			echo "BEAF will start a new run. 

How do you want to run BEAF?
-s: run Soft BEAF, a faster but simpler version of the BEAF software
-b: run full version, with all utilities"
			read ver
			Main $ver
		;;
	esac
else
	if [ -s config.kp ]
	then
		echo "A config.kp file from your previous run was detected, indicating that your last run of BEAF was interrupted before the program finished. Although BEAF couldn't detect at which point the program was interrupted (because the file 'CR.step' was removed), it was able to determine which OUTPUTs were already generated. BEAF is capable of continuing interrupted runs, but note that if you removed files from the BEAF folder after the program was abruptly interrupted, the new run can be disrupted. Do you want to continue from where the program stopped? (Y/N)"
		read SICLR # Should I Continue Last Run?
		case $SICLR in
			Y|y|Yes|yes|continue|Continue)
				CFLR="Y" # Continue From Last Run
				TimeLoss="Y"
				echo "1" > CR.step
				echo "BEAF will continue from last run."
				ver=`cat CR.mode`
				Main $ver
			;;
			*)
				CFLR="N" # Continue From Last Run
				rm -rf CR.step CR.mode; rm -rf *.kp; rm -rf Buckets
				echo "BEAF will start a new run. 
	
	How do you want to run BEAF?
	-s: run Soft BEAF, a faster but simpler version of the BEAF software
	-b: run full version, with all utilities"
				read ver
				Main $ver
			;;
		esac
	else
		echo "How do you want to run BEAF?
-s: run Soft BEAF, a faster but simpler version of the BEAF software
-b: run full version, with all utilities"
		read ver
		Main $ver
	fi
fi
