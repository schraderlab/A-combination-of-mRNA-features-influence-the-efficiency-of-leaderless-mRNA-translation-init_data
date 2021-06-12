# A-combination-of-mRNA-features-influence-the-efficiency-of-leaderless-mRNA-translation-init_data
******************************************************************************************************
*																									 *
*										Caulobacter crescentus										 *
*																									 *
******************************************************************************************************
1. The CDS were mapped using trypsin digestion and ribosome profiling (The Coding and Noncoding 
	Architecture of the Caulobacter crescentus Genome) available in the latest C. crescentus Genbank 
	entry (NC_011916.1). 


2. The operon annotation was downloaded from (The global regulatory architecture of transcription 
	during the Caulobacter cell cycle), and operons that were within 50 nts were combined into a single operon. 
	
3. Transcription start sites (TSSs) were mapped using 5’RACE and RNA-seq (The Coding and Noncoding 
	Architecture of the Caulobacter crescentus Genome).
	the following are the TSS files used:
	a. _6_TSS_plus
	b. _6_TSS_minus
	
4. The mapped TSSs were assigned to the operons if they were present within 300 nucleotides upstream from the 
	start codon of the first gene in the operon using built-in python script "retrieving_TSS_sites.py" using following input files:
	a. operon_list_plus.txt
	b. operon_list_minus.txt
	the output files were:
	a. plus_all_tss.txt
	b. minus_all_tss.txt
	
5. For all operons for which no mapped TSS was found within 300 nts upstream, they were copied into new files:
	a. operon_list_plus_noTSS.txt
	b. operon_list_minus_noTSS.txt
	The following RNA processed ends data was then used to assign potential TSSs to these operons using script 
	"retrieving_processed_TSS_sites.py":
	a. _rna-seq_5ends_f_nonzero
	b. _rna-seq_5ends_r_nonzero
	the output files were:
	a. plus_noTSS_all_potential_tss.txt
	b. minus_noTSS_all_potential_tss.txt
	
6. the outputs from previous two steps were further combined and processed in excel to create following two files:
	a. genes_plus.txt
	b. genes_minus.txt
	These files were used to retrieve gene sequences based on TSS location using built-in python script 
	"retrieving_entire_sequence.py".

7.	Further processing was done in excel for various analysis. Translation initiation regions (TIRs) for all genes 
	were retrieved using excel functions.

8.	dGunfold calculations were then for TIRS using the program found on github
	(https://github.com/schraderlab/dGunfold_program/tree/main/sequences_dGunfold) with following files as input:
	a. caulobacter_sequences.txt
	b. caulobacter_start_positions.txt
	
9. The final data is available as an excel file "C_crescentus.xlsx"







******************************************************************************************************
*																									 *
*										Mycobacterium species									     *
*																									 *
******************************************************************************************************
1. The following genes, TSS and TE data files were provided by Joe Wade lab: 
	a. _6_TSS_minus_H37Rv 
	b. _6_TSS_minus_smeg 
	c. _6_TSS_plus_H37Rv
	d. _6_TSS_plus_smeg
	e. M_tuberculosis_H37Rv_genes.txt
	f. M_smegmatis_genes.txt
	g. TE_Msmeg.txt
	
2. The predicted operon data for H37Rv and smegmatis were copied from OperonDB(http://operondb.ccb.jhu.edu/cgi-bin/taxon_list.cgi)
	into following text files:
	a. operon_data_H37Rv.txt
	b. operon_data_smegmatis.txt

3. The operon data was processed and prepared for further analysis using built-in python script 
	(operon_annotation_from_text_file_directly.py) into following files:
	a. operon_data_H37Rv_output.txt
	b. operon_data_smegmatis_output.txt
	
4.  The genes file (provide by Joe Wade lab) and operon data output files were combined in Excel to generate following files
	for further analysis:
	a. operon_list_H37Rv.txt
	b. operon_list_smeg.txt

5. The mapped TSSs were then assigned to the operons and sequences were retrieved for dGunfold analysis 
	using the following built-in python scripts: 
	a. retrieving_all_based_on_operons_for_H37Rv.py
	b. retrieving_all_based_on_operons_for_smeg.py
	the output files were:
	a. H37Rv_all.txt
	b. smegmatis_all.txt
	
6. The output files contain translation initiation region (TIR) for each sequence which is then used for 
	dGunfold calculations using the program found on github
	(https://github.com/schraderlab/dGunfold_program/tree/main/sequences_dGunfold) with following files as input:
	for tuberculosis/H37Rv:
	a. H37Rv_sequences.txt
	b. H37Rv_start_positions.txt
	for smegmatis:
	a. smegmatis_sequences.txt
	b. smegmatis_start_positions.txt
	
7. All the data is combined into an excel file for further analysis. the excel files are:
	a. M_tuberculosis.xlsx
	b. M_smegmatis.xlsx










******************************************************************************************************
*																									 *
*										Haloferax vocanii											 *
*																									 *
******************************************************************************************************
1. Following CDS/5'UTR distances/other RNAs mapped data were available:
	a. features_list.txt
	b. leaderless_supp.txt from the paper Genome-wide identification of transcriptional start sites in 
		the haloarchaeon Haloferax volcanii based on differential RNA-Seq (dRNA-Seq)

2. The following annotated ribosome profiling data file was available from DiRuggiero lab
	a. ribo.txt

3. The following operon data file was available from the paper Unbiased Map of Transcription Termination Sites in 
	Haloferax volcanii reveals manifold termination patterns
	a. operon.txt

4. The file "features_list.txt" was processed to remove duplicate entries using built-in python script "features.py".
	The output file was "locus_tag.txt"

5. All the above files were used (except "features_list.txt") to generate comprehensive and annotated map of H. volcanii
	and also retrieve TIR regions for dGunfold calculations using built-in python script "features_2.py".
	The output file is "features_2_output.txt"

6. dGunfold calculations were done using using the program found on github
	(https://github.com/schraderlab/dGunfold_program/tree/main/sequences_dGunfold) with following files as input:
	a. haloferax_sequences.txt
	b. haloferax_start_positions.txt

7. An excel file of final version is available for further analysis:
	a. haloferax_final.xlsx









******************************************************************************************************
*																									 *
*										Mouse mitochondria											 *
*																									 *
******************************************************************************************************
1. The "features.txt" file having mouse mitochondrial genome data was downloaded from Gene Expression Omnibus
	Sample GSM3764756

2. Since TSS data was not available, TIRs in-case of leadered mRNAs consisted of 50 nts such that 25 was upstream 
	and 25 downstream and for leaderless it comprised of upto 50 nts from the start codons.
	The sequences were retrieved using built-in python script "retrieving_all.py".
	
3. The TIRs (50 nts) sequences retrieved were used to computer dGunfold using program
	(https://github.com/schraderlab/dGunfold_program/tree/main/sequences_dGunfold) with following files as input:
	a. mouse_mito_sequences.txt
	b. mouse_mito_start_positions.txt

4. The final data used is provided in the excel file "mouse_final.xlsx" for further analysis.




	
******************************************************************************************************
*																									 *
*								C. crescentus reporters												 *
*																									 *
******************************************************************************************************
1. dGunfold of C. crescentus reporters was computed using program 
	(https://github.com/schraderlab/dGunfold_program/tree/main/sequences_dGunfold) with following files as input:
	a. repoter_sequences.txt
	b. repoter_start_positions.txt
	The annotated data can be found in Table S1 of the main paper. The order of the sequences is unchanged.



******************************************************************************************************
*																									 *
*										notes 												 		 *
*																									 *
******************************************************************************************************
1. all files for above mentioned sections can be found in their respective folders

2. The TIR sequences were extracted from all ORFS using 50 nts (25 nts upstream of start codon 
	and 25 nts downstream from start codon). 
	If the 5’ upstream untranslated region (UTR) was less than 25 nts, then 50 nts from transcription start 
	site was used for all TIR calculations. 
	For all the genes known to be present in the operons internally, 50 nts (25 nts upstream of start codon 
	and 25 nts downstream from start codon) were retrieved irrespective of TSS assigned or not.
	
	
3. In case of operons to which multiple TSS were assigned, if the TSS sites were very close to each other such that 
	the first gene in the operon was present at a distance of less than 25 nts from each of those TSS, multiple forms 
	were assigned to that gene. The farthest TSS was referred as form 1 followed by second farthest as form 2 and so on. 
	Sequence for each form was retrieved in such cases for ΔGunfold calculations.

4. For mRNAs labelled as leaderless but were also present in an operon internally, those were considered as both – in operon 
	and leaderless (processed) forms, and respective sequences were retrieved for ΔGunfold calculations.

5. For mouse mitochondria; since mapped TSS data was not available, 50 nts from start codon site for all the CDS labelled 
	as leadereless were retrieved and for all those labelled as leadered - 50 nts (25 nts upstream of start codon 
	and 25 nts downstream from start codon) were retrieved 





******************************************************************************************************
*																									 *
*										end 												 		 *
*																									 *
******************************************************************************************************