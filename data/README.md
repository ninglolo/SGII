# SGII
## Data sources 
+ Download lncRNA interaction data with other molecules from the NPInter database.

Download address: http://bigdata.ibp.ac.cn/npinter4/download/file/lncRNA_interaction.txt.gz

+ Download the sequences of lncRNAs from the NONCODE database.

Download address: http://www.noncode.org/datadownload/ncrna_NONCODE[v5.0].fasta.zip

+ Download protein and protein interaction data from the BIOGRID database.

Download address: https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.2.192/BIOGRID-ORGANISM-4.2.192.tab3.zip


## Files
+ eng.csv: Minimum free energy of secondary structure of sequences.
+ ncName_ncID_transID.csv: Name, ID, and transcript ID of lncRNAs.
+ transcripts_seq.fasta: Transcriptional sequence of lncRNAs.
+ LPI.csv: LncRNA-protein interaction (LPI) data.
+ PPI.csv: Protein-protein interaction (PPI) data.
+ LPPI.csv: LncRNA-protein-protein interaction (LPPI) data, composed of LPI.csv and PPI.csv. In order to prevent the protein and lncRNA from having the same name, we marked all the proteins with '<####>' as their name suffix when constructing LPPI.csv
