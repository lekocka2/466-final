# final project
# find cpg islands (proj 2)

# imports
import pandas as pd
# import csv


class genome():
    def __init__(self, file_FASTA, file_GFF):
        self.fastaFile = file_FASTA
        self.gffFile = file_GFF
        
        # from FASTA
        self.defLine = ""
        self.senseSeq = ""
        self.antisenseSeq = ""
        self.islandsDF = pd.DataFrame()
        
        # from GFF3
        self.annotDict = {}
        self.annotDF = pd.DataFrame()
        
        self.windowSize = 1000
        
        self.summaryDF = pd.DataFrame()
        self.overlapsDF = pd.DataFrame()
        
        self.isl_gene = 0
        self.isl_promoter = 0


    def readFASTA(self): # reads in the file
        file_handler = open(self.fastaFile)
        seq_file = file_handler.read()
        defLine = seq_file.split('\n')[0][0:] # split by new line, take only first line
        seq_name = defLine[0:defLine.index(' ')] # get the seq name at the beg of that line
        desc = defLine[len(seq_name):] # get the desc after the seq name on that first line
        sequence = ''.join(seq_file.split('\n')[1:]) # join all the other lines into one string
        self.defLine = defLine
        self.senseSeq = sequence
        
        
    def getAntisenseSeq(self): # gets the reverse complement of the sequence
        # find the compliments to the nucleotides in sense strand
        nuc_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'} # dictionary for getting compliment
        antisenseSeq = "" # initialize empty string
        for nuc in self.senseSeq:
            antisenseSeq += nuc_dict[nuc]
        
        self.antisenseSeq = antisenseSeq # save new sequence
        
        # TEST
        # print(self.senseSeq[int(len(self.senseSeq)/2)      :   int(len(self.senseSeq)     / 2)  +10])
        # print(antisenseSeq[int(len(self.antisenseSeq)/2)  :   int(len(self.antisenseSeq) / 2)  +10])


    def readGFF3(self):
        GFF_dict = {}
        with open(self.gffFile) as file_handler:
            
            for j in range(8):
                next(file_handler) # this skips the first 8 lines
            
            for line in file_handler:
                
                if line[0] == "#":
                    pass
                    
                else: 
                    l = line.split('\t') # split the line by tabs, returns list of elements in the line
                    subunit_desc = l[8]  # get the last item in the line
                    subunitID = subunit_desc.split(";")[0] # gets only the ID from the last item in the line
                    
                    if subunitID[0:7] == "ID=gene":
                        # subunitID = subunitID.split(":")[1]
                        subunit_dict = {}
                        
                        # create a subunit dict for the gene info
                        subunit_dict[subunitID] = [l[2], l[3], l[4], l[6], l[8]]
                        
                        # add to the unit dict
                        GFF_dict[subunitID] = [subunit_dict] # create list of values with first dict so far
                        
                    # all the subunits should always come after their respective gene, so this should work
                    if subunitID[0:6] == "Parent" or subunitID[0:4] == "ID=t":
                        subunit_dict = {}
                        # create a subunit dict for the subunit info
                        subunit_dict[subunitID] = [l[2], l[3], l[4], l[6], l[8]]
                        
                        lastKeyAdded = list(GFF_dict)[-1] # this gets the last key that was just added to the unit dict
                        GFF_dict[lastKeyAdded].append(subunit_dict) # append the new subunit_dict to the list of values
        
        self.annotDict = GFF_dict
                        
####### TESTING readGFF3 ############################################
        # print(GFF_dict["ID=gene:ENSG00000284986"])

# [{'ID=gene:ENSG00000284986': 
#   ['pseudogene', '52452', '73374', '-', 'ID=gene:ENSG00000284986;biotype=transcribed_unprocessed_pseudogene;description=phosphoglucomutase 5 (PGM5) pseudogene;gene_id=ENSG00000284986;logic_name=havana_homo_sapiens;version=1\n']}, 
# {'ID=transcript:ENST00000646046': 
#   ['pseudogenic_transcript', '52452', '73374', '-', 'ID=transcript:ENST00000646046;Parent=gene:ENSG00000284986;biotype=transcribed_unprocessed_pseudogene;tag=basic;transcript_id=ENST00000646046;version=1\n']}, 
# {'Parent=transcript:ENST00000646046': 
#   ['exon', '52452', '52614', '-', 'Parent=transcript:ENST00000646046;Name=ENSE00003824949;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00003824949;rank=2;version=1\n']}, 
# {'Parent=transcript:ENST00000646046': 
#   ['exon', '73144', '73374', '-', 'Parent=transcript:ENST00000646046;Name=ENSE00003820512;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00003820512;rank=1;version=1\n']}] 

## THIS WORKS!! YAY
#####################################################################


    def createAnnotDF(self):
        df =  pd.DataFrame((
            [subkey, key] + value
            for key, records in self.annotDict.items()
            for record in records
            for subkey, value in record.items()),
            columns=['subunit_ID', 'gene_ID', 'biotype', 'start_index', 'end_index', 'strand', 'desc'])
    
        df['gene_ID'] = df['gene_ID'].str.split(':').str[1]
        
        self.annotDF = df[ (df["biotype"] == "gene" )|
                           (df["biotype"] == "pseudogene") | 
                           (df["biotype"] == "three_prime_UTR") |
                           (df["biotype"] == "five_prime_UTR") |
                           (df["biotype"] == "ncRNA_gene") ]  
        return self.annotDF

####### TESTING #################################################################
        # print(df.iloc[0:5, 1:6]) # 0:5 is rows, 1:6 is columns
#                    gene_ID                 biotype start_index end_index strand
# 0  ID=gene:ENSG00000236875              pseudogene       12134     13783      +
# 1  ID=gene:ENSG00000236875  pseudogenic_transcript       12134     13783      +
# 2  ID=gene:ENSG00000236875                    exon       12134     12190      +
# 3  ID=gene:ENSG00000236875                    exon       12291     12340      +
# 4  ID=gene:ENSG00000236875                    exon       12726     12834      +


    def findGCRepeats(self, seq):
        for i in range(0, len(seq)):
            content = self.gcContent(seq[i:i+self.windowSize])
            return float(content)


    def gcContent(self, seq): 
        G_Count=0
        C_Count=0
        for i in range(0, len(seq), 1):
            if seq[i]=='G':
                G_Count=G_Count+1
            elif seq[i]=='C':
                C_Count=C_Count+1
        GC_Content=(G_Count + C_Count)/len(seq)
        float(GC_Content)
        GC_Content = GC_Content * 100 ## converts to %
        GC_Content = "{:.2f}".format(GC_Content)
        return(GC_Content)


    # sliding window code
    def slidingWindow(self):
        # sliding window on the forward strand
        df = pd.DataFrame(columns = ["start", "end", "CG-content", "strand"])
        
        for i in range(0,len(self.senseSeq) - self.windowSize + 1,1000):
            start = i
            end = i + self.windowSize
            SI = self.senseSeq[start:end+1] # get sequence between the indexes
            AI = self.antisenseSeq[start:end+1]
            
            # if ok, add to the data frame
            if self.findGCRepeats(SI) >= 55.0:
                # print("start=", start, "end=", end, "content=", self.findGCRepeats(SI), "+")
                    df.loc[len(df)] = [start, end, self.findGCRepeats(SI), "+"]
            if self.findGCRepeats(AI) >= 55.0:
                    df.loc[len(df)] = [start, end, self.findGCRepeats(AI), "-"]
        
        self.islandsDF = df
        
        # count islands in promoters and genes
        for i, row in self.islandsDF.iterrows():
            # row is now a series, subset like a list
            if "+" in row: self.isl_gene
            if "-" in row: self.isl_promoter
        
        return self.islandsDF
    
####### TESTING ###############################
        
#           start        end  CG-content strand
# 0         10000      13000       57.97      +
# 1         10000      13000       57.97      -
# 2         12000      15000       57.93      +
# 3         12000      15000       57.93      -
# 4         14000      17000       58.53      +
#         ...        ...         ...    ...
# 6013  138278000  138281000       62.67      -
# 6014  138280000  138283000       58.60      +
# 6015  138280000  138283000       58.60      -
# 6016  138282000  138285000       57.87      +
# 6017  138282000  138285000       57.87      -


    def CGcounts(self):
        # this is for sense strand, so for antisense it will be opposite C/G counts
        pairs_sense = pairs_anti = 0
        #check sense strand for pairs
        for i in range(0,len(self.senseSeq)-1,1):
                if self.senseSeq[i] == "C" and self.senseSeq[i+1] == "G": pairs_sense += 1
        #check antisense strand for pairs
        for j in range(0,len(self.antisenseSeq)-1,1):
                if self.antisenseSeq[j] == "C" and self.antisenseSeq[j+1] == "G": pairs_anti += 1
        
        return pairs_sense, pairs_anti


    def getSummaryStats(self):
        df = pd.DataFrame(columns = ["num CG pairs (sense)", "num CG pairs (antisense)", "total num islands", "avg island length"])
        
        
        df.loc[len(df)] = [self.CGcounts()[0], self.CGcounts()[1], len(self.islandsDF), self.windowSize]
        
        self.summaryDF = df
        return self.summaryDF


    def overlaps(self):
        
        # break into sense and antisense
        # for each unit ID in self.annotDF_sense and each island_sense, assign islands to list and put total num into df with that gene

        # islands_sense = pd.DataFrame()
        # islands_anti = pd.DataFrame()
        # annot_sense = pd.DataFrame()
        # annot_anti = pd.DataFrame()
                
        # islands_sense = self.islandsDF.loc[self.islandsDF.columns['strand'] == "+"] # this is a dataframe of only positive islands
        # islands_anti = self.islandsDF.loc[self.islandsDF.columns['strand'] == "-"] 
        # annot_sense = self.annotDF.loc[self.annotDF.columns['strand'] == "+"] # this is a dataframe of only positive genes/promoters
        # annot_anti = self.annotDF.loc[self.annotDF.columns['strand'] == "-"]
        
        # new df
        overlaps = pd.DataFrame(columns=["geneID", "start", "end", "strand", "island start", "island end"])
        
        # for index, island_row in islands_sense.iterrows():
        #     for index2, annot_row in annot_sense.iterrows():
        #         if (island_row["start"] >= annot_row["start"]) and (island_row["start"] <= annot_row["end"]):
        #             # add it
        #             overlaps.loc[len(overlaps)] = [annot_sense[index2]["biotype"], annot_sense[index2]["gene_ID"], annot_sense[index2]["strand"], islands_sense[index]["start"], islands_sense[index]["end"]]
           
        # for index, island_row in islands_anti.iterrows():
        #     for index2, annot_row in annot_anti.iterrows():
        #         if (island_row["start"] >= annot_row["start"]) and (island_row["start"] <= annot_row["end"]):
        #             overlaps.loc[len(overlaps)] = [annot_anti[index2]["biotype"], annot_anti[index2]["gene_ID"], annot_anti[index2]["strand"], islands_anti[index]["start"], islands_anti[index]["end"]]
        self.annotDF = self.annotDF.iloc[: , :-1]
        
        self.islandsDF['ints'] = [i for i in range(len(self.islandsDF))]
        isl_dict = self.islandsDF.set_index("ints").T.to_dict("list")
        annot_dict = self.annotDF.set_index("gene_ID").T.to_dict("list")
        
        # annotDF = columns=['subunit_ID', 'gene_ID', 'biotype', 'start_index', 'end_index', 'strand', 'desc'])
        # islands = df = pd.DataFrame(columns = ["start", "end", "CG-content", "strand"])
        for key, val in isl_dict.items():
            for key2, val2 in annot_dict.items():
                i_start = int(val[0])
                i_end = int(val[1])
                a_start = int(val2[2])
                a_end = int(val2[3])
                # island contained in the region
                if(i_start >= a_start) and (i_end <= a_end):
                    overlaps.loc[len(overlaps)] = [key2, val2[2], a_start, a_end, i_start, i_end]
                # region is contained in the island
                elif(i_start <= a_start) and (i_end >= a_end):
                    overlaps.loc[len(overlaps)] = [key2, val2[2], a_start, a_end, i_start, i_end]
                # island starts in the region and extends past the region
                elif(i_start >= a_start) and (i_start <= a_end):
                    overlaps.loc[len(overlaps)] = [key2, val2[2], a_start, a_end, i_start, i_end]
                # island ends inside the region
                elif(i_end >= a_start) and (i_end <= a_end):
                    overlaps.loc[len(overlaps)] = [key2, val2[2], a_start, a_end, i_start, i_end]
             
        
        # query this table by geneID to only show the input gene in the html
        # print("Overlaps", overlaps)
        self.overlapsDF = overlaps
        print(len(overlaps.index))
        # num_genes = num_promoters = 0
        # for k, row in overlaps.iterrows():
        #     if row["unitID"].startswith("I", 0,1):
        #         num_genes += 1
        #     else: num_promoters +=1
        
        # self.isl_genes = num_genes
        # self.isl_promoters = num_promoters
        

# ############### MAIN ################################################
def main():
    
    chrom9_genome = genome("Homo_sapiens.GRCh38.dna.chromosome.9.fa","Homo_sapiens.GRCh38.104.chromosome.9.gff3")
    chrom9_genome.readFASTA() # sets self.defLine and self.sequence
    chrom9_genome.getAntisenseSeq()
    chrom9_genome.readGFF3()  # sets annotDict
    annotDF = chrom9_genome.createAnnotDF() # sets annotDF
    
    # apply the sliding window
    islands_found = chrom9_genome.slidingWindow()
    
    # get counts
    pairs_sense,pairs_anti = chrom9_genome.CGcounts()
    # print("The sense strand has", pairs_sense, "CG pairs.")
    # print("The antisense strand has", pairs_anti, "CG pairs.")
    
    # get number of islands
    numIslands = chrom9_genome.islandsDF[chrom9_genome.islandsDF.columns[0]].count()
    # print("There were", numIslands, "islands detected in this chromosome using a sliding window of", chrom9_genome.windowSize, "and a step of 1000")
    # print("Since a sliding window was used, the length of each is the window size,", chrom9_genome.windowSize)
    
    # get the csv for the gene/island overlaps for the db
    chrom9_genome.overlaps()
    chrom9_genome.overlapsDF.to_csv("overlaps.csv", sep='\t', index=False)
    
    # get the csv's for SQL database
    islands_found.to_csv("islands.csv", sep='\t', index=False)
    annotDF.to_csv("annotation.csv", sep='\t', index=False)
    # get table for summary to put on website
    chrom9_genome.getSummaryStats()
    chrom9_genome.summaryDF.to_csv("summary.csv", sep='\t', index=False)


if __name__=='__main__':
    main()
    
    
    
    

