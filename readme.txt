Steps to use program

Must have Homo_sapiens.GRCh38.dna.chromosome.9.fa & Homo_sapiens.GRCh38.104.chromosome.9.gff3 in your PWD


From python IDE:
Run the entire python script to render files
- CSV files are created in same folder as script called islands.csv annotation.csv overlaps.csv and summary.csv

Run the createDB.py to create the tables
- I don't think this will work exactly as planned
- when importing tables into db manually, we chose to delete the first row of each, and we need to drop the "end" column from overlaps and rename "strand" as "gene_end".

In MySQL:
Navigate to bio466-f15.csi.miamioh.edu/kockale/cgi.test_copy.py into chrome browser
