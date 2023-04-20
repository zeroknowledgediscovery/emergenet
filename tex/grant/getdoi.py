from habanero import Crossref
import bibtexparser
from tqdm import tqdm
import sys

cr = Crossref()

with open(sys.argv[1]) as bibtex_file:
    bib_database = bibtexparser.load(bibtex_file)

count=0
count_present=0
count_total=0
for entry in tqdm(bib_database.entries):
    count_total=count_total+1
    if 'doi' not in entry.keys():
        if 'title' in entry.keys():
            x= cr.works(query = entry['title'], limit = 1)
            if 'DOI' in x['message']['items'][0].keys():
                entry['doi']=x['message']['items'][0]['DOI'] 
                count=count+1
    else:
        count_present=count_present+1
Count={'found':count,'already_present':count_present,'missing':count_total-count-count_present}
print(Count)

with open(sys.argv[1].replace('.bib','new.bib'), 'w') as bibtex_file:
    bibtexparser.dump(bib_database, bibtex_file)
