import re
import csv
import requests
from Bio import ExPASy
from Bio import SwissProt
from requests.adapters import HTTPAdapter, Retry

# compile a regular expression pattern to match the next link header
re_next_link = re.compile(r'<(.+)>; rel="next"')

# configure the retry logic for failed requests
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])

# create a session object with the retry logic and mount it to the https protocol
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

# define a function to extract the URL of the next page of search results from the headers
def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

# define a function to fetch batches of search results until there are no more pages
def get_batch(batch_url):
    while batch_url:
        # make a GET request to the batch URL
        response = session.get(batch_url)
        # raise an exception if the request fails
        response.raise_for_status()
        # get the total number of search results from the response headers
        total = response.headers["x-total-results"]
        # yield the batch of search results and the total number of results
        yield response, total
        # get the URL of the next page of search results from the headers
        batch_url = get_next_link(response.headers)

# define a list of search URLs to fetch
# define a list of search URLs to fetch
urls = [
    'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Ccc_interaction&format=tsv&query=cell%20AND%20%28reviewed%3Atrue%29&size=500',
    'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Ccc_interaction&format=tsv&query=mitochondria%20AND%20%28reviewed%3Atrue%29&size=500',
    'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Ccc_interaction&format=tsv&query=nucleus%20AND%20%28reviewed%3Atrue%29&size=500',
    'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Ccc_interaction&format=tsv&query=endoplasmic%20AND%20%28reviewed%3Atrue%29&size=500',
    'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Ccc_interaction&format=tsv&query=golgi%20AND%20%28reviewed%3Atrue%29&size=500',
    'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Ccc_interaction&format=tsv&query=lysosome%20AND%20%28reviewed%3Atrue%29&size=500',
    'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Ccc_interaction&format=tsv&query=peroxisome%20AND%20%28reviewed%3Atrue%29&size=500',
    'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Ccc_interaction&format=tsv&query=cytosol%20AND%20%28reviewed%3Atrue%29&size=500'
]

progress = 0
# open the output file for writing
with open('Cell-interactions.tsv', 'w') as f:
    for url in urls:
        for batch, total in get_batch(url):
            lines = batch.text.splitlines()
            print(lines)
            # write the header row to the output file on the first batch
            if not progress:
                print(lines[0], file=f)
            # write the data rows to the output file
            for line in lines[1:]:
                print(line, file=f)
            # update the progress counter and print the progress to the console
            progress += len(lines[1:])
            print(f'{progress} / {total}')
            

# open the TSV file for reading
with open('Cell-interactions.tsv', 'r') as tsvfile:
    # create a CSV reader object for the TSV file
    reader = csv.reader(tsvfile, delimiter='\t')
    # open a text file for writing the output
    with open('output2.txt', 'w') as txtfile:
        seen_ids = set()  # to keep track of unique IDs
        # loop through each row of the TSV file
        for row in reader:
            for cell in row:
                # split the cell contents on semicolons to get a list of IDs
                ids = cell.split(';')
                # loop through each ID in the list
                for id in ids:
                    # check if the ID is a 6-character string and hasn't been seen before
                    if len(id) == 6 and id not in seen_ids:
                        # write the ID to the output file
                        txtfile.write(id + '\n')
                        # add the ID to the set of seen IDs
                        seen_ids.add(id)
                        
                        
                        

protein_counter = 0  # initialize a counter for the number of proteins processed

# open a TSV file for writing the protein information
with open("protein_info.tsv", "w", newline="") as tsv_file:
    writer = csv.writer(tsv_file, delimiter="\t")
    # write the header row to the TSV file
    writer.writerow(["Entry", "Subcellular location [CC]", "Sequence"])
    
    # open the output2.txt file for reading the list of protein IDs
    with open('output2.txt', 'r') as f:
        lines = f.readlines()
        
    # loop through each protein ID in the list
    for protein_id in lines:
        # extract the 6-character protein ID from the line
        protein_id = str(protein_id[0:6])
        print(protein_id)

        # Use ExPASy to connect to the Uniprot server and retrieve the data
        try:
            handle = ExPASy.get_sprot_raw(protein_id)
            record = SwissProt.read(handle)
            
            # loop through the comments in the protein record to find the subcellular location
            for comment in record.comments:
                if comment.startswith("SUBCELLULAR LOCATION:"):
                    subcell_loc = comment.replace("SUBCELLULAR LOCATION:", "").strip()
                    #print("Subcellular location:", subcell_loc)
                    break

            # write the protein information to the TSV file
            writer = csv.writer(tsv_file, delimiter="\t")
            writer.writerow([record.entry_name, subcell_loc, record.sequence])
        except Exception as e:
            print("exception")
        
        protein_counter = protein_counter + 1  # increment the protein counter
        print(protein_counter)  # print the protein counter to the console
