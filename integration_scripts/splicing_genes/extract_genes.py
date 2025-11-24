import requests
import re
import csv
import time
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Input list
TERMS = [
    "snRNA Processing (GO:0016180)",
    "Transport Of Mature mRNA Derived From An Intronless Transcript R-HSA-159231",
    "Positive Regulation Of RNA Splicing (GO:0033120)",
    "Regulation Of Alternative mRNA Splicing, Via Spliceosome (GO:0000381)",
    "Alternative mRNA Splicing, Via Spliceosome (GO:0000380)",
    "Processing Of Intronless Pre-mRNAs R-HSA-77595",
    "Nuclear-Transcribed mRNA Catabolic Process, Exonucleolytic (GO:0000291)",
    "Positive Regulation Of mRNA Splicing, Via Spliceosome (GO:0048026)",
    "Positive Regulation Of mRNA Processing (GO:0050685)",
    "mRNA Splice Site Recognition (GO:0006376)",
    "mRNA Cis Splicing, Via Spliceosome (GO:0045292)",
    "Regulation Of mRNA Processing (GO:0050684)",
    "Regulation Of mRNA Splicing, Via Spliceosome (GO:0048024)",
    "Transport Of Mature mRNA Derived From An Intron-Containing Transcript R-HSA-159236",
    "Regulation Of RNA Splicing (GO:0043484)",
    "Regulation Of rRNA Processing (GO:2000232)",
    "Spliceosomal Complex Assembly (GO:0000245)",
    "RNA Processing (GO:0006396)",
    "Processing Of Capped Intron-Containing Pre-mRNA R-HSA-72203",
    "Processing Of Capped Intronless Pre-mRNA R-HSA-75067",
    "mRNA Splicing R-HSA-72172",
    "mRNA Splicing - Major Pathway R-HSA-72163",
    "tRNA Processing (GO:0008033)",
    "mRNA Splicing - Minor Pathway R-HSA-72165",
    "rRNA Processing (GO:0006364)"
]

REGEX_GO = r"(GO:\d+)"
REGEX_REACTOME = r"(R-HSA-\d+)"

def create_session():
    """Creates a requests session with retry logic."""
    session = requests.Session()
    # Retry on 500 errors and connection issues
    retry = Retry(total=5, backoff_factor=0.5, status_forcelist=[500, 502, 503, 504])
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session

session = create_session()

def get_reactome_genes(reactome_id):
    """Fetches reference entities from Reactome."""
    url = f"https://reactome.org/ContentService/data/participants/{reactome_id}/referenceEntities"
    print(f"Fetching Reactome data for {reactome_id}...")
    
    try:
        response = session.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()
        
        genes = []
        for entity in data:
            if 'displayName' in entity:
                name = entity['displayName']
                # Reactome names are often "UniProt:ID GENE_NAME"
                # We split by space and try to take the second part if it exists
                parts = name.split(' ')
                if len(parts) > 1:
                    symbol = parts[1]
                else:
                    symbol = parts[0]
                
                genes.append({
                    'gene_symbol': symbol,
                    'source_id': reactome_id,
                    'database': 'Reactome'
                })
        print(f"  -> Found {len(genes)} genes for Reactome ID {reactome_id}")
        return genes
    except Exception as e:
        print(f"Error fetching Reactome {reactome_id}: {e}")
        return []

def get_uniprot_go_genes(go_id):
    """
    Fetches genes associated with a GO term from UniProt.
    Query: go:<id> AND organism_id:9606 (Human)
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    print(f"Fetching UniProt (GO) data for {go_id}...")
    
    # Search for the GO term in Humans (9606)
    # Note: UniProt API expects just the number for GO terms in the query
    clean_go_id = go_id.replace("GO:", "")
    query = f'go:"{clean_go_id}" AND organism_id:9606'
    
    all_genes = []
    next_cursor = None
    
    while True:
        params = {
            'query': query,
            'fields': 'gene_primary',
            'format': 'json',
            'size': 500  # Batch size
        }
        if next_cursor:
            params['cursor'] = next_cursor

        try:
            response = session.get(url, params=params, timeout=15)
            response.raise_for_status()
            data = response.json()
            
            results = data.get('results', [])
            for item in results:
                # Extract primary gene name
                # Safety check: ensure 'genes' list is not empty
                genes_list = item.get('genes', [])
                if not genes_list:
                    continue
                    
                primary_gene = genes_list[0].get('geneName', {}).get('value')
                if primary_gene:
                    all_genes.append({
                        'gene_symbol': primary_gene,
                        'source_id': go_id,
                        'database': 'GO (via UniProt)'
                    })
            
            # Pagination handling
            if 'nextCursor' in response.headers: # UniProt uses headers for cursor sometimes
                next_cursor = response.headers['nextCursor']
            elif 'next' in response.links: # Or link headers
                 # The python requests link header parsing is handy here, 
                 # but UniProt often returns the cursor in the body for JSON
                 pass
            
            # For JSON format, check if 'nextCursor' is in the body wrapper if applicable,
            # usually UniProt V2 uses Link headers for pagination in stream, 
            # but search endpoint returns 'nextCursor' in the JSON body sometimes?
            # Let's stick to the simplest method: Check the link header which is standard.
            
            # Note: The UniProt REST API pagination can be tricky. 
            # If the result count is massive, it suggests using /stream, but for GO terms /search is usually fine.
            # We will check if we got fewer results than requested size, implying end of list.
            if len(results) < 500:
                break
                
            # If we need to paginate, UniProt requires the cursor.
            # As of 2024, UniProtKB search uses link headers.
            if 'next' in response.links:
                next_url = response.links['next']['url']
                # Extract cursor from URL
                match = re.search(r'cursor=([^&]+)', next_url)
                if match:
                    next_cursor = match.group(1)
                else:
                    break
            else:
                break
                
        except Exception as e:
            print(f"Error fetching UniProt {go_id}: {e}")
            break
            
    print(f"  -> Found {len(all_genes)} genes for GO ID {go_id}")
    return all_genes

def main():
    collected_genes = []
    
    print(f"Processing {len(TERMS)} terms...")
    
    for term_string in TERMS:
        # 1. Try Reactome
        r_match = re.search(REGEX_REACTOME, term_string)
        if r_match:
            r_id = r_match.group(1)
            collected_genes.extend(get_reactome_genes(r_id))
            continue

        # 2. Try GO
        go_match = re.search(REGEX_GO, term_string)
        if go_match:
            go_id = go_match.group(1)
            collected_genes.extend(get_uniprot_go_genes(go_id))
            continue
            
        print(f"Could not parse ID from: {term_string}")

    # Deduplicate by Gene Symbol
    unique_map = {}
    for entry in collected_genes:
        symbol = entry['gene_symbol']
        # Keep the first occurrence, or overwrite if you prefer specific sources
        if symbol not in unique_map:
            unique_map[symbol] = entry
            
    # Write to CSV
    outfile = "extracted_genes_final.csv"
    with open(outfile, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['gene_symbol', 'source_id', 'database'])
        writer.writeheader()
        for symbol in sorted(unique_map.keys()):
            writer.writerow(unique_map[symbol])

    print(f"\nProcessing complete.")
    print(f"Found {len(unique_map)} unique genes.")
    print(f"Saved to {outfile}")

if __name__ == "__main__":
    main()