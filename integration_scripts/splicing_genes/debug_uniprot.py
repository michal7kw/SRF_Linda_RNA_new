import requests
import json

def debug_uniprot(go_id):
    url = "https://rest.uniprot.org/uniprotkb/search"
    
    queries = [
        f'go:"{go_id.replace("GO:", "")}" AND organism_id:9606'
    ]
    
    for query in queries:
        print(f"\nTesting Query: {query}")
        params = {
            'query': query,
            'fields': 'gene_primary',
            'format': 'json',
            'size': 1
        }
        
        try:
            response = requests.get(url, params=params, timeout=15)
            response.raise_for_status()
            data = response.json()
            count = len(data.get('results', []))
            print(f"Results count: {count}")
            if count > 0:
                print("SUCCESS!")
                print(json.dumps(data['results'][0], indent=2))
                break
        except Exception as e:
            print(f"Error: {e}")

if __name__ == "__main__":
    debug_uniprot("GO:0016180")
