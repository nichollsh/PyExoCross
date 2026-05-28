# Import all what we need
import os
import json
import urllib3
import requests
import subprocess
from tqdm import tqdm
urllib3.disable_warnings()

# File Paths and Molecules
################## Could be changed ! ##################
# Directory that will hold the generated api__urls.txt file
url_dir = '/scratch/p321409/opacity/lbl/exomol/url/'
# Full path to the urls file (derived from url_dir)
url_path = os.path.join(url_dir, 'api__urls.txt')
file_path = '/scratch/p321409/opacity/lbl/exomol/'
# molecules = ['H2O', 'CO2', 'H2', 'H2S', 'N2', 'SiO']
# molecules = ['CH4', 'NH3', 'SO2', 'HCN', 'O3', 'N2O', 'O2']
molecules = ['CO']

# Preferred isotopologues per molecule (must match ExoMol API keys).
# Example values (replace with desired isotopologues):
preferred_isotopologues = {
    'H2O': ['(1H)2(16O)',],
    'CO2': ['(12C)(16O)2',],
    'H2': ['(1H)2',],
    'H2S': ['(1H)2(32S)',],
    'N2': ['(14N)2',],
    'SiO': ['(28Si)(16O)',],
    'CO': ['(12C)(16O)',],
    'CH4': ['(12C)(1H)4',],
    'NH3': ['(14N)(1H)3',],
    'SO2': ['(32S)(16O)2',],
    'HCN': ['(1H)(12C)(14N)',],
    'O3': ['(16O)3',],
    'N2O': ['(14N)2(16O)',],
    'O2': ['(16O)2',],
    'SiO2': ['(28Si)(16O)2',],
    'MgO': ['(24Mg)(16O)',],
}
########################################################

# Get API URLs
def get_api(molecules):
    molecule_str = []
    api_url = []
    for i in range(len(molecules)):
        molecule_str.append(molecules[i].replace('_p','+').split('__')[0].replace('+','_p'))
        api_url.append('https://exomol.com/api/?molecule=*&datatype=linelist'.replace('*',molecule_str[i]))
    return(api_url)

# Get Download Links with API
def get_urls(molecules, preferred_isotopologues):
    """Get the download url from API for preferred isotopologues only."""
    api_url = get_api(molecules)
    urls = []
    for i in tqdm(range(len(molecules))):
        preferred = preferred_isotopologues.get(molecules[i])
        if not preferred:
            raise ValueError(
                f"Preferred isotopologues not provided for {molecules[i]}. "
                "Populate preferred_isotopologues with desired API keys."
            )
        preferred_set = set(preferred)
        response = requests.get(api_url[i])
        if(response.status_code != 200):
            print('ExoMol API Error' + str(response.status_code))

        # If the obtained status code is 200, it is correct.
        else:
            content = response.text            # Get the relevant content.
            json_dict = json.loads(content)    # Convert json into dictionary.
            iso_formulas = list(json_dict.keys())
            found_preferred = set()
            for iso_formula in iso_formulas:
                if iso_formula not in preferred_set:
                    continue
                found_preferred.add(iso_formula)
                datasets = list(json_dict[iso_formula]['linelist'].keys())[1:]
                for dataset in datasets:
                    files_info = json_dict[iso_formula]['linelist'][dataset]
                    if files_info['recommended'] == True:
                        files_meta = files_info['files']
                        nfiles = len(files_meta)
                        trans_count = 0
                        trans_urls = []
                        for j in range(nfiles):
                            file_meta = files_meta[j]
                            url = "https://www." + file_meta.get('url')
                            if url.endswith('states.bz2'):
                                states_url = url.replace('_v1','')
                                def_url = states_url.replace('.states.bz2','.def.json')
                                pf_url = states_url.replace('.states.bz2','.pf')
                            elif url.endswith('trans.bz2'):
                                trans_urls.append(url)
                                trans_count += 1
                            else:
                                print('No line list files.')
                        start = len(urls)
                        urls.extend([def_url, pf_url, states_url])
                        urls.extend(trans_urls)
                        print(f'{molecules[i]} - {iso_formula} - {dataset}: {trans_count} trans file(s)')
                        for entry in urls[start:]:
                            print(entry)
            missing = preferred_set - found_preferred
            if missing:
                print(f"Warning: isotopologues not found for {molecules[i]}: {sorted(missing)}")
                        
    return(urls) 

# Download line list Files
# We write all the download URLs into a text file, name it as api__urls.txt. 
# In Linux, we use command:
# wget  -r -nH --cut-dirs=1 -P savePath -i PathOFapi__urls.txt
# Download line list files with urls and save them into correspoding folders.
def download_files(molecules, url_path, preferred_isotopologues):
    urls = get_urls(molecules, preferred_isotopologues)
    # Save all URLs to a text file
    os.makedirs(os.path.dirname(url_path), exist_ok=True)
    with open(url_path, "w", encoding="utf-8") as fh:
        for entry in urls:
            fh.write(f"{entry}\n")
    print('\nAll URLs have been saved to', url_path)
    command = f'wget -r -nH --cut-dirs=1 -P {file_path} -i {url_path}'
    subprocess.run(command, shell=True)
    print('\nAll files have been downloaded to', file_path, 'folder!')

download_files(molecules, url_path, preferred_isotopologues)