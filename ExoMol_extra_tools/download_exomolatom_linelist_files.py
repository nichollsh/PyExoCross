# Import all what we need
import os
import json
import re
import bz2
import shutil
import urllib3
import requests
import subprocess
from collections import defaultdict
from tqdm import tqdm
urllib3.disable_warnings()

# File Paths and Molecules
################## Could be changed ! ##################
# Select database: 'ExoMol' or 'ExoAtom'
database = 'ExoAtom'

# Directory that will hold the generated api__urls.txt file
url_dir = f'/scratch/p321409/opacity_lbl/exomol/url/'

# Full path to the urls file (derived from url_dir)
url_path = os.path.join(url_dir, 'api__urls.txt')
file_path = f'/scratch/p321409/opacity_lbl/{database.lower()}/'
# molecules = ['H2O', 'CO2', 'H2', 'H2S', 'N2', 'SiO']
# molecules = ['CH4', 'NH3', 'SO2', 'HCN', 'O3', 'N2O', 'O2']
molecules = ['Na']

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
    'Na': ['(23Na)',],
}

# Preferred datasets per atom (ExoAtom only). Example values:
preferred_datasets = {
    'Na': ['NIST', 'Kurucz'],
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

# Get ExoAtom URLs by parsing the ExoAtom search results page.
def get_exoatom_urls(species_list, preferred_datasets):
    base_url = "https://exomol.com"
    search_url = f"{base_url}/exoatom/"
    urls = []
    for species in tqdm(species_list):
        preferred = preferred_datasets.get(species)
        if not preferred:
            raise ValueError(
                f"Preferred datasets not provided for {species}. "
                "Populate preferred_datasets with desired dataset names."
            )
        preferred_set = set(preferred)
        response = requests.get(search_url, params={"qf": species})
        if response.status_code != 200:
            print(f"ExoAtom search error {response.status_code} for {species}")
            continue
        html = response.text
        matches = re.findall(r'href="(/exoatom/[^"]+)"', html)
        if not matches:
            print(f"Warning: no ExoAtom files found for {species}.")
            continue
        dataset_files = defaultdict(list)
        for path in matches:
            parts = path.strip("/").split("/")
            if len(parts) < 5:
                continue
            dataset = parts[3]
            if dataset not in preferred_set:
                continue
            dataset_files[dataset].append(base_url + path)
        for dataset, dataset_urls in dataset_files.items():
            print(f"{species} - {dataset}: {len(dataset_urls)} file(s)")
            for entry in dataset_urls:
                print(entry)
            urls.extend(dataset_urls)
        missing = preferred_set - set(dataset_files.keys())
        if missing:
            print(f"Warning: datasets not found for {species}: {sorted(missing)}")
    return urls

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
    if database.lower() == 'exoatom':
        urls = get_exoatom_urls(molecules, preferred_datasets)
    else:
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
    if database.lower() == 'exoatom':
        postprocess_exoatom_files(file_path, molecules)


def postprocess_exoatom_files(base_path, atoms):
    exoatom_dir = base_path
    exoatom_db_dir = os.path.join(exoatom_dir, 'db')
    
    if not os.path.isdir(exoatom_db_dir):
        print('No ExoAtom db folder found; skipping post-processing.')
        return

    for atom in atoms:
        src_root = os.path.join(exoatom_db_dir, atom)
        dst_root = os.path.join(exoatom_dir, atom)

        if not os.path.isdir(src_root):
            print(f'No ExoAtom files found to move for {atom}.')
            continue

        for root, dirs, files in os.walk(src_root):
            rel_root = os.path.relpath(root, src_root)
            target_root = os.path.join(dst_root, rel_root) if rel_root != '.' else dst_root
            os.makedirs(target_root, exist_ok=True)
            for filename in files:
                src_file = os.path.join(root, filename)
                dst_file = os.path.join(target_root, filename)
                if os.path.exists(dst_file):
                    os.remove(dst_file)
                shutil.move(src_file, dst_file)

        for root, dirs, files in os.walk(src_root, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))

        if os.path.isdir(src_root):
            os.rmdir(src_root)
            
        # compress_exoatom_files(dst_root)


def compress_exoatom_files(root_dir):
    for root, _dirs, files in os.walk(root_dir):
        for filename in files:
            if filename.endswith(('.states', '.trans')):
                src_path = os.path.join(root, filename)
                dst_path = src_path + '.bz2'
                if os.path.exists(dst_path):
                    continue
                with open(src_path, 'rb') as src, bz2.open(dst_path, 'wb') as dst:
                    shutil.copyfileobj(src, dst)
                os.remove(src_path)

download_files(molecules, url_path, preferred_isotopologues)