# Import all what we need
import os
import re
import json
import urllib3
import requests
import subprocess
from tqdm import tqdm
urllib3.disable_warnings()

# File Paths and Molecules
################## Could be changed ! ##################
# Directory that will hold the generated api__urls.txt file
url_dir = '/scratch/p321409/opacity_lbl/exomol/url/'

# Full path to the urls file (derived from url_dir)
url_path = os.path.join(url_dir, 'api__urls.txt')
file_path = '/scratch/p321409/opacity_lbl/exomol/'
all_isotopologues = {
    # 'MgH': {
    #     '24Mg-1H': {'wn_range': None},
    #     '25Mg-1H': {'wn_range': None},
    # },
    # 'H2O': {
    #     '1H2-16O': {'wn_range': [41000, 41200]},
    # },
    'H2O': {'1H2-16O':  {'wn_range': None}},
    'H2':  {'1H2':      {'wn_range': None}},
    'CO2': {'12C-16O2': {'wn_range': None}},
    'CO':  {'12C-16O':  {'wn_range': None}},
    'O2':  {'16O2':     {'wn_range': None}},
    'CH4': {'12C-1H4':  {'wn_range': None}},
    'SO2': {'32S-16O2': {'wn_range': None}},
    'SH':  {'32S-1H':   {'wn_range': None}},
    'S2':  {'32S2':     {'wn_range': None}},
    'H2S': {'1H2-32S':  {'wn_range': None}},
    'N2':  {'14N2':     {'wn_range': None}},
    'N2O': {'14N2-16O': {'wn_range': None}},
    'NH3': {'14N-1H3':  {'wn_range': None}},
    'HCN': {'1H-12C-14N': {'wn_range': None}},
    'OCS': {'16O-12C-32S': {'wn_range': None}},
    'CN':  {'12C-14N': {'wn_range': None}},

    'SiH4': {'28Si-1H4': {'wn_range': None}},
    'SiO': {'28Si-16O': {'wn_range': None}},
    'FeO': {'56Fe-16O': {'wn_range': None}},
    'TiO': {'48Ti-16O': {'wn_range': None}},
    'MgO': {'24Mg-16O': {'wn_range': None}},
}

# get one molecule and its isotopologues
molec = 'O2'
molecule_isotopologues = dict([(molec, all_isotopologues.get(molec))])
print('molecule_isotopologues:', molecule_isotopologues)

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
def normalize_molecule_isotopologues(molecule_isotopologues):
    molecules = list(molecule_isotopologues.keys())
    isotopologue_configs = []
    for molecule in molecules:
        molecule_isos = molecule_isotopologues[molecule]
        if molecule_isos in (None, ''):
            isotopologue_configs.append(None)
        elif isinstance(molecule_isos, str):
            isotopologue_configs.append({molecule_isos: {'wn_range': None}})
        elif isinstance(molecule_isos, dict):
            isotopologue_config = {}
            for isotopologue, config in molecule_isos.items():
                if isinstance(config, dict):
                    isotopologue_config[isotopologue] = config
                else:
                    isotopologue_config[isotopologue] = {'wn_range': config}
            isotopologue_configs.append(isotopologue_config)
        else:
            isotopologue_configs.append(
                {isotopologue: {'wn_range': None} for isotopologue in molecule_isos}
            )
    return molecules, isotopologue_configs


def get_wn_range(isotopologue_config):
    if isotopologue_config is None:
        return None
    return isotopologue_config.get('wn_range')


def strict_states_filename(isotopologue, dataset):
    return f'{isotopologue}__{dataset}.states.bz2'


def strict_trans_filename(isotopologue, dataset):
    return f'{isotopologue}__{dataset}.trans.bz2'


def strict_segmented_trans_pattern(isotopologue, dataset):
    return re.compile(
        rf'^{re.escape(isotopologue)}__{re.escape(dataset)}__(\d+)-(\d+)\.trans\.bz2$'
    )


def infer_iso_slug_from_url(url, dataset):
    '''
    Infer the isotopologue slug from the given URL and dataset name.
    Returns the isotopologue slug if it can be inferred, otherwise returns None.
    '''
    filename = os.path.basename(url)
    suffixes = [
        f'__{dataset}.states.bz2',
        f'__{dataset}.trans.bz2',
    ]
    for suffix in suffixes:
        if filename.endswith(suffix):
            slug = filename[:-len(suffix)]
            # print(f"       slug from URL -> {slug}")
            return slug

    # Segmented-by-wavenumber trans file, e.g. 1H2-16O__POKAZATEL__00000-00100.trans.bz2
    segmented_match = re.match(rf'^(.+)__{re.escape(dataset)}__\d+-\d+\.trans\.bz2$', filename)
    if segmented_match is not None:
        slug = segmented_match.group(1)
        # print(f"       slug from URL -> {slug}")
        return slug

    # Split-by-label trans file, e.g. 16O2__SWYT__M1.trans.bz2 / 16O2__SWYT__E2.trans.bz2
    labeled_match = re.match(rf'^(.+)__{re.escape(dataset)}__[A-Za-z0-9]+\.trans\.bz2$', filename)
    if labeled_match is not None:
        slug = labeled_match.group(1)
        # print(f"       slug from URL -> {slug}")
        return slug

    return None


def trans_split_kind(url, isotopologue, dataset):
    '''
    Classify a transition file URL against the given isotopologue/dataset.

    Some ExoMol datasets split their transitions across several files that
    are NOT wavenumber-segmented, e.g. by transition type
    (16O2__SWYT__M1.trans.bz2, 16O2__SWYT__E2.trans.bz2). These "labeled"
    files must all be downloaded regardless of wn_range, since wn_range
    filtering only makes sense for wavenumber-segmented files.

    Returns a tuple (kind, info):
        ('single', None)       - the one-file whole-range trans file
        ('range', (min, max))  - a wavenumber-range-segmented file
        ('labeled', label)     - a non-wavenumber labeled split (e.g. 'M1')
        (None, None)           - does not match this isotopologue/dataset
    '''
    filename = os.path.basename(url)
    if filename == strict_trans_filename(isotopologue, dataset):
        return 'single', None

    range_match = strict_segmented_trans_pattern(isotopologue, dataset).match(filename)
    if range_match is not None:
        return 'range', (int(range_match.group(1)), int(range_match.group(2)))

    labeled_match = re.match(
        rf'^{re.escape(isotopologue)}__{re.escape(dataset)}__([A-Za-z0-9]+)\.trans\.bz2$',
        filename,
    )
    if labeled_match is not None:
        return 'labeled', labeled_match.group(1)

    return None, None


def trans_url_in_wn_range(url, isotopologue, dataset, wn_range):
    kind, info = trans_split_kind(url, isotopologue, dataset)
    if kind is None:
        return False
    if kind in ('single', 'labeled'):
        return True

    # kind == 'range'
    if wn_range in (None, []):
        return True
    wn_min, wn_max = wn_range
    file_min, file_max = info
    return file_min >= wn_min and file_max <= wn_max


def get_urls(molecule_isotopologues):
    """Get the download url from API."""
    molecules, isotopologue_configs = normalize_molecule_isotopologues(molecule_isotopologues)
    for molecule, isotopologue_config in zip(molecules, isotopologue_configs):
        if isotopologue_config is None:
            raise ValueError(
                f"Isotopologues not provided for {molecule}. "
                "Populate molecule_isotopologues with desired isotopologues."
            )
        print(f"Getting download URLs for {molecule} with isotopologues: {list(isotopologue_config.keys())}")
    api_url = get_api(molecules)
    urls = []
    for i in tqdm(range(len(molecules))):
        target_isotopologue_config = isotopologue_configs[i]
        if target_isotopologue_config is None:
            raise ValueError(
                f"\nIsotopologues not provided for {molecules[i]}. "
                "Populate molecule_isotopologues with desired isotopologues."
            )
        response = requests.get(api_url[i], timeout=60)
        if(response.status_code != 200):
            print('ExoMol API Error' + str(response.status_code))

        # If the obtained status code is 200, it is correct.
        else:
            print(f"\nSuccessfully retrieved API response for {molecules[i]}")
            content = response.text            # Get the relevant content.
            json_dict = json.loads(content)    # Convert json into dictionary.
            found_isotopologues = set()
            for iso_formula, iso_info in json_dict.items():
                linelist_info = iso_info.get('linelist', {})
                print(f"     isotopologue {iso_formula} for molecule {molecules[i]}")
                for dataset, files_info in linelist_info.items():
                    # skip to the dictionary item
                    if not isinstance(files_info, dict):
                        continue

                    # get recommended files for the dataset
                    if files_info.get('recommended'):
                        files_meta = files_info.get('files', [])
                        nfiles = len(files_meta)
                        print("       recommended dataset", dataset, "has", nfiles, "file(s)")
                        trans_count = 0
                        trans_urls = []
                        states_url = None
                        iso_slug = None
                        current_wn_range = None

                        # loop through files in the request
                        for j in range(nfiles):
                            # construct url
                            file_meta = files_meta[j]
                            url = "https://www." + file_meta.get('url')
                            filename = os.path.basename(url)

                            # print(f"           file {j+1}/{nfiles}: {file_meta.get('url')}")

                            # get the isotopologue slug from the URL or filename
                            inferred_iso_slug = infer_iso_slug_from_url(url, dataset)
                            if inferred_iso_slug is not None:
                                iso_slug = inferred_iso_slug
                            if iso_slug is None:
                                print(f"Warning: Could not infer isotopologue slug from URL {url}. Skipping.")
                                continue

                            if iso_slug not in target_isotopologue_config:
                                # Isotopologue is not in target isotopologue config for molecule
                                print("        skipping unrequested isotopologue")
                                continue

                            # get wavenumber range
                            current_wn_range = (get_wn_range(target_isotopologue_config[iso_slug]))

                            # get the states and trans files
                            if filename == strict_states_filename(iso_slug, dataset):
                                states_url = url
                                def_url = states_url.replace('.states.bz2','.def.json')
                                pf_url = states_url.replace('.states.bz2','.pf')
                            elif trans_url_in_wn_range(url, iso_slug, dataset, current_wn_range):
                                trans_urls.append(url)
                                trans_count += 1
                        
                        # we didn't pick up a states file for this isotopologue and dataset
                        if states_url is None:
                            # Not one of the requested isotopologues (or no slug could be
                            # inferred at all) - don't fall back to guessing a states URL,
                            # otherwise every isotopologue in the API response gets pulled in.
                            if iso_slug is None or iso_slug not in target_isotopologue_config:
                                continue

                            print(f'{molecules[i]} - {iso_slug} - {dataset}: no strict states file found.')

                            # try looking manually...
                            states_url = f"https://exomol.com/db/{molecules[i]}/{iso_slug}/{dataset}/{iso_slug}__{dataset}.states.bz2"
                            if requests.head(states_url).status_code == 200:
                                def_url = states_url.replace('.states.bz2','.def.json')
                                pf_url = states_url.replace('.states.bz2','.pf')
                                print(f"       Found states file manually: {states_url}")
                            else:
                                # no luck...
                                continue

                        # we didn't pick up any trans files for this isotopologue and dataset

                        # record found isotopologue in the set
                        found_isotopologues.add(iso_slug)
                        start = len(urls)
                        urls.extend([def_url, pf_url, states_url])
                        urls.extend(trans_urls)
                        print(f'        {molecules[i]} - {iso_slug} - {dataset}: {trans_count} file(s) to download')
                        for entry in urls[start:]:
                            print(f"       {entry}")
            
            if target_isotopologue_config is not None:
                missing_isotopologues = [
                    iso for iso in target_isotopologue_config
                    if iso not in found_isotopologues
                ]
                for missing_isotopologue in missing_isotopologues:
                    print(f'{molecules[i]} - {missing_isotopologue}: recommended line list not found in ExoMol API response.')
               
    return(urls) 

# Download line list Files
# We write all the download URLs into a text file, name it as api__urls.txt. 
# In Linux, we use command:
# wget  -r -nH --cut-dirs=1 -P savePath -i PathOFapi__urls.txt
# Download line list files with urls and save them into correspoding folders.
def download_files(molecule_isotopologues, url_path):

    # get urls for this linelist
    urls = get_urls(molecule_isotopologues)

    # check total size of files to download
    total_size = 0
    for url in urls:
        response = requests.head(url)
        if response.status_code == 200:
            total_size += int(response.headers.get('Content-Length', 0))
    total_size_gb = total_size / (1024 ** 3)
    print(f"\nTotal size of {len(urls)} files to download: {total_size_gb:.2f} GB")

    # Save all URLs to a text file
    os.makedirs(os.path.dirname(url_path), exist_ok=True)
    with open(url_path, "w", encoding="utf-8") as fh:
        for entry in urls:
            fh.write(f"{entry}\n")
    print('\nTargeted URLs have been saved to', url_path)

    # Download all files using wget command
    command = f'wget -r  -N -nH --show-progress --cut-dirs=1 -P {file_path} -i {url_path}'
    subprocess.run(command, shell=True)

    # Print success message
    print('\nAll files have been downloaded to', file_path, 'folder!')

if __name__ == '__main__':
    download_files(molecule_isotopologues, url_path)
