import requests
from bs4 import BeautifulSoup
import re
import wget
from pathlib import Path
import glob
import numpy as np
import xarray as xr

def find_netcdf_links(url):
    """
    Searches a web page for NetCDF file links and returns wget commands.

    Args:
        url (str): URL of the webpage to scrape.

    Returns:
        list: A list of wget commands for downloading the NetCDF files.
    """
    try:
        # Send a GET request to the URL
        response = requests.get(url)
        response.raise_for_status()  # Raise an error for bad HTTP status codes

        # Parse the page content using BeautifulSoup
        soup = BeautifulSoup(response.text, 'html.parser')

        # Find all links on the page
        links = soup.find_all('a', href=True)

        # Filter links for NetCDF files (e.g., .nc extension)
        netcdf_links = [link['href'] for link in links if link['href'].endswith('.nc')]

        # Convert relative URLs to absolute if needed
        absolute_links = [
            link if link.startswith('http') else requests.compat.urljoin(url, link)
            for link in netcdf_links
        ]

        return absolute_links

    except requests.RequestException as e:
        print(f"Error fetching the webpage: {e}")
        return []
    


if __name__ == "__main__":
    url = 'http://dsrs.atmos.umd.edu/DATA/soda3.15.2/REGRIDED/ocean/'
    out_dir = '/space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/soda/raw'

    Path(out_dir).mkdir(parents=True, exist_ok=True)
    downloaded = glob.glob(out_dir + '*.nc')

    wget_links_list = find_netcdf_links(url)

    years = np.unique([item.split('_')[-3] for item in wget_links_list])
    for year in years:
        year_links = [ls for ls in wget_links_list if year in ls]
        for dl_link in year_links:
            if '_mn_ocean' in dl_link:
                print(f'\n downloading {dl_link.split("/")[-1]} .... \n')
                try:
                    if (out_dir + dl_link.split('/')[-1]) not in downloaded:
                        wget.download(dl_link, out=out_dir + dl_link.split('/')[-1])
                    else:
                        print('Already downloaded.')
                except:
                    print('file not existant')
