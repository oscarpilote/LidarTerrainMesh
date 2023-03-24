import os
import sys
import requests
from io import BytesIO
from zipfile import ZipFile
from subprocess import Popen

alti_years = dict()
surf_years = dict()

def raster_name(x, y):
    return f"{x}_{y}.tif"

def raster_request(x, y):
    if ((x,y) not in alti_years):
        return None
    year = alti_years[(x,y)]
    return \
            "https://data.geo.admin.ch/ch.swisstopo.swissalti3d/swissalti3d_" \
            f'{year}_{x}-{y}/swissalti3d_{year}_{x}-{y}_2_2056_5728.tif'

def las_name(x, y):
    return f"{x}_{y}.las"

def las_request(x, y):
    if ((x,y) not in surf_years):
        return None
    year = surf_years[(x,y)]
    return \
            "https://data.geo.admin.ch/ch.swisstopo.swisssurface3d/" \
            f'swisssurface3d_{year}_{x}-{y}/swisssurface3d_{year}_{x}-{y}' \
            "_2056_5728.las.zip"

def build_data_dicts():
    f = open("swissalti3d.csv", "r")
    for line in f.readlines():
        tmp = line.split("_")
        year = int(tmp[1])
        x = int(tmp[2][0:4])
        y = int(tmp[2][5:9])
        alti_years[(x,y)] = year
    f.close()

    f = open("swisssurface3d.csv", "r")
    for line in f.readlines():
        tmp = line.split("_")
        year = int(tmp[1])
        x = int(tmp[2][0:4])
        y = int(tmp[2][5:9])
        surf_years[(x,y)] = year
    f.close()

def download_tile_data(s, x, y, verbose = True):
    # LAS
    fname = os.path.join(base_dir, las_name(x, y))
    if (not os.path.exists(fname)):
        url = las_request(x, y)
        if (url):
            if verbose:
                print("Downloading ", las_name(x,y))
            r = s.get(url)
            if (r.status_code == requests.codes.ok):
                zarch = ZipFile(BytesIO(r.content))
                fin = zarch.open(las_name(x,y))
                fout = open(fname, "wb")
                fout.write(fin.read())
                fin.close()
                fout.close()
            else:
                if verbose:
                    printf("Could not get : ", url)
        else:
            if verbose:
                print("Skipping non available ", las_name(x,y))
    else:
        if verbose:
            print("Recycling ", las_name(x, y))
    
    # Raster
    fname = os.path.join(base_dir, raster_name(x, y))
    if (not os.path.exists(fname)):
        url = raster_request(x, y)
        if (url):
            if verbose:
                print("Downloading ", raster_name(x,y))
            r = s.get(url)
            if (r.status_code == requests.codes.ok):
                f = open(fname, "wb")
                f.write(r.content)
                f.close()
            else:
                if verbose:
                    print("Could not get : ", url)
        else:
            if verbose:
                print("Skipping non available ", raster_name(x,y))
    else:
        if verbose:
            print("Recycling ", raster_name(x, y))

if __name__ == '__main__':

    xmin = int(sys.argv[1])
    xmax = int(sys.argv[2])
    ymin = int(sys.argv[3])
    ymax = int(sys.argv[4])
    
    build_data_dicts()
    base_dir = "data/"
    out_dir = "test/"

    s = requests.Session()
    
    # download bottom row and the one below
    print("Donwload bottom data")
    for y in (ymin - 1, ymin):
        for x in range(xmin - 1, xmax + 2):
            download_tile_data(s, x, y)

    # build rows
    for y in range(ymin, ymax + 1):
        print(f"Donwload row {y}")
        # download row above
        for x in range(xmin - 1, xmax + 2):
            download_tile_data(s, x, y + 1)
        # record commands
        commands = []
        for x in range(xmin, xmax + 1):
            commands.append(f"python3 swiss_lidar.py {x} {y} data/ test/ > {x}_{y}.log")
        # launch commands
        print("Launching commands")
        procs = [Popen(cmd, shell = True) for cmd in commands]
        # wait for completions
        for proc in procs:
            proc.wait()
        # remove unnecessary data
        print("Removing data")
        os.system(f"rm -f test/*_{y}.points.ply test/*_{y}.poisson.ply")
        os.system(f"rm -f data/*_{y - 1}.las data/*_{y - 1}.tif")
    # remove remaining data
    os.system(f"rm -f data/*_{ymax}.las data/*_{ymax}.tif")
    os.system(f"rm -f data/*_{ymax + 1}.las data/*_{ymax + 1}.tif")







     
