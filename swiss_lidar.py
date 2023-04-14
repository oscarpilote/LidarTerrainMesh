import os
import sys
import requests
from io import BytesIO
from zipfile import ZipFile

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

if __name__ == '__main__':

    x0 = int(sys.argv[1])
    y0 = int(sys.argv[2])
    build_data_dicts()
    
    if ((x0,y0) not in surf_years and (x0,y0) not in alti_years):
        print("No data for this location.")
        sys.exit()
    
    base_dir = sys.argv[3]
    out_dir  = sys.argv[4]

    s = requests.Session()

    # Get data
    for x in (x0, x0 - 1, x0 + 1):
        for y in (y0, y0 - 1, y0 + 1):
            
            # LAS
            fname = os.path.join(base_dir, las_name(x, y))
            if (not os.path.exists(fname)):
                url = las_request(x, y)
                if (url):
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
                        printf("Could not get : ", url)
                else:
                    print("Skipping non available ", las_name(x,y))
            else:
                print("Recycling ", las_name(x, y))
            
            # Raster
            fname = os.path.join(base_dir, raster_name(x, y))
            if (not os.path.exists(fname)):
                url = raster_request(x, y)
                if (url):
                    print("Downloading ", raster_name(x,y))
                    r = s.get(url)
                    if (r.status_code == requests.codes.ok):
                        f = open(fname, "wb")
                        f.write(r.content)
                        f.close()
                    else:
                        printf("Could not get : ", url)
                else:
                    print("Skipping non available ", raster_name(x,y))
            else:
                print("Recycling ", raster_name(x, y))

    # Launch reconstruction
    cmd = f"./build/release/src/swiss_lidar {x0} {y0} {base_dir} {out_dir} " \
            "10 4 0.1 5"
    print(cmd);
    os.system(cmd)





     
