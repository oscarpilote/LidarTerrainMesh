import os
import sys
import requests

links = dict()

def raster_name(x, y):
    return f"{x:0>{4}}_{y:0>{4}}.tif"

def raster_request(x, y):
    return None

def las_name(x, y):
    return f"{x:0>{4}}_{y:0>{4}}.copc.laz"

def las_request(x, y):
    if ((x,y) not in links):
        return None
    return links[(x,y)]

def build_data_dicts():
    f = open("french_lidar.csv", "r")
    for line in f.readlines():
        line = line[:-1]
        tmp = line.split("FXX_")
        x = int(tmp[1][0:4])
        y = int(tmp[1][5:9])
        links[(x,y)] = line;
    f.close()

if __name__ == '__main__':

    x0 = int(sys.argv[1])
    y0 = int(sys.argv[2])
    build_data_dicts()
    
    if ((x0,y0) not in links):
        print("No data for this location.")
        sys.exit()
    
    base_dir = sys.argv[3]
    out_dir  = sys.argv[4]

    if not os.path.exists(base_dir):
        print(f"Creating directory : {base_dir}")
        os.makedirs(base_dir)
    if not os.path.exists(out_dir):
        print(f"Creating directory : {out_dir}")
        os.makedirs(out_dir)


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
                        fout = open(fname, "wb")
                        fout.write(r.content)
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
    if ((len(sys.argv) >= 6) and int(sys.argv[5]) == 11):
        cmd = f"./build/release/src/french_lidar {x0} {y0} " \
              f"{base_dir} {out_dir} 11 4 0.05 5"
    else:
        cmd = f"./build/release/src/french_lidar {x0} {y0} " \
              f"{base_dir} {out_dir} 10 4 0.1 5"

    print(cmd);
    os.system(cmd)





     
