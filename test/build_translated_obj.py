# Chamonix test : x0 = 1000 y0 = 6536

import os
import sys

if __name__ == '__main__':

    ply_dir = sys.argv[1]
    obj_dir = sys.argv[2]
    x0 = float(sys.argv[3])
    y0 = float(sys.argv[4])
    if (sys.argv[5]):
        delta = float(sys.argv[5])


    if not os.path.exists(obj_dir):
        print(f"Creating directory : {obj_dir}")
        os.makedirs(obj_dir)

    for ply_base in os.listdir(ply_dir):
        if ("final.ply" not in ply_base):
             continue
        x = float(ply_base[0:4])
        y = float(ply_base[5:9])
        if (sys.argv[5] and y - x != delta):
            continue
        ply_name = os.path.join(ply_dir, ply_base)
        obj_name = os.path.join(obj_dir, ply_base.replace("final.ply", "obj"))
        if (os.path.exists(obj_name)):
            continue
        cmd = f"./build/release/src/ply_to_translated_obj {ply_name} {obj_name} " \
              f"{x - x0} {y - y0}"
        print(cmd)
        os.system(cmd)






     

