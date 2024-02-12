import numpy as np
import h5py as h5
import os
import sys

def run(subID):
    name = 'subhalo' + str(subID)
    sub_dir = name
    
    # load in the tracer ids previously selected
    snap_list = np.arange(100)
    files = []
    for i in snap_list:
        f = h5.File(sub_dir + '/props_' + name + '_snap' + str(i).zfill(3) + '.h5', mode='r')
        files.append(f)
    
    keys = list(files[0].keys())
    out = {}
    for key in keys:
        out[key] = []
    
    header_keys = list(files[0]['Header'].attrs.keys())
    header_out = {}
    for hkey in header_keys:
        header_out[hkey] = []
    
    for f in files:
        for key in keys:
            if key == 'Header':
                for hkey in header_keys:
                    header_out[hkey].append(f['Header'].attrs[hkey])
            else:
                out[key].append(f[key][:])
    
    for key in keys:
        out[key] = np.array(out[key])
    
    for hkey in header_keys:
        header_out[hkey] = np.array(header_out[hkey])
        
    for f in files:
        f.close()
        
    fout = h5.File(sub_dir + '/props_' + name + '.h5', mode='w')
    for key in out.keys():
        fout.create_dataset(key, data=out[key])
    for hkey in header_out.keys():
        fout.create_dataset(hkey, data=header_out[hkey])
    fout.close()
    
    return None
    
if __name__ == '__main__':
    subID = int(sys.argv[1])

    name = 'subhalo' + str(subID)
    assert os.path.exists(name)
    
    run(subID)
