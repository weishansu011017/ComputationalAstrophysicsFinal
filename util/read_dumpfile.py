import numpy as np
import pandas as pd
import h5py

class ParticlesTable:
    def __init__(self):
        self.Table = None
        self.params = {}

    @staticmethod
    def read_dumpfile(filepath):
        obj = ParticlesTable()
        with h5py.File(filepath, "r") as f:
            # Table data
            particle_index = f["Table/particle_index"][:]
            x = f["Table/x"][:]
            y = f["Table/y"][:]
            z = f["Table/z"][:]
            vx = f["Table/vx"][:]
            vy = f["Table/vy"][:]
            vz = f["Table/vz"][:]
            m = f["Table/m"][:]
            h = f["Table/h"][:]
            dt = f["Table/dt"][:]
            obj.Table = pd.DataFrame({
                "particle_index": particle_index,
                "x": np.float64(x),
                "y": np.float64(y),
                "z": np.float64(z),
                "vx": np.float64(vx),
                "vy": np.float64(vy),
                "vz": np.float64(vz),
                "m": np.float64(m),
                "h": np.float64(h),
                "dt": np.float64(dt),
            })
            # params data
            param_group = f["params"]
            for k in param_group.keys():
                val = param_group[k][()]
                obj.params[k] = val

        return obj
    
    def transfer_cgs(self, year = True):
        udist = np.float64(self.params["udist"])
        umass = np.float64(self.params["umass"])
        utime = np.float64(self.params["utime"])
        uv = udist/utime

        self.Table["x"] *= udist
        self.Table["y"] *= udist
        self.Table["z"] *= udist
        self.Table["vx"] *= uv
        self.Table["vy"] *= uv
        self.Table["vz"] *= uv
        self.Table["m"] *= umass
        self.Table["h"] *= udist
        self.Table["dt"] *= utime

        self.params["Mtot"] *= umass
        self.params["t"] *= utime
        if year:
            self.params["t"] /= 31556926 
            self.Table["dt"] /= 31556926 

