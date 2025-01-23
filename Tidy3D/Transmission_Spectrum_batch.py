import numpy as np
import tidy3d as td 
from tidy3d.constants import C_0 
from tidy3d import web
import matplotlib.pyplot as plt 

### Geometrical parameters 
a = 1       # Period (micrometer)
freq_range_dimless = (0.24,0.26) # Dimensionless unit fa/c

batch_051_055 = web.Batch.from_file("delta_0.51_0.55/batch_data.json")
batch_data_051_055 = batch_051_055.load(path_dir = "data_batch_051_055")
sim_data_051 = batch_data_051_055["sim_0"]
sim_data_052 = batch_data_051_055["sim_1"]
sim_data_053 = batch_data_051_055["sim_2"]
sim_data_054 = batch_data_051_055["sim_3"]
sim_data_055 = batch_data_051_055["sim_4"]

fig,ax = plt.subplots(1,figsize=(10,8),tight_layout=True)
plt.show()