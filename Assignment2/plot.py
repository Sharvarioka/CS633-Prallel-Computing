#!/usr/bin/env python
# coding: utf-8

# In[47]:

import csv
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

sns.set()

demo_input_format = pd.DataFrame.from_dict({
    "D": [],
    "P": [],
    "ppn": [],
    "mode": [],  # 1 --> optimized, 0 --> standard
    "time": [],
})

file_name=["plot_Bcast.csv","plot_Gather.csv","plot_Reduce.csv"]
for name in file_name:
    temp=pd.read_csv(name,usecols=['Type','Time'])
    i=0
    for execution in range(10):
        for P in [4, 16]:
            for ppn in [1, 8]:
                for D in [16, 256, 2048]:
                    demo_input_format = demo_input_format.append({
                        "D": D, "P": P, "ppn": ppn, "mode": 'Default', "time": temp['Time'][i]
                    }, ignore_index=True)
                    i+=1
                    demo_input_format = demo_input_format.append({
                        "D": D, "P": P, "ppn": ppn, "mode": 'Optimized', "time": temp['Time'][i]
                    }, ignore_index=True)
                    i+=1

                    
    demo_input_format["(P, ppn)"] = list(map(lambda x, y: ("(" + x + ", " + y + ")"), map(str, demo_input_format["P"]), map(str, demo_input_format["ppn"])))

    print(demo_input_format)

    tkt=sns.catplot(x="(P, ppn)", y="time", data=demo_input_format, kind="bar", col="D", hue="mode")
    
    tkt.set_xlabels('P and ppn', fontsize=10) 
    tkt.set_ylabels('Time in seconds', fontsize=10)
    plt.savefig(name.strip(".csv")+".jpg")
    plt.show()




