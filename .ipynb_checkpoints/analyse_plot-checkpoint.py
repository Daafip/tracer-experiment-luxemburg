import sys
import subprocess
import time

# # implement pip as a subprocess:
# packages = ["pandas","matplotlib"]
# for package in packages:
#     subprocess.check_call([sys.executable, '-m', 'pip', 'install', package])
time.sleep(1)
import pandas as pd
import matplotlib.pyplot as plt
import glob
import scipy

solution = 0.1 # g/l
solution_ug = solution * 10**6 # ug/l
solution_ug_ml = solution_ug / 10**3 # ug/ml

files = glob.glob("*.tab")

print(f"found: {files}")

def plot(index):
    df = pd.read_csv(files[index],delimiter="\t",skiprows=1,index_col=0, parse_dates=True)
    fig, ax = plt.subplots(1,2, figsize=(12,5))
    zero = df['Rhodamine (ug/L)'].iloc[:20].mean()
    df['Rhodamine (ug/L)']  = df['Rhodamine (ug/L)'] - zero
    df[['Rhodamine (ug/L)']].plot(ax=ax[0])
    df[['Temp (C)']].plot(ax=ax[1])

    title = files[index][:files[index].find(".tab")]
    
    ax[0].set_ylabel('Rhodamine (ug/L)')
    
    ax[1].set_ylabel('Temp (C)')
    
    area = scipy.integrate.simpson(df['Rhodamine (ug/L)'].to_numpy(),dx=2)
    
    try:
        with open(files[index]) as fin:
            for i, row in enumerate(fin):
                if i == 0:
                    val = row[row.find("-")+2:].strip()[:-2]
                    ml_added = int(val)
    except Exception as e:
        ml_added = int(input("What was the input ml (int)"))
    
    ug_added = solution_ug_ml * ml_added 
    discharge = ug_added/area 
    
    ax[0].set_title(title + f' - Q:{discharge:.2f} l/s with {ml_added}ml')
    fig.savefig(f"{title}.png")
    return df, fig, val


for index, file in enumerate(files):
    df, fig, val = plot(index)
    
time.sleep(5)