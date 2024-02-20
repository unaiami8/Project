#!/usr/bin/python3

import requests
import subprocess

API_URL = "https://dummyjson.com/todos?limit=1"

response = requests.get(API_URL)

if response.status_code == 200:
    data = response.json()
    for todo in data["todos"]:
        arg1 = todo['todo']
        subprocess.run(["module load slurm"], shell=True, capture_output=True, 
text=True)
        out = subprocess.run(["sbatch job.sbatch '" + arg1 + "'"], shell=True, 
capture_output=True, text=True)
        print(out)
else:
    print(f"API Error: {response.status_code}")
