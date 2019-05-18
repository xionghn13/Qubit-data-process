import tkinter as tk
from tkinter import filedialog
import os

DataPath = 'E:/Projects\Fluxonium\data_process/Fluxonium042619/'

root = tk.Tk()
root.withdraw()


file_path = filedialog.askopenfilenames(initialdir=DataPath)
for file in file_path:
    print("\t'%s'," % os.path.basename(file))