import tkinter as tk
from tkinter import filedialog
import os

DataPath = 'C:\\Users\\admin\Labber\Data\\2019\\09\Data_0910\\figures'

root = tk.Tk()
root.withdraw()


file_path = filedialog.askopenfilenames(initialdir=DataPath)
for file in file_path:
    print("\t'%s'," % os.path.basename(file))