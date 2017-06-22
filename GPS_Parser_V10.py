import os.path # Import os.path to search file paths
import tkinter as tk # to select files
from tkinter.filedialog import askopenfilename

'''
%=====================================================================================
% U S E R    P A R A M E T E R S
%=====================================================================================
'''
Vnum = 'V10'
PssblMsg = [3, 3500, 3501, 3502, 3512, 3623]
GPSModes = ['(1) Test ','(2) Initialization ','(3) *not used* ',
            '(4) Fine Alignment ','(5) Air Alignment ','(6) Transfer Alignment ',
            '(7) Air Navigation ','(8) *not used* ','(9) GPS Only ']

# Select directory and list file
tk.Tk().withdraw()
file = askopenfilename(title='GPS_Parser '+ Vnum + ' SELECT GPS or NAV FILE')
datadir, filename = os.path.split(file)
print(datadir)
print(filename)

# want figure? 
