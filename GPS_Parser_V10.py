import os.path # Import os.path to search file paths
import tkinter as tk # to select files
from tkinter.filedialog import askopenfilename
import tkinter.messagebox
from GPS_Parser_V10func import GPS_Parser_V10func
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
root = tk.Tk()
root.withdraw()
file = askopenfilename(title='GPS_Parser '+ Vnum + ' SELECT GPS or NAV FILE')
datadir, filename = os.path.split(file)
print(datadir)
print(filename)

# want figure? 
qstring = 'Output figure with GPS time from Mssg ID 3s ?'
choice = tkinter.messagebox.askquestion('FIGURE OUTPUT', qstring, icon='question')
if choice == 'yes': GPSFigs = True;
else: GPSFigs = False;

# extract GPS data?
qstring = 'Extract Msg3501 GPS data to CSV file?'
choice = tkinter.messagebox.askquestion('FIGURE OUTPUT', qstring, icon='question')
if choice == 'yes': PAR3501 = True;
else: PAR3501 = False;

# function generates the figures and parse file depending on user choice as 
# well as return data needed for Check_Flight script
[MSGCOUNT, REPCOUNT, GPSTimSec] = GPS_Parser_V10func(filename, datadir, GPSFigs, PAR3501)
