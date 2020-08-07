Good day....

The main files in the folder are:
1.) BARTOLO_Final.m
	- this is the code for simulating shunt faults.
	- just follow carefully the instructions asked during the command line input because it's not designed to handle erroneous input. 
	  For example choosing phase B for change in symmetry must be equivalent to entering '2'
	- inputs should be real number only. No letter or other characters.
	- if the program crashes(probably due to user inputs) just re-run the program.
	- the default spreadsheet data has only 9 buses(hence 1 to 9 only for a bus number input). If user wants more buses he must add additional branches in the spreadsheet.
	- the default arrangement of the bus numbers and their possible base kV can be found in Test System.PNG
	--> Sample Input to command line:
		What is the base MVA?: 100
		What is the Bus Number for Base kV? 9
		What is the Base Line-Line kV for bus no.9 ? :138
		What is the faulted Bus Number? : 4
		Type of fault(1:3P, 2:LG, 3:LL, 4:DLG)? : 2
		What is the fault impedance(Ohms)?: 0
		What is the FAULTED Phase(1:PhaseA, 2:PhaseB, 3:PhaseC)?: 1

2.) FaultStudyData.xls
	- This is the spreadsheet input where the user can edit the data of the system
	- Do not change the name or file type of the spreadsheet.
	- Do not change the sheet name or headings of the worksheet
	- The MATLAB file can't handle erroneous data in the spreadhsheet. Hence user must put only valid data in the spreadhsheet
		- counting numbers only for bus numbers, FOC(0-11), transformer connection(1-Wye to ground, 2-Delta... no wye connection yet).
		- Positive real numbers only for kV/MVA rating.
		- Per unit impedances can be real or complex numbers.
	
	- for motors and generators only the steady value of the reactance will take effect in fault calculation(as of the moment)


Other files in the folder are:
1.) Test System - numbering of the buses for the default data in the spreadsheet.

2.) Bartolo_EE251_MP1_Documentation.pdf

3.) ETAP_files folder - contains the necessary files for running the 9bus Test System in the ETAP sofware.

4.) ETAP_[3phase/DLG/SLG/LL]_results.PNG - captured ETAP results for short cicuit analysis of the 9 bus Test System.

5.) ReadMe.txt - you are reading it. LOL.

Thank you.

