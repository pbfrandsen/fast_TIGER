TIGER-rates
===========

The TIGER-rates software calculates the TIGER site rates introduced 
in (http://dx.doi.org/10.1093/sysbio/syr064) for a multiple sequence alignment.

The TIGER-rates program has been tested on Mac and Linux. It successfully
compiled using the clang++ as well as the g++ (version 4.2.1 and 4.7.3)
compilers.

There are three recommended ways to compile the TIGER-rates programm:
      (i)   Navigate to the "DAAD_project" directory and compile the program
            by typing the command
	    "clang++ -o TIGER-rates main.cpp"
	     or
	    "g++ -o TIGER-rates main.cpp"
            After this, you should have a new binary called "TIGER-rates"
	    in your directory.
      (ii)  Navigate to the main source folder and type make. After this,
            you should have a new binary called "TIGER-rates" in your directory.
      (iii) On Mac OSX, double click on the DAAD_project.xcodeproj file.
            This opens XCODE (if installed). Click on "Build and debug"
            to compile the TIGER-rates program. The executable can then be
	    found in the "build/Debug" folder.

The program works with phylip alignments. For an alignment called
"alignment.phy", type the command "./TIGER-rates dna alignment.phy"
to run the program.
The program will output a file called "alignment.phy_r8s.txt" that contains
the rates for each site separated by line breaks.

