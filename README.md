#Silent Mutator
This perl script is designed to look for sites in a DNA sequence
which can be altered to either remove or introduce a restriction
site without altering the translated sequence.

The program can either be run as a single command or interactively.
To run the program from the command line enter:

	$silent_mutator.pl [OPTIONS]

To run the program interactively enter:

	$silent_mutator.pl -i

Then follow the on screen instructions.

Running the program without any options will prompt you to enter
a sequence. This will then be analysed for silent mutations using
the user's list of enzymes.

##Options

-b --biglist  : Use all of the enzymes in the big list. Better than pDRAW32!
-u --userlist : DEFAULT. Use only the enzymes in the user's list.
-o --oneenzyme: Use only the specified enzyme, eg -o EcoRI.

-v --view     : View the enzymes in the user list.
-a --add      : Add an enzyme to the user list FOR THIS RUN ONLY.
		Add multiple enzymes by separating them with a space,
		eg -a FokI AhaII BcgI

-j --justcut  : No mutations, just a list of enzymes cutting at least once.

-s --sequence : Enter the sequence IN FRAME, eg -s atgtgtagg. 

-h --help     : Print out this usage statement. Again.

##Output

All of the output can be redirected to a file of your choice with '>', eg:

	$silent_mutator.pl -s atgttgtgatagatga > newfile.txt


The output takes the form:
1. Sequence broken into codons
2. Translated codons
3. Enzymes which cut and the recognition site
4. One (or zero) examples of a removable site
5. Potential silent mutations
				
Regarding item 3, this list only includes the first site found for each enzyme.
Multiple sites per enzyme are not listed. The sequence is shown with a break 
highlighting the bases recognised by each enzyme. It is important to note that
some enzmyes do not cut at the site they recognise.

Items 4 and 5 use the same 'break' notation as item 3 to help identify the cut
site.


##License
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

A copy of the license is include in license.txt.

This license was chosen to encourage forks of the project to merge any updates or improvements.
