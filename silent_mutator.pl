#!/usr/bin/perl

=begin license
Silent Mutator: Perl script to detect potential silent mutation sites in DNA sequences.
Copyright (C) 2012  Peter Hastie

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
=end license
=cut

use strict;
use warnings;
use Getopt::Long;

my $relative_location=&remove_path($0);
my $usage=
"

Silent Mutator - Help

This perl script is designed to look for sites in a DNA sequence
which can be altered to either remove or introduce a restriction
site without altering the translated sequence.

The program can either be run as a single command or interactively.
To run the program from the command line enter:

	$relative_location [OPTIONS]

To run the program interactively enter:

	$relative_location -i

Then follow the on screen instructions.

Running the program without any options will prompt you to enter
a sequence. This will then be analysed for silent mutations using
the user's list of enzymes.

##############
Options -

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

##############
Output -

All of the output can be redirected to a file of your choice with '>', eg:

	$relative_location -s atgttgtgatagatga > newfile.txt

The output takes the form 	1 - Sequence broken into codons
				2 - Translated codons
				3 - Enzymes which cut and the recognition site
				4 - One (or zero) examples of a removable site
				5 - Potential silent mutations
				
Regarding item 3, this list only includes the first site found for each enzyme.
Multiple sites per enzyme are not listed. The sequence is shown with a break 
highlighting the bases recognised by each enzyme. It is important to note that
some enzmyes do not cut at the site they recognise.

Items 4 and 5 use the same 'break' notation as item 3 to help identify the cut
site.

END OF HELP FILE
";


#Pete's restriction enzyme calculator
print "Silent Mutator Copyright (C) 2012 Peter Hastie\nThis program comes with ABSOLUTELY NO WARRANTY. This is free software, and you are welcome to redistribute it under certain conditions\n\n";

my %enzymes=(	"AaaI"=>"cggccg",
"AagI"=>"atcgat",
"AarI"=>"cacctgc",
"AasI"=>"gacnnnnnngtc",
"AatI"=>"aggcct",
"AatII"=>"gacgtc",
"AauI"=>"tgtaca",
"AbaI"=>"tgatca",
"AbeI"=>"cctcagc",
"AbrI"=>"ctcgag",
"AccI"=>"gtmkac",
"AccII"=>"cgcg",
"AccIII"=>"tccgga",
"Acc16I"=>"tgcgca",
"Acc36I"=>"acctgc",
"Acc65I"=>"ggtacc",
"Acc113I"=>"agtact",
"AccB1I"=>"ggyrcc",
"AccB2I"=>"rgcgcy",
"AccB7I"=>"ccannnnntgg",
"AccBSI"=>"ccgctc",
"AccEBI"=>"ggatcc",
"AceI"=>"gcwgc",
"AceII"=>"gctagc",
"AceIII"=>"cagctc",
"AciI"=>"ccgc",
"AclI"=>"aacgtt",
"AclNI"=>"actagt",
"AclWI"=>"ggatc",
"AcpI"=>"ttcgaa",
"AcpII"=>"ccannnnntgg",
"AcrII"=>"ggtnacc",
"AcsI"=>"raatty",
"AcuI"=>"ctgaag",
"AcvI"=>"cacgtg",
"AcyI"=>"grcgyc",
"AdeI"=>"cacnnngtg",
"AeuI"=>"ccwgg",
"AfaI"=>"gtac",
"Afa22MI"=>"cgatcg",
"Afa16RI"=>"cgatcg",
"AfeI"=>"agcgct",
"AflI"=>"ggwcc",
"AflII"=>"cttaag",
"AflIII"=>"acrygt",
"AgeI"=>"accggt",
"AglI"=>"ccwgg",
"AhaI"=>"ccsgg",
"AhaII"=>"grcgyc",
"AhaIII"=>"tttaaa",
"AhaB8I"=>"ggtacc",
"AhdI"=>"gacnnnnngtc",
"AhlI"=>"actagt",
"AhyI"=>"cccggg",
"AitI"=>"agcgct",
"AjnI"=>"ccwgg",
"AjoI"=>"ctgcag",
"AleI"=>"cacnnnngtg",
"AlfI"=>"gcannnnnntgc",
"AliI"=>"ggatcc",
"AliAJI"=>"ctgcag",
"AloI"=>"gaacnnnnnntcc",
"AluI"=>"agct",
"AlwI"=>"ggatc",
"Alw21I"=>"gwgcwc",
"Alw26I"=>"gtctc",
"Alw44I"=>"gtgcac",
"AlwNI"=>"cagnnnctg",
"AlwXI"=>"gcagc",
"Ama87I"=>"cycgrg",
"AocI"=>"cctnagg",
"AocII"=>"gdgchc",
"AorI"=>"ccwgg",
"Aor13HI"=>"tccgga",
"Aor51HI"=>"agcgct",
"AosI"=>"tgcgca",
"AosII"=>"grcgyc",
"ApaI"=>"gggccc",
"ApaBI"=>"gcannnnntgc",
"ApaCI"=>"ggatcc",
"ApaLI"=>"gtgcac",
"ApaORI"=>"ccwgg",
"ApeKI"=>"gcwgc",
"ApiI"=>"ctgcag",
"ApoI"=>"raatty",
"ApyI"=>"ccwgg",
"AquI"=>"cycgrg",
"AscI"=>"ggcgcgcc",
"AseI"=>"attaat",
"AseII"=>"ccsgg",
"AsiI"=>"ggatcc",
"AsiAI"=>"accggt",
"AsiSI"=>"gcgatcgc",
"AsnI"=>"attaat",
"AspI"=>"gacnnngtc",
"Asp700I"=>"gaannnnttc",
"Asp713I"=>"ctgcag",
"Asp718I"=>"ggtacc",
"Asp745I"=>"ggwcc",
"AspAI"=>"ggtnacc",
"AspA2I"=>"cctagg",
"AspEI"=>"gacnnnnngtc",
"AspHI"=>"gwgcwc",
"Asp10HI"=>"ttcgaa",
"Asp10HII"=>"ccannnnntgg",
"Asp26HI"=>"gaatgc",
"Asp27HI"=>"gaatgc",
"Asp35HI"=>"gaatgc",
"Asp36HI"=>"gaatgc",
"Asp40HI"=>"gaatgc",
"Asp50HI"=>"gaatgc",
"AspLEI"=>"gcgc",
"AspMI"=>"aggcct",
"AspMDI"=>"gatc",
"AspNI"=>"ggnncc",
"AspS9I"=>"ggncc",
"AssI"=>"agtact",
"AstWI"=>"grcgyc",
"AsuI"=>"ggncc",
"AsuII"=>"ttcgaa",
"AsuIII"=>"grcgyc",
"AsuC2I"=>"ccsgg",
"AsuHPI"=>"ggtga",
"AsuNHI"=>"gctagc",
"AtsI"=>"gacnnngtc",
"AvaI"=>"cycgrg",
"AvaII"=>"ggwcc",
"AvcI"=>"ggncc",
"AviII"=>"tgcgca",
"AvrII"=>"cctagg",
"AvrBII"=>"cctagg",
"AxyI"=>"cctnagg",
"Bac36I"=>"ggncc",
"BaeI"=>"acnnnngtayc",
"BalI"=>"tggcca",
"Bal228I"=>"ggncc",
"BamHI"=>"ggatcc",
"BamNxI"=>"ggwcc",
"BanI"=>"ggyrcc",
"BanII"=>"grgcyc",
"BanIII"=>"atcgat",
"BanAI"=>"ggcc",
"BasI"=>"ccannnnntgg",
"BauI"=>"cacgag",
"BavI"=>"cagctg",
"BavAI"=>"cagctg",
"BavAII"=>"ggncc",
"BavBI"=>"cagctg",
"BavBII"=>"ggncc",
"BavCI"=>"atcgat",
"BbeI"=>"ggcgcc",
"BbiII"=>"grcgyc",
"Bbi24I"=>"acgcgt",
"BbrI"=>"aagctt",
"Bbr7I"=>"gaagac",
"BbrPI"=>"cacgtg",
"BbsI"=>"gaagac",
"BbuI"=>"gcatgc",
"BbvI"=>"gcagc",
"BbvII"=>"gaagac",
"Bbv12I"=>"gwgcwc",
"Bbv16II"=>"gaagac",
"BbvAI"=>"gaannnnttc",
"BbvAII"=>"atcgat",
"BbvAIII"=>"tccgga",
"BbvBI"=>"ggyrcc",
"BbvCI"=>"cctcagc",
"Bca77I"=>"wccggw",
"BccI"=>"ccatc",
"Bce4I"=>"gcnnnnnnngc",
"Bce22I"=>"ggncc",
"Bce83I"=>"cttgag",
"Bce243I"=>"gatc",
"Bce751I"=>"ggatcc",
"BceAI"=>"acggc",
"BceBI"=>"cgcg",
"BceCI"=>"gcnnnnnnngc",
"BcefI"=>"acggc",
"BcgI"=>"cgannnnnntgc",
"Bci29I"=>"atcgat",
"BciBI"=>"atcgat",
"BciBII"=>"ccwgg",
"BciVI"=>"gtatcc",
"BclI"=>"tgatca",
"BcmI"=>"atcgat",
"BcnI"=>"ccsgg",
"BcoI"=>"cycgrg",
"Bco5I"=>"ctcttc",
"Bco27I"=>"ccgg",
"Bco116I"=>"ctcttc",
"Bco118I"=>"rccggy",
"BcoAI"=>"cacgtg",
"BcoKI"=>"ctcttc",
"BcuI"=>"actagt",
"BcuAI"=>"ggwcc",
"BdiI"=>"atcgat",
"BdiSI"=>"ctryag",
"BecAII"=>"ggcc",
"BepI"=>"cgcg",
"BetI"=>"wccggw",
"BfaI"=>"ctag",
"BfiI"=>"actggg",
"Bfi57I"=>"gatc",
"Bfi89I"=>"yggccr",
"BflI"=>"ccnnnnnnngg",
"BfmI"=>"ctryag",
"BfrI"=>"cttaag",
"BfrBI"=>"atgcat",
"BfuI"=>"gtatcc",
"BfuAI"=>"acctgc",
"BfuCI"=>"gatc",
"BglI"=>"gccnnnnnggc",
"BglII"=>"agatct",
"BimI"=>"ttcgaa",
"Bim19I"=>"ttcgaa",
"Bim19II"=>"ggcc",
"BinI"=>"ggatc",
"BlfI"=>"tccgga",
"Bli41I"=>"atcgat",
"Bli86I"=>"atcgat",
"Bli736I"=>"ggtctc",
"BliAI"=>"atcgat",
"BliHKI"=>"cctnagg",
"BliRI"=>"atcgat",
"BlnI"=>"cctagg",
"BloHI"=>"rgatcy",
"BloHII"=>"ctgcag",
"BlpI"=>"gctnagc",
"BluI"=>"ctcgag",
"Bme12I"=>"gatc",
"Bme18I"=>"ggwcc",
"Bme142I"=>"rgcgcy",
"Bme216I"=>"ggwcc",
"Bme361I"=>"ggcc",
"Bme585I"=>"cccgc",
"Bme1390I"=>"ccngg",
"Bme1580I"=>"gkgcmc",
"BmgBI"=>"cacgtc",
"BmrI"=>"actggg",
"BmtI"=>"gctagc",
"BmyI"=>"gdgchc",
"BnaI"=>"ggatcc",
"BoxI"=>"gacnnnngtc",
"BpcI"=>"ctryag",
"BpiI"=>"gaagac",
"BplI"=>"gagnnnnnctc",
"BpmI"=>"ctggag",
"BpoAI"=>"attaat",
"BptI"=>"ccwgg",
"BpuI"=>"grgcyc",
"Bpu10I"=>"cctnagc",
"Bpu14I"=>"ttcgaa",
"Bpu95I"=>"cgcg",
"Bpu1102I"=>"gctnagc",
"BpuAI"=>"gaagac",
"BpuAmI"=>"gagctc",
"BpuB5I"=>"cgtacg",
"BpuDI"=>"cctnagc",
"BpuEI"=>"cttgag",
"BpuSI"=>"gggac",
"BsaI"=>"ggtctc",
"Bsa29I"=>"atcgat",
"BsaAI"=>"yacgtr",
"BsaBI"=>"gatnnnnatc",
"BsaHI"=>"grcgyc",
"BsaJI"=>"ccnngg",
"BsaMI"=>"gaatgc",
"BsaOI"=>"cgrycg",
"BsaWI"=>"wccggw",
"BsaXI"=>"acnnnnnctcc",
"BscI"=>"atcgat",
"Bsc4I"=>"ccnnnnnnngg",
"Bsc91I"=>"gaagac",
"Bsc107I"=>"ccnnnnnnngg",
"BscAI"=>"gcatc",
"BscBI"=>"ggnncc",
"BscCI"=>"gaatgc",
"BscFI"=>"gatc",
"Bse1I"=>"actgg",
"Bse8I"=>"gatnnnnatc",
"Bse15I"=>"cycgrg",
"Bse16I"=>"ccwgg",
"Bse17I"=>"ccwgg",
"Bse21I"=>"cctnagg",
"Bse24I"=>"ccwgg",
"Bse64I"=>"ggtnacc",
"Bse118I"=>"rccggy",
"Bse634I"=>"rccggy",
"BseAI"=>"tccgga",
"BseBI"=>"ccwgg",
"BseCI"=>"atcgat",
"BseDI"=>"ccnngg",
"Bse3DI"=>"gcaatg",
"BseGI"=>"ggatg",
"BseJI"=>"gatnnnnatc",
"BseKI"=>"gcagc",
"BseLI"=>"ccnnnnnnngg",
"BseMI"=>"gcaatg",
"BseMII"=>"ctcag",
"BseNI"=>"actgg",
"BsePI"=>"gcgcgc",
"BseQI"=>"ggcc",
"BseRI"=>"gaggag",
"BseSI"=>"gkgcmc",
"BseT9I"=>"ggtnacc",
"BseT10I"=>"ggtnacc",
"BseXI"=>"gcagc",
"BseX3I"=>"cggccg",
"BseYI"=>"cccagc",
"BseZI"=>"ctcttc",
"BsgI"=>"gtgcag",
"BshI"=>"ggcc",
"Bsh45I"=>"gwgcwc",
"Bsh1236I"=>"cgcg",
"Bsh1285I"=>"cgrycg",
"Bsh1365I"=>"gatnnnnatc",
"BshFI"=>"ggcc",
"BshGI"=>"ccwgg",
"BshKI"=>"ggncc",
"BshNI"=>"ggyrcc",
"BshTI"=>"accggt",
"BsiI"=>"cacgag",
"BsiBI"=>"gatnnnnatc",
"BsiCI"=>"ttcgaa",
"BsiEI"=>"cgrycg",
"BsiHKAI"=>"gwgcwc",
"BsiHKCI"=>"cycgrg",
"BsiKI"=>"ggtnacc",
"BsiLI"=>"ccwgg",
"BsiMI"=>"tccgga",
"BsiQI"=>"tgatca",
"BsiSI"=>"ccgg",
"BsiWI"=>"cgtacg",
"BsiXI"=>"atcgat",
"BsiYI"=>"ccnnnnnnngg",
"BsiZI"=>"ggncc",
"BslI"=>"ccnnnnnnngg",
"BslFI"=>"gggac",
"BsmI"=>"gaatgc",
"BsmAI"=>"gtctc",
"BsmBI"=>"cgtctc",
"BsmFI"=>"gggac",
"BsmSI"=>"ccwwgg",
"Bso31I"=>"ggtctc",
"BsoBI"=>"cycgrg",
"BsoCI"=>"gdgchc",
"BsoFI"=>"gcngc",
"BsoMAI"=>"gtctc",
"Bsp6I"=>"gcngc",
"Bsp13I"=>"tccgga",
"Bsp19I"=>"ccatgg",
"Bsp24I"=>"gacnnnnnntgg",
"Bsp50I"=>"cgcg",
"Bsp63I"=>"ctgcag",
"Bsp67I"=>"gatc",
"Bsp68I"=>"tcgcga",
"Bsp98I"=>"ggatcc",
"Bsp105I"=>"gatc",
"Bsp106I"=>"atcgat",
"Bsp119I"=>"ttcgaa",
"Bsp120I"=>"gggccc",
"Bsp123I"=>"cgcg",
"Bsp143I"=>"gatc",
"Bsp143II"=>"rgcgcy",
"Bsp211I"=>"ggcc",
"Bsp423I"=>"gcagc",
"Bsp519I"=>"grgcyc",
"Bsp1286I"=>"gdgchc",
"Bsp1407I"=>"tgtaca",
"Bsp1720I"=>"gctnagc",
"Bsp1894I"=>"ggncc",
"Bsp2095I"=>"gatc",
"Bsp4009I"=>"ggatcc",
"BspAI"=>"gatc",
"BspA2I"=>"cctagg",
"Bsp153AI"=>"cagctg",
"BspAAI"=>"ctcgag",
"BspAAII"=>"tctaga",
"BspAAIII"=>"ggatcc",
"BspANI"=>"ggcc",
"BspBI"=>"ctgcag",
"BspBII"=>"ggncc",
"BspBRI"=>"ggcc",
"BspBS31I"=>"gaagac",
"BspCI"=>"cgatcg",
"BspCNI"=>"ctcag",
"BspDI"=>"atcgat",
"BspEI"=>"tccgga",
"BspFI"=>"gatc",
"BspF4I"=>"ggncc",
"BspHI"=>"tcatga",
"BspIS4I"=>"gaagac",
"BspJI"=>"gatc",
"BspJII"=>"atcgat",
"BspKI"=>"ggcc",
"BspKT5I"=>"ctgaag",
"BspKT6I"=>"gatc",
"BspKT8I"=>"aagctt",
"BspLI"=>"ggnncc",
"BspLAI"=>"gcgc",
"BspLAII"=>"ttcgaa",
"BspLAIII"=>"aagctt",
"BspLS2I"=>"gdgchc",
"BspLU4I"=>"cycgrg",
"BspLU11I"=>"acatgt",
"BspLU11III"=>"gggac",
"BspMI"=>"acctgc",
"BspMII"=>"tccgga",
"BspM39I"=>"cagctg",
"BspM90I"=>"gtatac",
"BspMAI"=>"ctgcag",
"BspMKI"=>"gtcgac",
"BspNI"=>"ccwgg",
"BspO4I"=>"cagctg",
"BspOVI"=>"gacnnnnngtc",
"BspOVII"=>"atcgat",
"BspPI"=>"ggatc",
"BspRI"=>"ggcc",
"BspR7I"=>"cctnagg",
"BspST5I"=>"gcatc",
"BspTI"=>"cttaag",
"BspT104I"=>"ttcgaa",
"BspT107I"=>"ggyrcc",
"BspTNI"=>"ggtctc",
"BspTS514I"=>"gaagac",
"BspWI"=>"gcnnnnnnngc",
"BspXI"=>"atcgat",
"BspXII"=>"tgatca",
"BspZEI"=>"atcgat",
"BsrI"=>"actgg",
"BsrAI"=>"ggwcc",
"BsrBI"=>"ccgctc",
"BsrBRI"=>"gatnnnnatc",
"BsrDI"=>"gcaatg",
"BsrFI"=>"rccggy",
"BsrGI"=>"tgtaca",
"BsrSI"=>"actgg",
"BssAI"=>"rccggy",
"BssECI"=>"ccnngg",
"BssHI"=>"ctcgag",
"BssHII"=>"gcgcgc",
"BssIMI"=>"gggtc",
"BssKI"=>"ccngg",
"BssNAI"=>"gtatac",
"BssSI"=>"cacgag",
"BssT1I"=>"ccwwgg",
"BstI"=>"ggatcc",
"Bst1I"=>"ccwgg",
"Bst2I"=>"ccwgg",
"Bst6I"=>"ctcttc",
"Bst11I"=>"actgg",
"Bst12I"=>"gcagc",
"Bst19I"=>"gcatc",
"Bst19II"=>"gatc",
"Bst28I"=>"atcgat",
"Bst38I"=>"ccwgg",
"Bst40I"=>"ccgg",
"Bst71I"=>"gcagc",
"Bst98I"=>"cttaag",
"Bst100I"=>"ccwgg",
"Bst1107I"=>"gtatac",
"BstACI"=>"grcgyc",
"BstAPI"=>"gcannnnntgc",
"BstAUI"=>"tgtaca",
"BstBI"=>"ttcgaa",
"Bst2BI"=>"cacgag",
"BstBAI"=>"yacgtr",
"BstBSI"=>"gtatac",
"BstB7SI"=>"rccggy",
"BstBS32I"=>"gaagac",
"BstBZ153I"=>"gcgcgc",
"Bst4CI"=>"acngt",
"BstC8I"=>"gcnngc",
"BstD102I"=>"ccgctc",
"BstDEI"=>"ctnag",
"BstDSI"=>"ccrygg",
"BstEII"=>"ggtnacc",
"BstENI"=>"cctnnnnnagg",
"BstENII"=>"gatc",
"BstEZ359I"=>"gttaac",
"BstFI"=>"aagctt",
"BstF5I"=>"ggatg",
"BstFNI"=>"cgcg",
"BstFZ438I"=>"cccgc",
"BstGZ53I"=>"cgtctc",
"BstH2I"=>"rgcgcy",
"BstH9I"=>"ggatc",
"BstHHI"=>"gcgc",
"BstHPI"=>"gttaac",
"BstHZ55I"=>"ccannnnnntgg",
"BstIZ316I"=>"cacnnngtg",
"BstJZ301I"=>"ctnag",
"BstKTI"=>"gatc",
"BstM6I"=>"ccwgg",
"BstMBI"=>"gatc",
"BstMCI"=>"cgrycg",
"BstMWI"=>"gcnnnnnnngc",
"BstMZ611I"=>"ccngg",
"BstNI"=>"ccwgg",
"Bst31NI"=>"ccgctc",
"BstNSI"=>"rcatgy",
"BstNZ169I"=>"atcgat",
"BstOI"=>"ccwgg",
"BstOZ616I"=>"gggac",
"BstPI"=>"ggtnacc",
"BstPAI"=>"gacnnnngtc",
"BstPZ418I"=>"ggatg",
"BstPZ740I"=>"cttaag",
"BstRZ246I"=>"atttaaat",
"BstSI"=>"cycgrg",
"BstSCI"=>"ccngg",
"BstSFI"=>"ctryag",
"BstSNI"=>"tacgta",
"BstSWI"=>"atttaaat",
"BstT7I"=>"tgatca",
"BstT9I"=>"ggtnacc",
"BstT10I"=>"ggtnacc",
"Bst31TI"=>"ggatc",
"BstTS5I"=>"gaagac",
"BstUI"=>"cgcg",
"Bst2UI"=>"ccwgg",
"BstVI"=>"ctcgag",
"BstV1I"=>"gcagc",
"BstV2I"=>"gaagac",
"BstXI"=>"ccannnnnntgg",
"BstX2I"=>"rgatcy",
"BstYI"=>"rgatcy",
"BstZI"=>"cggccg",
"BstZ17I"=>"gtatac",
"Bsu6I"=>"ctcttc",
"Bsu15I"=>"atcgat",
"Bsu23I"=>"tccgga",
"Bsu36I"=>"cctnagg",
"Bsu54I"=>"ggncc",
"Bsu1532I"=>"cgcg",
"Bsu1854I"=>"grgcyc",
"BsuBI"=>"ctgcag",
"BsuFI"=>"ccgg",
"BsuRI"=>"ggcc",
"BsuTUI"=>"atcgat",
"BteI"=>"ggcc",
"BtgI"=>"ccrygg",
"BtgZI"=>"gcgatg",
"BthAI"=>"ggwcc",
"BthCI"=>"gcngc",
"BthDI"=>"ccwgg",
"BthEI"=>"ccwgg",
"BtkI"=>"cgcg",
"BtkII"=>"gatc",
"BtrI"=>"cacgtc",
"BtsI"=>"gcagtg",
"BveI"=>"acctgc",
"BvuI"=>"grgcyc",
"BvuBI"=>"cgtacg",
"CacI"=>"gatc",
"Cac8I"=>"gcnngc",
"CaiI"=>"cagnnnctg",
"CauI"=>"ggwcc",
"CauII"=>"ccsgg",
"CauB3I"=>"tccgga",
"CbiI"=>"ttcgaa",
"CboI"=>"ccgg",
"CbrI"=>"ccwgg",
"CciNI"=>"gcggccgc",
"CcoI"=>"gccggc",
"CcrI"=>"ctcgag",
"CcuI"=>"ggncc",
"CcyI"=>"gatc",
"CdiI"=>"catcg",
"CelI"=>"ggatcc",
"CelII"=>"gctnagc",
"CeqI"=>"gatatc",
"CflI"=>"ctgcag",
"CfoI"=>"gcgc",
"CfrI"=>"yggccr",
"Cfr6I"=>"cagctg",
"Cfr9I"=>"cccggg",
"Cfr10I"=>"rccggy",
"Cfr13I"=>"ggncc",
"Cfr42I"=>"ccgcgg",
"CfrA4I"=>"ctgcag",
"CfrBI"=>"ccwwgg",
"CfrJ4I"=>"cccggg",
"CfuI"=>"gatc",
"CfuII"=>"ctgcag",
"ChaI"=>"gatc",
"CjeI"=>"ccannnnnngt",
"CjePI"=>"ccannnnnnntc",
"ClaI"=>"atcgat",
"CltI"=>"ggcc",
"CpfI"=>"gatc",
"CpoI"=>"cggwccg",
"CscI"=>"ccgcgg",
"CsiAI"=>"accggt",
"CsiBI"=>"gcggccgc",
"CspI"=>"cggwccg",
"Csp6I"=>"gtac",
"Csp45I"=>"ttcgaa",
"CspAI"=>"accggt",
"CspBI"=>"gcggccgc",
"CspCI"=>"caannnnngtgg",
"Csp68KI"=>"ggwcc",
"Csp68KII"=>"ttcgaa",
"Csp68KIII"=>"atgcat",
"Csp68KVI"=>"cgcg",
"CspKVI"=>"cgcg",
"CstI"=>"ctgcag",
"CstMI"=>"aaggag",
"CthII"=>"ccwgg",
"CviAI"=>"gatc",
"CviAII"=>"catg",
"CviBI"=>"gantc",
"CviJI"=>"rgcy",
"CviQI"=>"gtac",
"CviRI"=>"tgca",
"CviRII"=>"gtac",
"CviTI"=>"rgcy",
"CvnI"=>"cctnagg",
"DdeI"=>"ctnag",
"DmaI"=>"cagctg",
"DpaI"=>"agtact",
"DpnI"=>"gatc",
"DpnII"=>"gatc",
"DraI"=>"tttaaa",
"DraII"=>"rggnccy",
"DraIII"=>"cacnnngtg",
"DrdI"=>"gacnnnnnngtc",
"DriI"=>"gacnnnnngtc",
"DsaI"=>"ccrygg",
"DsaII"=>"ggcc",
"DsaIII"=>"rgatcy",
"DsaIV"=>"ggwcc",
"DsaV"=>"ccngg",
"DseDI"=>"gacnnnnnngtc",
"EacI"=>"ggatc",
"EaeI"=>"yggccr",
"Eae46I"=>"ccgcgg",
"EaeAI"=>"cccggg",
"EagI"=>"cggccg",
"EagBI"=>"cgatcg",
"EagMI"=>"ggwcc",
"Eam1104I"=>"ctcttc",
"Eam1105I"=>"gacnnnnngtc",
"EarI"=>"ctcttc",
"EcaI"=>"ggtnacc",
"EciI"=>"ggcgga",
"Eci125I"=>"ggtnacc",
"EclI"=>"cagctg",
"Ecl136II"=>"gagctc",
"EclHKI"=>"gacnnnnngtc",
"EclRI"=>"cccggg",
"EclXI"=>"cggccg",
"Ecl18kI"=>"ccngg",
"Ecl37kI"=>"ctgcag",
"Ecl2zI"=>"ctgcag",
"Eco24I"=>"grgcyc",
"Eco31I"=>"ggtctc",
"Eco32I"=>"gatatc",
"Eco47I"=>"ggwcc",
"Eco47III"=>"agcgct",
"Eco52I"=>"cggccg",
"Eco56I"=>"gccggc",
"Eco57I"=>"ctgaag",
"Eco64I"=>"ggyrcc",
"Eco72I"=>"cacgtg",
"Eco78I"=>"ggcgcc",
"Eco81I"=>"cctnagg",
"Eco88I"=>"cycgrg",
"Eco91I"=>"ggtnacc",
"Eco105I"=>"tacgta",
"Eco130I"=>"ccwwgg",
"Eco147I"=>"aggcct",
"Eco255I"=>"agtact",
"Eco1831I"=>"ccsgg",
"EcoA4I"=>"ggtctc",
"EcoHI"=>"ccsgg",
"EcoHK31I"=>"yggccr",
"EcoICRI"=>"gagctc",
"Eco75KI"=>"grgcyc",
"Eco57MI"=>"ctgrag",
"EcoNI"=>"cctnnnnnagg",
"EcoO44I"=>"ggtctc",
"EcoO65I"=>"ggtnacc",
"EcoO109I"=>"rggnccy",
"EcoO128I"=>"ggtnacc",
"EcoRI"=>"gaattc",
"EcoRII"=>"ccwgg",
"EcoRV"=>"gatatc",
"EcoT14I"=>"ccwwgg",
"EcoT22I"=>"atgcat",
"EcoT38I"=>"grgcyc",
"EcoVIII"=>"aagctt",
"Eco13kI"=>"ccngg",
"Eco21kI"=>"ccngg",
"Eco27kI"=>"cycgrg",
"Eco29kI"=>"ccgcgg",
"Eco53kI"=>"gagctc",
"Eco137kI"=>"ccngg",
"EgeI"=>"ggcgcc",
"EheI"=>"ggcgcc",
"ErhI"=>"ccwwgg",
"ErhB9I"=>"cgatcg",
"ErhB9II"=>"ccwwgg",
"ErpI"=>"ggwcc",
"EsaBC3I"=>"tcga",
"EsaBC4I"=>"ggcc",
"EspI"=>"gctnagc",
"Esp3I"=>"cgtctc",
"Esp4I"=>"cttaag",
"Esp1396I"=>"ccannnnntgg",
"FalI"=>"aagnnnnnctt",
"FalII"=>"cgcg",
"FaqI"=>"gggac",
"FatI"=>"catg",
"FauI"=>"cccgc",
"FauBII"=>"cgcg",
"FauNDI"=>"catatg",
"FbaI"=>"tgatca",
"FblI"=>"gtmkac",
"FbrI"=>"gcngc",
"FdiI"=>"ggwcc",
"FdiII"=>"tgcgca",
"FgoI"=>"ctag",
"FmuI"=>"ggncc",
"FnuAI"=>"gantc",
"FnuCI"=>"gatc",
"FnuDI"=>"ggcc",
"FnuDII"=>"cgcg",
"FnuDIII"=>"gcgc",
"FnuEI"=>"gatc",
"Fnu4HI"=>"gcngc",
"FokI"=>"ggatg",
"FriOI"=>"grgcyc",
"FseI"=>"ggccggcc",
"FsiI"=>"raatty",
"FspI"=>"tgcgca",
"FspII"=>"ttcgaa",
"Fsp1604I"=>"ccwgg",
"FspAI"=>"rtgcgcay",
"FspBI"=>"ctag",
"Fsp4HI"=>"gcngc",
"FspMSI"=>"ggwcc",
"FssI"=>"ggwcc",
"FunI"=>"agcgct",
"FunII"=>"gaattc",
"GalI"=>"ccgcgg",
"GceI"=>"ccgcgg",
"GceGLI"=>"ccgcgg",
"GdiI"=>"aggcct",
"GdiII"=>"cggccr",
"GstI"=>"ggatcc",
"GsuI"=>"ctggag",
"HacI"=>"gatc",
"HaeI"=>"wggccw",
"HaeII"=>"rgcgcy",
"HaeIII"=>"ggcc",
"HaeIV"=>"gaynnnnnrtc",
"HalI"=>"gaattc",
"HalII"=>"ctgcag",
"HapII"=>"ccgg",
"HgaI"=>"gacgc",
"HgiI"=>"grcgyc",
"HgiAI"=>"gwgcwc",
"HgiBI"=>"ggwcc",
"HgiCI"=>"ggyrcc",
"HgiCII"=>"ggwcc",
"HgiCIII"=>"gtcgac",
"HgiDI"=>"grcgyc",
"HgiDII"=>"gtcgac",
"HgiEI"=>"ggwcc",
"HgiGI"=>"grcgyc",
"HgiHI"=>"ggyrcc",
"HgiHII"=>"grcgyc",
"HgiHIII"=>"ggwcc",
"HgiJI"=>"ggwcc",
"HgiJII"=>"grgcyc",
"HgiS22I"=>"ccsgg",
"HhaI"=>"gcgc",
"HhaII"=>"gantc",
"Hin1I"=>"grcgyc",
"Hin1II"=>"catg",
"Hin2I"=>"ccgg",
"Hin4I"=>"gaynnnnnvtc",
"Hin6I"=>"gcgc",
"HinJCI"=>"gtyrac",
"HinP1I"=>"gcgc",
"HincII"=>"gtyrac",
"HindII"=>"gtyrac",
"HindIII"=>"aagctt",
"HinfI"=>"gantc",
"HjaI"=>"gatatc",
"HpaI"=>"gttaac",
"HpaII"=>"ccgg",
"HphI"=>"ggtga",
"Hpy8I"=>"gtnnac",
"Hpy51I"=>"gtsac",
"Hpy99I"=>"cgwcg",
"Hpy178III"=>"tcnnga",
"Hpy188I"=>"tcnga",
"Hpy188III"=>"tcnnga",
"HpyAV"=>"ccttc",
"HpyBI"=>"gtac",
"HpyBII"=>"gtnnac",
"HpyCI"=>"gatatc",
"HpyC1I"=>"ccatc",
"HpyCH4I"=>"catg",
"HpyCH4III"=>"acngt",
"HpyCH4IV"=>"acgt",
"HpyCH4V"=>"tgca",
"HpyF10VI"=>"gcnnnnnnngc",
"HpyF44III"=>"tgca",
"HsoI"=>"gcgc",
"Hsp92I"=>"grcgyc",
"Hsp92II"=>"catg",
"HspAI"=>"gcgc",
"HsuI"=>"aagctt",
"ItaI"=>"gcngc",
"I-CeuI"=>"taactataacggtcctaaggtagcga",
"KasI"=>"ggcgcc",
"Kaz48kI"=>"rggnccy",
"KoxII"=>"grgcyc",
"KpnI"=>"ggtacc",
"Kpn2I"=>"tccgga",
"Kpn378I"=>"ccgcgg",
"Kpn2kI"=>"ccngg",
"Kpn49kI"=>"gaattc",
"Kpn49kII"=>"ccsgg",
"KspI"=>"ccgcgg",
"Ksp22I"=>"tgatca",
"Ksp632I"=>"ctcttc",
"KspAI"=>"gttaac",
"Kzo9I"=>"gatc",
"Kzo49I"=>"ggwcc",
"LcaI"=>"atcgat",
"LlaAI"=>"gatc",
"LlaBI"=>"ctryag",
"LlaCI"=>"aagctt",
"LlaG2I"=>"gctagc",
"Lmu60I"=>"cctnagg",
"LplI"=>"atcgat",
"LpnI"=>"rgcgcy",
"LspI"=>"ttcgaa",
"LweI"=>"gcatc",
"MabI"=>"accwggt",
"MaeI"=>"ctag",
"MaeII"=>"acgt",
"MaeIII"=>"gtnac",
"MaeK81I"=>"cgtacg",
"MaeK81II"=>"ggncc",
"MamI"=>"gatnnnnatc",
"MavI"=>"ctcgag",
"MbiI"=>"ccgctc",
"MboI"=>"gatc",
"MboII"=>"gaaga",
"MchI"=>"ggcgcc",
"MchAI"=>"gcggccgc",
"MchAII"=>"ggcc",
"McrI"=>"cgrycg",
"MfeI"=>"caattg",
"MflI"=>"rgatcy",
"MfoAI"=>"ggcc",
"Mgl14481I"=>"ccsgg",
"MgoI"=>"gatc",
"MhaAI"=>"ctgcag",
"MhlI"=>"gdgchc",
"MkrAI"=>"gatc",
"MlaI"=>"ttcgaa",
"MlaAI"=>"ctcgag",
"MlsI"=>"tggcca",
"MltI"=>"agct",
"MluI"=>"acgcgt",
"Mlu23I"=>"ggatcc",
"Mlu31I"=>"tggcca",
"MluB2I"=>"tcgcga",
"MluNI"=>"tggcca",
"MlyI"=>"gagtc",
"Mly113I"=>"ggcgcc",
"MmeI"=>"tccrac",
"MnlI"=>"cctc",
"MnoI"=>"ccgg",
"Mph1103I"=>"atgcat",
"MroI"=>"tccgga",
"MroNI"=>"gccggc",
"MroXI"=>"gaannnnttc",
"MscI"=>"tggcca",
"MseI"=>"ttaa",
"MslI"=>"caynnnnrtg",
"MspI"=>"ccgg",
"Msp17I"=>"grcgyc",
"Msp20I"=>"tggcca",
"Msp67I"=>"ccngg",
"MspA1I"=>"cmgckg",
"MspB4I"=>"ggyrcc",
"MspCI"=>"cttaag",
"MspR9I"=>"ccngg",
"MspSWI"=>"atttaaat",
"MspV281I"=>"gwgcwc",
"MspYI"=>"yacgtr",
"MssI"=>"gtttaaac",
"MstI"=>"tgcgca",
"MstII"=>"cctnagg",
"MthZI"=>"ctag",
"MunI"=>"caattg",
"MvaI"=>"ccwgg",
"Mva1269I"=>"gaatgc",
"MvnI"=>"cgcg",
"MvrI"=>"cgatcg",
"MwoI"=>"gcnnnnnnngc",
"MxaI"=>"gagctc",
"NaeI"=>"gccggc",
"NarI"=>"ggcgcc",
"NblI"=>"cgatcg",
"NciI"=>"ccsgg",
"NcoI"=>"ccatgg",
"NcrI"=>"agatct",
"NcuI"=>"gaaga",
"NdaI"=>"ggcgcc",
"NdeI"=>"catatg",
"NdeII"=>"gatc",
"NgoAIII"=>"ccgcgg",
"NgoAIV"=>"gccggc",
"NgoMIV"=>"gccggc",
"NgoPII"=>"ggcc",
"NgoPIII"=>"ccgcgg",
"NheI"=>"gctagc",
"NlaII"=>"gatc",
"NlaIII"=>"catg",
"NlaIV"=>"ggnncc",
"Nli3877I"=>"cycgrg",
"NmeCI"=>"gatc",
"NmeRI"=>"cagctg",
"NmuCI"=>"gtsac",
"NopI"=>"gtcgac",
"NotI"=>"gcggccgc",
"NphI"=>"gatc",
"NruI"=>"tcgcga",
"NruGI"=>"gacnnnnngtc",
"NsbI"=>"tgcgca",
"NsiI"=>"atgcat",
"NsiCI"=>"gatatc",
"NspI"=>"rcatgy",
"NspII"=>"gdgchc",
"NspIII"=>"cycgrg",
"NspIV"=>"ggncc",
"NspV"=>"ttcgaa",
"Nsp7121I"=>"ggncc",
"Nsp29132II"=>"ggatcc",
"NspBII"=>"cmgckg",
"NspHI"=>"rcatgy",
"NspLKI"=>"ggcc",
"NspMACI"=>"agatct",
"NspSAI"=>"cycgrg",
"NspSAII"=>"ggtnacc",
"NspSAIV"=>"ggatcc",
"NunII"=>"ggcgcc",
"OfoI"=>"cycgrg",
"OkrAI"=>"ggatcc",
"OliI"=>"cacnnnngtg",
"OxaNI"=>"cctnagg",
"PacI"=>"ttaattaa",
"Pac25I"=>"cccggg",
"PaeI"=>"gcatgc",
"PaeAI"=>"ccgcgg",
"PaeBI"=>"cccggg",
"PaeHI"=>"grgcyc",
"PaePI"=>"ctgcag",
"PaeQI"=>"ccgcgg",
"PaeR7I"=>"ctcgag",
"Pae2kI"=>"agatct",
"Pae5kI"=>"ccgcgg",
"Pae14kI"=>"ccgcgg",
"Pae17kI"=>"cagctg",
"Pae18kI"=>"agatct",
"PagI"=>"tcatga",
"PalI"=>"ggcc",
"PamI"=>"tgcgca",
"PamII"=>"grcgyc",
"PanI"=>"ctcgag",
"PasI"=>"cccwggg",
"PauI"=>"gcgcgc",
"PauAI"=>"rcatgy",
"PauAII"=>"tttaaa",
"PceI"=>"aggcct",
"PciI"=>"acatgt",
"PctI"=>"gaatgc",
"Pde12I"=>"ggncc",
"Pde133I"=>"ggcc",
"Pde137I"=>"ccgg",
"PdiI"=>"gccggc",
"PdmI"=>"gaannnnttc",
"PfaAI"=>"ggyrcc",
"PfaAII"=>"catatg",
"PfaAIII"=>"gcatgc",
"PfeI"=>"gawtc",
"Pfl8I"=>"ggatcc",
"Pfl21I"=>"ctgcag",
"Pfl23II"=>"cgtacg",
"Pfl27I"=>"rggwccy",
"PflBI"=>"ccannnnntgg",
"PflFI"=>"gacnnngtc",
"PflKI"=>"ggcc",
"PflMI"=>"ccannnnntgg",
"PfoI"=>"tccngga",
"PgaI"=>"atcgat",
"PhaI"=>"gcatc",
"PhoI"=>"ggcc",
"PinAI"=>"accggt",
"PinBI"=>"atgcat",
"PinBII"=>"tccgga",
"PlaI"=>"ggcc",
"PlaII"=>"ttcgaa",
"PlaAI"=>"cycgrg",
"PlaAII"=>"gtac",
"PleI"=>"gagtc",
"Ple19I"=>"cgatcg",
"PmaCI"=>"cacgtg",
"PmeI"=>"gtttaaac",
"Pme55I"=>"aggcct",
"PmlI"=>"cacgtg",
"PpaAI"=>"ttcgaa",
"PpaAII"=>"tcga",
"PpeI"=>"gggccc",
"PpiI"=>"gaacnnnnnctc",
"PpsI"=>"gagtc",
"Ppu10I"=>"atgcat",
"Ppu111I"=>"gaattc",
"PpuAI"=>"cgtacg",
"PpuMI"=>"rggwccy",
"PpuXI"=>"rggwccy",
"PshAI"=>"gacnnnngtc",
"PshBI"=>"attaat",
"PsiI"=>"ttataa",
"Psp03I"=>"ggwcc",
"Psp5II"=>"rggwccy",
"Psp6I"=>"ccwgg",
"Psp23I"=>"ctgcag",
"Psp1406I"=>"aacgtt",
"PspAI"=>"cccggg",
"PspALI"=>"cccggg",
"Psp124BI"=>"gagctc",
"PspCI"=>"cacgtg",
"PspEI"=>"ggtnacc",
"PspGI"=>"ccwgg",
"PspLI"=>"cgtacg",
"PspN4I"=>"ggnncc",
"PspOMI"=>"gggccc",
"PspPI"=>"ggncc",
"PspPPI"=>"rggwccy",
"PspXI"=>"vctcgagb",
"PsrI"=>"gaacnnnnnntac",
"PssI"=>"rggnccy",
"PstI"=>"ctgcag",
"PstNHI"=>"gctagc",
"PsuI"=>"rgatcy",
"Psu161I"=>"cgatcg",
"PsuAI"=>"yacgtr",
"PsyI"=>"gacnnngtc",
"PtaI"=>"tccgga",
"Pun14627I"=>"tgcgca",
"Pun14627II"=>"cagctg",
"PunAI"=>"cycgrg",
"PunAII"=>"rcatgy",
"PvuI"=>"cgatcg",
"PvuII"=>"cagctg",
"Pvu84II"=>"cagctg",
"PI-SceI"=>"atctatgtcgggtgcggagaaagaggtaatgaaatggca",
"RalF40I"=>"gatc",
"RcaI"=>"tcatga",
"RflFI"=>"gtcgac",
"RflFII"=>"agtact",
"RleAI"=>"cccaca",
"RmaI"=>"ctag",
"Rme21I"=>"atcgat",
"RsaI"=>"gtac",
"RshI"=>"cgatcg",
"RspLKI"=>"gcatgc",
"RspLKII"=>"ggatcc",
"RspXI"=>"tcatga",
"RsrI"=>"gaattc",
"RsrII"=>"cggwccg",
"Rsr2I"=>"cggwccg",
"RtrI"=>"gtcgac",
"Rtr63I"=>"gtcgac",
"SacI"=>"gagctc",
"SacII"=>"ccgcgg",
"SacNI"=>"grgcyc",
"SalI"=>"gtcgac",
"SalPI"=>"ctgcag",
"SanDI"=>"gggwccc",
"SapI"=>"gctcttc",
"SarI"=>"aggcct",
"SatI"=>"gcngc",
"SauI"=>"cctnagg",
"Sau96I"=>"ggncc",
"Sau3239I"=>"ctcgag",
"Sau3AI"=>"gatc",
"SauBMKI"=>"gccggc",
"SauHPI"=>"gccggc",
"SauLPI"=>"gccggc",
"SauLPII"=>"ctcgag",
"SauMI"=>"gatc",
"SauNI"=>"gccggc",
"SauSI"=>"gccggc",
"SbfI"=>"cctgcagg",
"Sbi68I"=>"ctcgag",
"Sbo13I"=>"tcgcga",
"SbvI"=>"ggcc",
"ScaI"=>"agtact",
"SceIII"=>"gccggc",
"SchI"=>"gagtc",
"SchZI"=>"ccgcgg",
"SciI"=>"ctcgag",
"SciNI"=>"gcgc",
"ScrFI"=>"ccngg",
"SdaI"=>"cctgcagg",
"SdiI"=>"ggccnnnnnggcc",
"SduI"=>"gdgchc",
"SecI"=>"ccnngg",
"SelI"=>"cgcg",
"SenPT16I"=>"cggccg",
"SenPT14bI"=>"ccgcgg",
"SepI"=>"atgcat",
"SexAI"=>"accwggt",
"SexBI"=>"ccgcgg",
"SexCI"=>"ccgcgg",
"SfaI"=>"ggcc",
"SfaNI"=>"gcatc",
"SfcI"=>"ctryag",
"SfeI"=>"ctryag",
"SfiI"=>"ggccnnnnnggcc",
"SflI"=>"ctgcag",
"SfoI"=>"ggcgcc",
"Sfr274I"=>"ctcgag",
"Sfr303I"=>"ccgcgg",
"SfuI"=>"ttcgaa",
"SgfI"=>"gcgatcgc",
"SgrAI"=>"crccggyg",
"SgrBI"=>"ccgcgg",
"SimI"=>"gggtc",
"SinI"=>"ggwcc",
"SlaI"=>"ctcgag",
"SleI"=>"ccwgg",
"Slu1777I"=>"gccggc",
"SmaI"=>"cccggg",
"SmiI"=>"atttaaat",
"SmiMI"=>"caynnnnrtg",
"SmlI"=>"ctyrag",
"SmuI"=>"cccgc",
"SmuEI"=>"ggwcc",
"SnaBI"=>"tacgta",
"SniI"=>"ccwgg",
"SnoI"=>"gtgcac",
"SolI"=>"ggatcc",
"Sol10179I"=>"ctcgag",
"SpaHI"=>"gcatgc",
"SpeI"=>"actagt",
"SphI"=>"gcatgc",
"SplI"=>"cgtacg",
"SpmI"=>"atcgat",
"SpoI"=>"tcgcga",
"SpuI"=>"ccgcgg",
"SrfI"=>"gcccgggc",
"SrlI"=>"gccggc",
"Srl5DI"=>"ctgcag",
"Srl32DII"=>"gaattc",
"Srl55DI"=>"gaattc",
"Srl56DI"=>"ctryag",
"SruI"=>"tttaaa",
"Sru4DI"=>"attaat",
"Sru30DI"=>"aggcct",
"SsbI"=>"aagctt",
"SscL1I"=>"gantc",
"Sse9I"=>"aatt",
"Sse232I"=>"cgccggcg",
"Sse1825I"=>"gggwccc",
"Sse8387I"=>"cctgcagg",
"Sse8647I"=>"aggwcct",
"SseAI"=>"ggcgcc",
"SseBI"=>"aggcct",
"SshAI"=>"cctnagg",
"SsiI"=>"ccgc",
"SsiAI"=>"gatc",
"SsiBI"=>"gatc",
"SslI"=>"ccwgg",
"SsoI"=>"gaattc",
"SsoII"=>"ccngg",
"SspI"=>"aatatt",
"Ssp1I"=>"ttcgaa",
"Ssp4800I"=>"tgtaca",
"Ssp5230I"=>"gacgtc",
"Ssp27144I"=>"atcgat",
"SspAI"=>"ccwgg",
"SspBI"=>"tgtaca",
"SspCI"=>"gccggc",
"SspD5I"=>"ggtga",
"SspD5II"=>"atgcat",
"SspRFI"=>"ttcgaa",
"SsrI"=>"gttaac",
"SstI"=>"gagctc",
"SstII"=>"ccgcgg",
"Sst12I"=>"ctgcag",
"SteI"=>"aggcct",
"SthI"=>"ggtacc",
"Sth117I"=>"ccwgg",
"Sth132I"=>"cccg",
"Sth134I"=>"ccgg",
"Sth368I"=>"gatc",
"StrI"=>"ctcgag",
"StsI"=>"ggatg",
"StuI"=>"aggcct",
"StyI"=>"ccwwgg",
"StyD4I"=>"ccngg",
"SuaI"=>"ggcc",
"SunI"=>"cgtacg",
"SurI"=>"ggatcc",
"SviI"=>"ttcgaa",
"SwaI"=>"atttaaat",
"TaaI"=>"acngt",
"TaiI"=>"acgt",
"TaqI"=>"tcga",
"TaqII"=>"gaccgacaccca",
"Taq52I"=>"gcwgc",
"TaqXI"=>"ccwgg",
"TasI"=>"aatt",
"TatI"=>"wgtacw",
"TauI"=>"gcsgc",
"TelI"=>"gacnnngtc",
"TfiI"=>"gawtc",
"ThaI"=>"cgcg",
"TliI"=>"ctcgag",
"Tru1I"=>"ttaa",
"Tru9I"=>"ttaa",
"Tru201I"=>"rgatcy",
"TscI"=>"acgt",
"TseI"=>"gcwgc",
"Tsp1I"=>"actgg",
"Tsp32I"=>"tcga",
"Tsp32II"=>"tcga",
"Tsp45I"=>"gtsac",
"Tsp49I"=>"acgt",
"Tsp509I"=>"aatt",
"TspBI"=>"ccrygg",
"Tsp4CI"=>"acngt",
"TspDTI"=>"atgaa",
"TspEI"=>"aatt",
"Tsp8EI"=>"gccnnnnnggc",
"TspGWI"=>"acgga",
"TspRI"=>"castgnn",
"Tth111I"=>"gacnnngtc",
"Tth111II"=>"caarca",
"TthHB8I"=>"tcga",
"Uba4009I"=>"ggatcc",
"Uba153AI"=>"cagctg",
"UbaM39I"=>"cagctg",
"UnbI"=>"ggncc",
"UthSI"=>"cccggg",
"Uur960I"=>"gcngc",
"Van91I"=>"ccannnnntgg",
"Vha464I"=>"cttaag",
"VneI"=>"gtgcac",
"VpaK32I"=>"gctcttc",
"VpaK11AI"=>"ggwcc",
"VpaK11BI"=>"ggwcc",
"VspI"=>"attaat",
"XagI"=>"cctnnnnnagg",
"XapI"=>"raatty",
"XbaI"=>"tctaga",
"XcaI"=>"gtatac",
"XceI"=>"rcatgy",
"XciI"=>"gtcgac",
"XcmI"=>"ccannnnnnnnntgg",
"XcyI"=>"cccggg",
"XhoI"=>"ctcgag",
"XhoII"=>"rgatcy",
"XmaI"=>"cccggg",
"XmaIII"=>"cggccg",
"XmaCI"=>"cccggg",
"XmaJI"=>"cctagg",
"XmiI"=>"gtmkac",
"XmnI"=>"gaannnnttc",
"XorII"=>"cgatcg",
"XpaI"=>"ctcgag",
"XspI"=>"ctag",
"YenI"=>"ctgcag",
"ZanI"=>"ccwgg",
"ZhoI"=>"atcgat",
"ZraI"=>"gacgtc",
"ZrmI"=>"agtact",
"Zsp2I"=>"atgcat",

		);


my @user_list=('AflII','AgeI','ApaLI','AseI','BamHI','BglII','BsaBI','BspHI','BstBI','ClaI','EagI','EcoRI','EcoRV','HindIII','HpaI','KpnI','MfeI','NcoI','NdeI','NheI','NotI','PacI','PstI','SacI','SalI','ScaI','SacII','SmaI','SpeI','SphI','XbaI','XhoI');
my $use_list='u'; #want to set u as the default
my $input_sequence='';
my $use_interactive=0;
my @enzymes_to_add=();
my $justcut=0;
my $one_enzyme='';

#allow use of Getopt::Long
my $getopt_variable = GetOptions("userlist" => sub{$use_list='u'}, 
				 "biglist"=> sub{$use_list='b'},
				 "sequence=s" => \$input_sequence,
				 "interactive" => \$use_interactive,
				 "view" => sub {print "Enzymes in user list:\n\n",join ("\n",@user_list),"\n\n"},
				 "add:s{,}" => \@enzymes_to_add, #this will group zero or more arguments following -a
				 "justcut" => \$justcut,
				 "oneenzyme=s" => \$one_enzyme,
				 "help" => sub{print $usage,"\n\n"},
				 );

#add enzymes if necessary
if (scalar(@enzymes_to_add)>0){
	print "Adding to enzymes in user list:\n\n";
	foreach my $add_to_user_list(@enzymes_to_add){
		print "$add_to_user_list $enzymes{$add_to_user_list}\n";
		push (@user_list,$add_to_user_list);
	}
	print "\n";
}

#fall through to interactive mode if use_interactive chosen
my $rerun_main_options=1;

if ($use_interactive==1){
	
	while($rerun_main_options>0){
		&main_options(\@user_list);
	}
	if ($input_sequence eq ''){
		print "\nPlease enter a sequence: \n";
		$input_sequence=<STDIN>;
	}
}

#give a chance to enter sequence without going into the full menu
if ($input_sequence eq ''){
	
	print "Please enter a sequence: \n";
	$input_sequence=<STDIN>;
}


chomp($input_sequence);
chomp($use_list);
die "Enter at least three codons\n" if (length($input_sequence)<9);
$input_sequence=lc($input_sequence);
die "Non atgc characters in sequence\n" if ($input_sequence=~m/[^atgc]/);

chomp($one_enzyme);
		
		
my %codon_table = (tca => 'ser',tcg => 'ser',tcc => 'ser',tct => 'ser',
			   ttt => 'phe',ttc => 'phe',tta => 'leu',ttg => 'leu',
			   tat => 'tyr',tac => 'tyr',taa => 'stop',tag => 'stop',
			   tgt => 'cys',tgc => 'cys',tga => 'stop',tgg => 'trp',
			   cta => 'leu',ctg => 'leu',ctc => 'leu',ctt => 'leu',
			   cca => 'pro',ccg => 'pro',ccc => 'pro',cct => 'pro',
			   cat => 'his',cac => 'his',caa => 'gln',cag => 'gln',
			   cga => 'arg',cgg => 'arg',cgc => 'arg',cgt => 'arg',
			   att => 'ile',atc => 'ile',ata => 'ile',atg => 'met',
			   aca => 'thr',acg => 'thr',acc => 'thr',act => 'thr',
			   aat => 'asn',aac => 'asn',aaa => 'lys',aag => 'lys',
			   agt => 'ser',agc => 'ser',aga => 'arg',agg => 'arg',
			   gta => 'val',gtg => 'val',gtc => 'val',gtt => 'val',
			   gca => 'ala',gcg => 'ala',gcc => 'ala',gct => 'ala',
			   gat => 'asp',gac => 'asp',gaa => 'glu',gag => 'glu',
			   gga => 'gly',ggg => 'gly',ggc => 'gly',ggt => 'gly');
		
my @ala=('gct','gcc','gca','gcg'); 
my @leu=('tta','ttg','ctt','ctc','cta','ctg'); 
my @arg=('cgt','cgc','cga','cgg','aga','agg'); 
my @lys=('aaa','aag'); 
my @asn=('aat','aac'); 
my @met=('atg'); 
my @asp=('gat','gac'); 
my @phe=('ttt','ttc'); 
my @cys=('tgt','tgc'); 
my @pro=('cct','ccc','cca','ccg'); 
my @gln=('caa','cag'); 
my @ser=('tct','tcc','tca','tcg','agt','agc'); 
my @glu=('gaa','gag'); 
my @thr=('act','acc','aca','acg'); 
my @gly=('ggt','ggc','gga','ggg'); 
my @trp=('tgg');
my @his=('cat','cac'); 
my @tyr=('tat','tac'); 
my @ile=('att','atc','ata'); 
my @val=('gtt','gtc','gta','gtg'); 
my @stop=('tag','tga','taa');

#have the codons as a hash where each element of the hash is an array

my %codons=(	'ala' => \@ala,
		'leu' => \@leu,
		'arg' => \@arg,
		'lys' => \@lys,
		'asn' => \@asn,
		'met' => \@met,
		'asp' => \@asp,
		'phe' => \@phe,
		'cys' => \@cys,
		'pro' => \@pro,
		'gln' => \@gln,
		'ser' => \@ser,
		'glu' => \@glu,
		'thr' => \@thr,
		'gly' => \@gly,
		'trp' => \@trp,
		'his' => \@his,
		'tyr' => \@tyr,
		'ile' => \@ile,
		'val' => \@val,
		'stop'=> \@stop,
		);
		
#as a test, loop through the hash and print out everything
#my @hash_keys=keys(%codons);
#foreach my $key(@hash_keys){
#	print $key,": ";
#	
#	#now print all of the contents of that hash value
#	print join (', ',@{$codons{$key}}),"\n";	
#}
		

#copy the input_sequence and change it into a 3 codon version
my $input2codon=$input_sequence;
$input2codon=~s/(...)/$1 /g;
print "\nInput:      ",$input2codon,"\n\n";

#translate the input sequence
my @translated_sequence=&translate_in_frame($input_sequence);
print "Translated: ",join(' ',@translated_sequence),"\n\n";

#convert degenerate enzymes into proper regular expressions
&enzyme_to_regex;

#find existing cut sites
print "Currently contains the following cut sites:\n";
my @existing_cutsites=&find_existing_cuts($input_sequence);
foreach my $cut(@existing_cutsites){
	if ($justcut==0){
		print join(" ",@{$cut}),"\n";
	}
	else{
		print $cut->[0],"\n";
	}
}
print"\n";

if ($justcut==0){
#try to destroy those cut sites
my @lost_site=&destroy_cutsites(@existing_cutsites);
#if there are no cut sites say that
if (scalar(@lost_site)==0){
	print "No existing cut sites found\n\n";
}
else{
	#make it obvious what the codons are to highlight change
	my $new_site_codons=$lost_site[1].$lost_site[2].$lost_site[3];
	$new_site_codons=~s/(...)/$1 /g;
	print "Lose $lost_site[0] cut site if it is changed to $new_site_codons\n\n";
}

#make mutations to find new cut sites
my %cutsites_from_alternatives=&make_new_cutsites($input_sequence);
print "Found these silent mutations:\n\n";
my $number_mutations_found=0;
#if any of these alternative cutsites are novel, display that
my %cutsites_from_original=();
foreach my $cut (@existing_cutsites){
	$cutsites_from_original{$cut->[0]}='';
}
my @cuts_from_alternatives=keys(%cutsites_from_alternatives);
foreach my $cut (@cuts_from_alternatives){
	unless (exists $cutsites_from_original{$cut}){
		print "$cut $cutsites_from_alternatives{$cut}\n\n";
		$number_mutations_found++;
	}
}
if ($number_mutations_found==0){
	print"No silent mutations found\n";
}
}


#####################################################################
sub translate_in_frame{
	my ($input_sequence_passed)=@_;
	$input_sequence_passed=~s/(...)/$codon_table{$1} /g;
	my @translated=split / /,$input_sequence_passed;
	return @translated;
}
	
		
#####################################################################				

sub find_existing_cuts{
	my ($input_sequence_passed)=@_;
	
	my @existing_cuts;
	my @enzyme_names;
	if ($use_list eq 'b'){
		@enzyme_names=keys(%enzymes);
	}
	if ($use_list eq 'u'){
		@enzyme_names=@user_list;
	}
	if ($one_enzyme ne ''){
		@enzyme_names=$one_enzyme;
	}
	foreach my $enzyme(@enzyme_names){
			
		if($input_sequence_passed=~m/$enzymes{$enzyme}/){
			my @cutter=($enzyme,$`,$&,$');
			#push a reference to the @cutter into @existing_cuts
			push(@existing_cuts,\@cutter);
		}
	}
	return @existing_cuts;
}

######################################################################

sub destroy_cutsites{

	my @existing_cuts=@_;
	
	foreach my $cut (@existing_cuts){
		my $cut_name=$cut->[0];
		my $cut_before=$cut->[1];
		my $cut_site=$cut->[2];
	
		#work out what the smallest workable fragment is
		my $firstcodon=int(length($cut_before)/3);
		my $firstbase=$firstcodon*3;
			
		my $lastcodon_length=(length($cut_before)+length($cut_site));
		my $lastcodon;
		if (($lastcodon_length)%3==0){
			$lastcodon=$lastcodon_length/3;
		}
		else{
			$lastcodon=(int($lastcodon_length/3))+1;
		}
		my $lastbase=$lastcodon*3;
		
		my $sequence_around_site=substr($input_sequence,$firstbase,($lastbase-$firstbase));
		#print "sequence_around_site $sequence_around_site\n";
		
		#now try the alternatives for each codon until a match is found which destroys the cutsite
		my $codons_around_site=$sequence_around_site;
		$codons_around_site=~s/(...)/$1 /g;
		my @codons_around_site_array=split(/ /,$codons_around_site);
		
		#Track which codon is being altered so sequence can be rebuilt
		my @codons_looked_at=();
		my $codon_being_looked_at='';
		my @codons_to_look_at=@codons_around_site_array;
		
		my $size_of_codons_to_look_at = scalar(@codons_to_look_at);
				
		for (my $i=0 ; $i<$size_of_codons_to_look_at;$i++){
			$codon_being_looked_at=shift(@codons_to_look_at);
			my $residue=$codon_table{$codon_being_looked_at};
			my @alternatives=@{$codons{$residue}};
			#so now go through each of the alternatives, put the flanking sequences in and see if cutsite is gone
			foreach my $alternative(@alternatives){
				my $reconstituted_sequence=join('',@codons_looked_at).$alternative.join('',@codons_to_look_at);
				unless ($reconstituted_sequence=~m/$cut_site/){
					my @lost_cut_site=($cut_name, join('',@codons_looked_at),$alternative,join('',@codons_to_look_at));
					return @lost_cut_site;
				}
			}
			push(@codons_looked_at,$codon_being_looked_at);
		}
		
		
		
			 
	}
}

#############################################################################

sub make_new_cutsites{
	
	my ($input_sequence)=@_;
	
	my %new_cutsites=(); #for output
	
	$input_sequence=~s/(...)/$1 /g;
	my @codons_to_mutate=split(/ /,$input_sequence);
	#if the last entry is not a full codon, pop it off
	unless (length($codons_to_mutate[-1])==3){
		my $lost_codon_fragment=pop (@codons_to_mutate);
	}
	
	#Track which codon is being altered so sequence can be rebuilt
	my @codons_looked_at=();
	my @codons_being_looked_at=();
	
	my $size_of_codons_to_mutate = scalar(@codons_to_mutate);
	
	#will be looking at three codons in a row
	#shift the first two residues outside the loop, then just shift one at a time as normal
	my $codon_being_shifted=shift(@codons_to_mutate);
	push (@codons_being_looked_at,$codon_being_shifted);
	$codon_being_shifted=shift(@codons_to_mutate);
	push (@codons_being_looked_at,$codon_being_shifted);
	
	for (my $i=0; $i<($size_of_codons_to_mutate-2);$i++){
		$codon_being_shifted=shift(@codons_to_mutate);
		push (@codons_being_looked_at,$codon_being_shifted);
		#print join(' ',@codons_being_looked_at),"\n";
		
		my $residue=$codon_table{$codons_being_looked_at[0]};
		my $next_residue=$codon_table{$codons_being_looked_at[1]};
		my $nextnext_residue=$codon_table{$codons_being_looked_at[2]};
		
		my @alternatives=@{$codons{$residue}};
		my @next_alternatives=@{$codons{$next_residue}};
		my @nextnext_alternatives=@{$codons{$nextnext_residue}};
		
		#go through alternatives, reconstitute, if a cutsite appears that's not already in @existing_cutsites, make a note of it
		foreach my $alternative1(@alternatives){
			foreach my $alternative2(@next_alternatives){
				foreach my $alternative3(@nextnext_alternatives){
					
					my $reconstituted_sequence=join('',@codons_looked_at).$alternative1.$alternative2.$alternative3.join('',@codons_to_mutate);	
					#now use the reconstituted_sequence as input for the find_existing_cuts subroutine
					my @cutsites_in_alternative=&find_existing_cuts($reconstituted_sequence);
						
					#if anything in @cutsites_in_alternative is different from @existing_cutsites (and hasn't already been added), add it to @new_cutsites
			
					foreach my $cut(@cutsites_in_alternative){
						#check that the new_cutsites hash doesn't already contain the enzyme and add it
						unless (exists $new_cutsites{$cut->[0]}){
							$new_cutsites{$cut->[0]}=$cut->[1].' '.$cut->[2].' '.$cut->[3];	
						}
				
					}
				}	
			}
		}
		my $completed_codon=shift(@codons_being_looked_at);
		push(@codons_looked_at,$completed_codon);
	}
	return %new_cutsites;
}

############################################################
sub main_options{
	my ($ref_to_user_list)=@_;
			
	print "\nTo use enzymes already in the user list enter 'u'\nTo use enzymes from the big list type 'b'\nTo view enzymes in the user list type 'v'\nTo add an enzyme to the user list type 'a'\n";
	print "To use only one specified enzyme enter 'o'\n";
	print "To only produce a list of enzymes that cut at least once enter 'j'\n";
	print "For help enter 'h'\n";
	if ($input_sequence eq ''){print "When ready to enter sequence type 'done'\n\n";}
	else {print "When ready to run type 'done'\n\n";}
	my $option=<STDIN>;
	chomp($option);
	
	if ($option eq 'u'){$use_list='u'}
	if ($option eq 'b'){$use_list='b'}
	if ($option eq 'v'){print join ("\n",@$ref_to_user_list),"\n";}
	if ($option eq 'a'){&add_enzyme_to_user_list($ref_to_user_list)}
	if ($option eq 'j'){$justcut=1}
	if ($option eq 'o'){&set_one_enzyme}
	if ($option eq 'h'){print $usage,"\n\n"}
	if ($option eq 'done'){$rerun_main_options=0}
	
}

#############################################################
sub add_enzyme_to_user_list{
	my ($ref_to_user_list)=@_;
	
	print "Add enzyme name without spaces, eg EcoRI or type Done\n";
	my $new_enzyme=<STDIN>;
	chomp($new_enzyme);
	
	if (exists $enzymes{$new_enzyme}){
		push @$ref_to_user_list, $new_enzyme;
		print "$new_enzyme: $enzymes{$new_enzyme} temporarily added to the user list\n";
	}
	else{
		print "Sorry, $new_enzyme wasn't recognised. Check the enzyme name or add it to %enzymes in the code\n";
	}
}

################################
sub enzyme_to_regex{
	my @enzyme_names;
	if ($use_list eq 'b'){
		@enzyme_names=keys(%enzymes);
	}
	if ($use_list eq 'u'){
		@enzyme_names=@user_list;
	}
	foreach my $enzyme(@enzyme_names){
		$enzymes{$enzyme}=~s/r/\(a|g\)/g;
		$enzymes{$enzyme}=~s/w/\(a|t\)/g;
		$enzymes{$enzyme}=~s/y/\(c|t\)/g;
		$enzymes{$enzyme}=~s/n/\(c|t|a|g\)/g;
		$enzymes{$enzyme}=~s/m/\(a|c\)/g;
		$enzymes{$enzyme}=~s/k/\(g|t\)/g;
		$enzymes{$enzyme}=~s/s/\(c|g\)/g;
		$enzymes{$enzyme}=~s/v/\(a|g|c\)/g;
		$enzymes{$enzyme}=~s/h/\(a|c|t\)/g;
		$enzymes{$enzyme}=~s/d/\(a|g|t\)/g;
		$enzymes{$enzyme}=~s/b/\(g|c|t\)/g;
	}
}

##################################
sub set_one_enzyme{
	print "Enter name of enzyme to use:\n";
	$one_enzyme=<STDIN>;
}

##################################
sub remove_path{
	my ($location) = @_;
	my @location_fields=split ('/',$location);
	my $relative_location=$location_fields[-1];
	return $relative_location;
}
