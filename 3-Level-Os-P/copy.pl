#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 4.0 patchlevel 0
#    	last modified Thu Apr 15 14:44:22 CEST 2004
#    	System: Linux 2.6.18-1.2239.fc5
#    
#    	Copyright (C) 1986 - 1993, 1998, 2004
#    	Thomas Williams, Colin Kelley and many others
#    
#    	This is gnuplot version 4.0.  Please refer to the documentation
#    	for command syntax changes.  The old syntax will be accepted
#    	throughout the 4.0 series, but all save files use the new syntax.
#    
#    	Type `help` to access the on-line reference manual.
#    	The gnuplot FAQ is available from
#    		http://www.gnuplot.info/faq/
#    
#    	Send comments and requests for help to
#    		<gnuplot-info@lists.sourceforge.net>
#    	Send bugs, suggestions and mods to
#    		<gnuplot-bugs@lists.sourceforge.net>
#
set terminal postscript eps enhanced defaultplex \
   leveldefault color colortext \
   dashed dashlength 1.0 linewidth 1.0 butt \
   palfuncparam 2000,0.003 \
   "Times" 20
set output '3L-NoRWA-P.eps'
set multiplot
set size 1.0,0.3333
#
set origin 0.0,0.6666
set xtics border mirror norotate autofreq font "Times,20" 
set ytics border mirror norotate autofreq font "Times,20"
load 'population-NoRWA-3La-P.pl'
unset key
#
set origin 0.0,0.3333
set size 1,0.3333
set xtics border mirror norotate autofreq font "Times,20"
set ytics border mirror norotate autofreq font "Times,20"
load 'population-NoRWA-3Lb-P.pl'
#
set origin 0.0,0.0
set size 1,0.3333
set xtics border mirror norotate autofreq font "Times,20"
set ytics border mirror norotate autofreq font "Times,20"
load 'optfield-NoRWA-3L-P.pl'
#
unset multiplot
#EOF

