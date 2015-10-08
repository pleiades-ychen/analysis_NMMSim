Made an extra class of data structure, NMMEvent, to add some useful quantities
for furthure analysis. second_proc process the pre-processed root with NMMEvent structure. 

To use second_proc() with NMMEvent class, follow either way below. 


1) compilation

rootcint -f NMMEventDict.cc -c NMMEvent.hh  LinkDef.h

g++ -O2 `$ROOTSYS/bin/root-config --libs` -I$ROOTSYS/include *.cc -o second_proc

2) run as a script

root -l
[root] .L NMMEvent.cc+	// + for compiling to a library .so
[root] .L second_proc.cc
[root] second_proc(infilename, outfilename)
