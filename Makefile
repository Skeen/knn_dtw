IDIR =includes/

all: clf clus cmp

clf:
	g++ -fopenmp -o clf.run -I$(IDIR) -O3 $(IDIR)ColMajorCell.cpp $(IDIR)FastDTW.cpp $(IDIR)WarpPath.cpp $(IDIR)SearchWindow.cpp $(IDIR)cmdline.c clf.cpp

clus:
	g++ -fopenmp -o clus_clf.run -I$(IDIR) -O3 $(IDIR)ColMajorCell.cpp $(IDIR)FastDTW.cpp $(IDIR)WarpPath.cpp $(IDIR)SearchWindow.cpp $(IDIR)cmdline.c clus_clf.cpp

cmp:
	g++ -fopenmp -o compare.run -I$(IDIR) -O3 $(IDIR)ColMajorCell.cpp $(IDIR)FastDTW.cpp $(IDIR)WarpPath.cpp $(IDIR)SearchWindow.cpp $(IDIR)cmdline.c compare.cpp

clean:
	rm *.run
