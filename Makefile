IDIR =includes/

all: clf clus cmp

clf:
	g++ -std=c++1y -fopenmp -o clf.run -I$(IDIR) -O3 $(IDIR)ColMajorCell.cpp $(IDIR)FastDTW.cpp $(IDIR)WarpPath.cpp $(IDIR)SearchWindow.cpp $(IDIR)cmdline.c clf.cpp

clus:
	g++ -std=c++1y -fopenmp -o clus_clf.run -I$(IDIR) -O3 $(IDIR)ColMajorCell.cpp $(IDIR)FastDTW.cpp $(IDIR)WarpPath.cpp $(IDIR)SearchWindow.cpp $(IDIR)cmdline.c clus_clf.cpp

cmp:
	g++ -std=c++1y -fopenmp -o compare.run -I$(IDIR) -O3 $(IDIR)ColMajorCell.cpp $(IDIR)FastDTW.cpp $(IDIR)WarpPath.cpp $(IDIR)SearchWindow.cpp $(IDIR)cmdline.c compare.cpp

clean:
	rm *.run
