IDIR=includes/
CXX_FLAGS=-std=c++1y -fopenmp -O3

CPP_FILES := $(wildcard includes/*.cpp)
C_FILES := $(wildcard includes/*.c)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o) $(C_FILES:.c=.o)))
# We also depend on this generated file
OBJ_FILES += obj/cmdline.o

all: clf.run clus.run cmp.run

obj/%.o: includes/%.cpp
	@mkdir -p obj/
	$(CXX) $(CXX_FLAGS) -c -o $@ $<

obj/%.o: includes/%.c
	@mkdir -p obj/
	$(CXX) $(CXX_FLAGS) -c -o $@ $<

# Special generated file
includes/cmdline.c: includes/args_parser.gengetopt
	cd includes && gengetopt --input=args_parser.gengetopt --include-getopt 

clf.run: $(OBJ_FILES)
	$(CXX) $(CXX_FLAGS) -o $@ $^ clf.cpp -I$(IDIR)

clus.run: $(OBJ_FILES)
	$(CXX) $(CXX_FLAGS) -o $@ $^ clus_clf.cpp -I$(IDIR)

cmp.run: $(OBJ_FILES)
	$(CXX) $(CXX_FLAGS) -o $@ $^ compare.cpp -I$(IDIR)

clean:
	rm *.run
	rm -rf obj/
	rm includes/cmdline.*

.PHONY: clean all
