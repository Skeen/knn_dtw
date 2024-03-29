IDIR=includes/
CXX_FLAGS=-std=c++1y -fopenmp -O3

CPP_FILES := $(wildcard includes/*.cpp)
C_FILES := $(wildcard includes/*.c)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o) $(C_FILES:.c=.o)))
# We also depend on this generated file
OBJ_FILES += obj/cmdline.o
# We depend on object files
DEP_FILES := $(OBJ_FILES:.o=.d)
# We also depend on executable dependencies
DEP_FILES += obj/clf.d

all: clf.run

# Compile cpp files
obj/%.o: includes/%.cpp
	@mkdir -p obj/
	$(CXX) $(CXX_FLAGS) -MM -MT $@ -MF $(patsubst %.o,%.d,$@) $<
	$(CXX) $(CXX_FLAGS) -c -o $@ $<

# Compile c files
obj/%.o: includes/%.c
	@mkdir -p obj/
	$(CXX) $(CXX_FLAGS) -MM -MT $@ -MF $(patsubst %.o,%.d,$@) $<
	$(CXX) $(CXX_FLAGS) -c -o $@ $<

# Compile executable files
obj/%.o: %.cpp
	@mkdir -p obj/
	$(CXX) $(CXX_FLAGS) -MM -MT $@ -MF $(patsubst %.o,%.d,$@) $< -I$(IDIR)
	$(CXX) $(CXX_FLAGS) -c -o $@ $< -I$(IDIR)

# Special generated file
includes/cmdline.c: includes/args_parser.gengetopt
	cd includes && gengetopt --input=args_parser.gengetopt --include-getopt 

## Generate executables by linking
clf.run: $(OBJ_FILES) obj/clf.o
	$(CXX) $(CXX_FLAGS) -o $@ $^

run: clf.run
	./clf.run --query_filename=data/qry.job --reference_filename=data/ref.job

run_format: clf.run
	./clf.run --query_filename=data/qry.job --reference_filename=data/ref.job | python -m json.tool

clean:
	rm -f *.run
	rm -rf obj/
	rm -f includes/cmdline.*

-include $(DEP_FILES)

.PHONY: clean all run run_format
