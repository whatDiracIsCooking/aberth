# compiler
CXX := g++

# flags
CPPFLAGS := -std=c++17 -O2 -g

# directories
BUILDDIR := ./build
OBJDIR := $(BUILDDIR)/obj
TARGETDIR := $(BUILDDIR)/bin

# source files
SRC := main.cpp 		\
	   aberth.cpp		\
	   aberth_api.cpp	\
	   tests.cpp

OBJ := $(patsubst %.cpp, $(OBJDIR)/%.o, $(filter %.cpp, $(SRC)))

default : clean $(TARGETDIR)/aberth

clean :
	rm -f $(TARGETDIR)/* $(OBJDIR)/*

# create objects
$(OBJDIR)/%.o : %.cpp
	@mkdir -p $(dir $@)
	@$(CXX) -c $(CPPFLAGS) $< -o $@

# executable
$(TARGETDIR)/aberth : $(OBJ)
	@mkdir -p $(dir $@)
	@$(CXX) $(OBJ) -o $@
