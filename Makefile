CXX := g++

NVCC := nvcc

CPPFLAGS := -std=c++17 -O2 -g

BUILDDIR := ./build
OBJDIR := $(BUILDDIR)/obj
TARGETDIR := $(BUILDDIR)/bin

SRC := main.cpp 		\
	   aberth.cpp		\
	   aberth_api.cpp	\
	   tests.cpp

OBJ := $(patsubst %.cu,  $(OBJDIR)/%.o, $(filter %.cu,  $(SRC))) \
       $(patsubst %.cpp, $(OBJDIR)/%.o, $(filter %.cpp, $(SRC)))

default : clean $(TARGETDIR)/aberth

clean :
	rm -f $(TARGETDIR)/* $(OBJDIR)/*

$(OBJDIR)/%.o : %.cpp
	@mkdir -p $(dir $@)
	@$(CXX) -c $(CPPFLAGS) $< -o $@

$(TARGETDIR)/aberth : $(OBJ)
	@mkdir -p $(dir $@)
	@$(CXX) $(OBJ) $(LDFLAGS) -o $@
