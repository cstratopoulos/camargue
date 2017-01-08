# Adapted from http://stackoverflow.com/a/27794283/6516346
#Compiler and Linker
CC          := g++

#The Target Binary Program
TARGET      := camargue

#The Directories, Source, Includes, Objects, Binary
SRCDIR      := source
INCDIR      := includes
CPXDIR      := /opt/ibm/ILOG/CPLEX_Studio127/cplex/include/ilcplex/
PROGDIR     := /home/christos/Documents/School/research/programs/
BUILDDIR    := objects
TARGETDIR   := .
SRCEXT      := cpp
DEPEXT      := d
OBJEXT      := o

#Flags, Libraries and Includes
CFLAGS      := -Wall -O3 -pedantic -fopenmp \
-Wno-missing-braces -Wno-sign-compare  -Wno-long-long\
-Wno-variadic-macros\
-std=c++11
LIB         := /opt/ibm/ILOG/CPLEX_Studio127/cplex/lib/x86-64_linux/static_pic/libcplex.a \
/home/christos/Documents/School/research/programs/concorde/concorde.a \
-lm -lpthread -fopenmp 
INC         := -I$(INCDIR)  -I$(CPXDIR) -I$(PROGDIR)
INCDEP      := -I$(INCDIR)  -I$(CPXDIR) -I$(PROGDIR)

#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------
SOURCES     := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))

#Defauilt Make
all: $(TARGET)

#Remake
remake: clean all

#Make the Directories
directories:
	@mkdir -p $(TARGETDIR)
	@mkdir -p $(BUILDDIR)

clean:
	@$(RM) -f $(BUILDDIR)/*.o
	@$(RM) -f $(BUILDDIR)/*.d
	@$(RM) -f $(BUILDDIR)/tests/*.o
	@$(RM) -f $(BUILDDIR)/tests/*.d
	@$(RM) -f $(TARGET)

# deftest:
# 	sed -e 's/#undef PSEP_DO_TESTS/#define PSEP_DO_TESTS/' \
# -i.back includes/tests.hpp && rm includes/tests.hpp.back

# undeftest:
# 	sed -e 's/#define PSEP_DO_TESTS/#undef PSEP_DO_TESTS/' \
# -i.back includes/tests.hpp && rm includes/tests.hpp.back

# test: deftest all undeftest

#Pull in dependency info for *existing* .o files
-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

#Link
$(TARGET): $(OBJECTS)
	$(CC) -o $(TARGETDIR)/$(TARGET) $^ $(LIB)

#Compile
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<
	@$(CC) $(CFLAGS) $(INCDEP) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp

#Non-File Targets
.PHONY: all remake clean #deftest undeftest test
