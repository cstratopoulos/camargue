# Adapted from http://stackoverflow.com/a/27794283/6516346
#Compiler and Linker
CC          := /usr/local/Cellar/gcc49/4.9.3/bin/g++-4.9

#The Target Binary Program
TARGET      := PSEP

#The Directories, Source, Includes, Objects, Binary
SRCDIR      := source
INCDIR      := includes
BOOSTDIR    := /usr/local/boost_1_61_0
CPXDIR      := \
/Users/christos/Applications/IBM/ILOG/CPLEX_Studio1261/cplex/include/ilcplex/
PROGDIR     := /Users/christos/Dropbox/school/research/programs/
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
LIB         := /Users/christos/Applications/IBM/ILOG/CPLEX_Studio1261/cplex/lib/x86-64_osx/static_pic/libcplex.a \
/Users/christos/Dropbox/school/research/programs/concorde/concorde.a \
-lm -lpthread -fopenmp 
INC         := -I$(INCDIR) -I$(BOOSTDIR) -I$(CPXDIR) -I$(PROGDIR)
INCDEP      := -I$(INCDIR) -I$(BOOSTDIR) -I$(CPXDIR) -I$(PROGDIR)

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
	@$(RM) -f $(TARGET)

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
.PHONY: all remake clean
