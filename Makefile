# Adapted from http://stackoverflow.com/a/27794283/6516346
#Compiler and Linker
CC          := g++

#The Target Binary Program
TARGET      := PSEP

#The Directories, Source, Includes, Objects, Binary and Resources
SRCDIR      := source
INCDIR      := includes
BOOSTDIR    := /usr/local/boost_1_61_0
CPXDIR      := \
/Users/christos/Applications/IBM/ILOG/CPLEX_Studio1261/cplex/include/ilcplex/
PROGDIR     := /Users/christos/Dropbox/school/research/programs/
BUILDDIR    := objects
TARGETDIR   := .
#RESDIR      := 
SRCEXT      := cpp
DEPEXT      := d
OBJEXT      := o

#Flags, Libraries and Includes
CFLAGS      := -Wall -O3 -pedantic -Wno-missing-braces -Wno-long-long\
-Wno-variadic-macros -Wno-gnu-zero-variadic-macro-arguments\
-Wno-gnu-statement-expression -std=c++11
LIB         := /Users/christos/Applications/IBM/ILOG/CPLEX_Studio1261/cplex/lib/x86-64_osx/static_pic/libcplex.a \
/Users/christos/Dropbox/school/research/programs/concorde/concorde.a \
-lm -lpthread
INC         := -I$(INCDIR) -I$(BOOSTDIR) -I$(CPXDIR) -I$(PROGDIR)
INCDEP      := -I$(INCDIR) -I$(BOOSTDIR) -I$(CPXDIR) -I$(PROGDIR)

#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------
SOURCES     := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))

#Defauilt Make
all: resources $(TARGET)

#Remake
remake: cleaner all

#Copy Resources from Resources Directory to Target Directory
#resources: directories
#    @cp $(RESDIR)/* $(TARGETDIR)/

#Make the Directories
directories:
	@mkdir -p $(TARGETDIR)
	@mkdir -p $(BUILDDIR)

#Clean only Objecst
clean:
	@$(RM) -rf $(BUILDDIR)

#Full Clean, Objects and Binaries
cleaner: clean
	@$(RM) -rf $(TARGETDIR)

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
.PHONY: all remake clean cleaner resources
