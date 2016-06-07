CC = g++
BFGS = bfgs.bin
BRENT = brent.bin
GENERIC = generic.bin

CFLAGS   = -std=c++11 -O3 -Wall

LINKER   = g++ -o

SRCDIR   = src
INCDIR   = include
OBJDIR   = obj
BINDIR   = bin

#SOURCES  := $(SRCDIR)/Optimization/Unrestricted/Multivariable/BFGS.cpp \
#			$(SRCDIR)/Optimization/Unrestricted/Singlevariable/Brent.cpp
INCLUDES := -I./$(INCDIR)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

all: $(BINDIR)/$(BFGS) $(BINDIR)/$(BRENT) $(BINDIR)/$(GENERIC)

$(BINDIR)/$(BFGS): $(OBJECTS)
	$(LINKER) $@ $(LFLAGS) $(INCLUDES) $(OBJECTS) examples/Example_BFGS.cpp
	@echo "Linking complete"

$(BINDIR)/$(BRENT): $(OBJECTS)
	$(LINKER) $@ $(LFLAGS) $(INCLUDES) $(OBJECTS) examples/Example_Brent.cpp
	@echo "Linking complete"

$(BINDIR)/$(GENERIC): $(OBJECTS)
	$(LINKER) $@ $(LFLAGS) $(INCLUDES) $(OBJECTS) examples/Example_Any.cpp
	@echo "Linking complete"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@
	@echo "Compiled "$<" successfully"