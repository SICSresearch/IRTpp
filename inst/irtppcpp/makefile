SRC = src/type/ghquads.cpp \
	src/utils/Input.cpp \
	src/type/dataset.cpp \
	src/estimation/estep.cpp \
	src/estimation/mstep.cpp \
        src/utils/asa111.cpp
OBJ = $(SRC:.cpp=.o)
INCL = -I./src/ -I./include/SPGO/include/
CFLAGS = -std=c++11 -Wunused-function -O3
TINCL = -I./tests/ -I./include/
TSRC = tests/matrix.cpp \
		tests/test.cpp


all: SICS

$(OBJ): %.o : %.h

.cpp.o:
src/%.o: src/%.cpp
	g++ $(CFLAGS) $(INCL) -Wall -c -I./src/ -o $@ $<

SICS: $(OBJ)
	g++ -Wall $(CFLAGS) $(INCL) $^ src/main.cpp -o $@
#SICS:
#	g++ -Wall $(CFLAGS) $(INCL) $^ src/main.cpp -o $@

#debug:
#	g++ $(CFLAGS) -ggdb -I./src/ -o SICS_dbg $(SRC) src/main.cpp

test:
	g++ $(TINCL) -o $@1  $(TSRC)
	./$@1
	rm $@1

clean:
	rm SICS $(OBJ)
	rm test1 SICS $(OBJ)

citest:
	g++-4.8 $(TINCL) $(CFLAGS) -o $@1  $(TSRC)
	./$@1
