SRCS = src/type/ghquads.cpp \
	src/utils/Input.cpp \
	src/type/dataset.cpp \
	src/estimation/estep.cpp \
	src/estimation/mstep.cpp \
	src/utils/asa111.cpp

LIBRARY = irtppcpp
SRC_DIR = src
CPPFLAGS = -std=c++11 -march=native -g3 -Wall -fPIC
INCLUDES = -I./$(SRC_DIR) -I./include/SPGO/include/
OBJS = $(SRCS:.cpp=.o)

all : lib$(LIBRARY).a

$(OBJS): %.o : %.h

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	g++ $(CPPFLAGS) -c $(INCLUDES) -o $@ $<

lib$(LIBRARY).a : $(OBJS)
	ar rsv $@ $^

clean:
	$(RM) $(OBJS)
	$(RM) lib$(LIBRARY).a
