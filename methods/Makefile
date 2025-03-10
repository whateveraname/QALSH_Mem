# ------------------------------------------------------------------------------
#  Compile with C++ 11
# ------------------------------------------------------------------------------
SRCS=random.cc pri_queue.cc util.cc RSHJ.cpp
OBJS=${SRCS:.cc=.o}

CXX=g++ -std=c++11
CPPFLAGS=-w -O3 -DDO_PREFETCH

.PHONY: clean

all: ${OBJS}
	${CXX} ${CPPFLAGS} -o rshj ${OBJS}

clean:
	-rm ${OBJS} rshj
