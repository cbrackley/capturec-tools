

CXX = c++

CFLAGS = -O3 -std=c++0x -Wall -g


src = $(wildcard *.cc)
obj = $(src:.cc=.o)
deps = $(obj:.o=.d)

directionality_SRC =	directionality.cc	\
			bedfiles.cc


all: directionality

directionality: $(directionality_SRC:.cc=.o)
	$(CXX) $(CFLAGS) -o $@ $^


%.o : %.cc
	$(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@ -MMD -MF $(@:.o=.d)

-include $(deps)

.PHONY: clean

clean:
	rm -f $(obj) $(data_to_beads_EXE) $(wildcard *.d)
