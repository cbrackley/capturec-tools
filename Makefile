

CFLAGS = -O3 -std=c++0x -Wall -g

SRCDIR   = src
OBJDIR   = obj

SRC = $(wildcard $(SRCDIR)/*.cc)
OBJ = $(SRC:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)
DEPS = $(obj:.o=.d)

directionality_SRC =	directionality.cc	\
			bedfiles.cc

direct_derivative_SRC =	direct_derivative.cc	\
			bedfiles.cc

#prpn_in_window_SRC =	prpn_in_window.cc	\
#			bedfiles.cc


log_reads_v_separation_SRC = 	log_reads_v_separation.cc	\
				bedfiles.cc

local_v_long_SRC = 	local_v_long.cc	\
			bedfiles.cc

find_aretfacts_SRC =	find_aretfacts.cc	\
			bedfiles.cc

read_stats_SRC =	read_stats.cc	\
			bedfiles.cc

executables = directionality log_reads_v_separation local_v_long find_aretfacts direct_derivative read_stats #prpn_in_window

all: $(executables)


directionality: $(directionality_SRC:%.cc=$(OBJDIR)/%.o)
	$(CXX) $(CFLAGS) -o $@ $^

log_reads_v_separation: $(log_reads_v_separation_SRC:%.cc=$(OBJDIR)/%.o)
	$(CXX) $(CFLAGS) -o $@ $^

local_v_long: $(local_v_long_SRC:%.cc=$(OBJDIR)/%.o)
	$(CXX) $(CFLAGS) -o $@ $^

find_aretfacts: $(find_aretfacts_SRC:%.cc=$(OBJDIR)/%.o)
	$(CXX) $(CFLAGS) -o $@ $^

direct_derivative: $(direct_derivative_SRC:%.cc=$(OBJDIR)/%.o)
	$(CXX) $(CFLAGS) -o $@ $^

#prpn_in_window: $(prpn_in_window_SRC:%.cc=$(OBJDIR)/%.o)
#	$(CXX) $(CFLAGS) -o $@ $^

read_stats: $(read_stats_SRC:%.cc=$(OBJDIR)/%.o)
	$(CXX) $(CFLAGS) -o $@ $^


$(OBJDIR)/%.o : $(SRCDIR)/%.cc
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@ -MMD -MF $(@:.o=.d)

-include $(DEPS)

.PHONY: clean

clean:
	rm -f $(OBJ) $(executables) $(wildcard $(OBJDIR)/*.d) 
	rmdir $(OBJDIR)
