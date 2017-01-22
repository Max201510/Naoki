bin_PROGRAMS = recalckeys gathergrads getpatches
#pcasift_SOURCES = pcasift.cc image.cc keypoint.cc keypoint.h image.h
dist_data_DATA = gpcavects.txt docs/README Doxyfile
recalckeys_SOURCES = recalckeys.cc image.cc keypoint.cc keypoint.h image.h
gathergrads_SOURCES = gathergrads.cc keypoint.cc keypoint.h image.cc image.h
getpatches_SOURCES = getpatches.cc image.cc keypoint.cc image.h keypoint.h

AM_LDFLAGS = -lnetpbm

#MEGAOPT= -march=pentium4 -funroll-loops -fomit-frame-pointer -mfpmath=sse
AM_CXXFLAGS = -Wall -O3 $(MEGAOPT)

# Memory check, but slow
#CXXFLAGS += -DMALLOC_CHECK_=2

# Profiling
#CXXFLAGS += -pg

# No asserts
AM_CXXFLAGS += -DNDEBUG

# Add electric fence
#CXXFLAGS += -g
#LDFLAGS += -lefence
