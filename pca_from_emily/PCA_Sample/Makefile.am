bin_PROGRAMS = pcasift recalckeys gathergrads getpatches
pcasift_SOURCES = pcasift.cc image.cc keypoint.cc keypoint.h image.h
dist_data_DATA = gpcavects.txt docs/README LICENSE.DL Doxyfile
#mytests_SOURCES = mytests.cc image.cc keypoint.cc
recalckeys_SOURCES = recalckeys.cc image.cc keypoint.cc keypoint.h image.h
#testkeys_SOURCES = testkeys.cc image.cc keypoint.cc
gathergrads_SOURCES = gathergrads.cc keypoint.cc keypoint.h image.cc image.h
getpatches_SOURCES = getpatches.cc image.cc keypoint.cc image.h keypoint.h

LDFLAGS = -lnetpbm

MEGAOPT= -march=pentium4 -funroll-loops -fomit-frame-pointer -mfpmath=sse
CXXFLAGS = -Wall -O3 $(MEGAOPT)

# Memory check, but slow
#CXXFLAGS += -DMALLOC_CHECK_=2

# Profiling
#CXXFLAGS += -pg

# No asserts
CXXFLAGS += -DNDEBUG

# Add electric fence
#CXXFLAGS += -g
#LDFLAGS += -lefence
