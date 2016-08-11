PROJECT := hecate

CONFIG_FILE := Makefile.config

# Explicitly check for the config file, otherwise make -k will proceed anyway.
ifeq ($(wildcard $(CONFIG_FILE)),)
	$(error $(CONFIG_FILE) not found. See $(CONFIG_FILE).example.)
endif
include $(CONFIG_FILE)

BUILD_DIR_LINK := $(BUILD_DIR)
ifeq ($(RELEASE_BUILD_DIR),)
	RELEASE_BUILD_DIR := .$(BUILD_DIR)_release
endif
ifeq ($(DEBUG_BUILD_DIR),)
	DEBUG_BUILD_DIR := .$(BUILD_DIR)_debug
endif

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	BUILD_DIR := $(DEBUG_BUILD_DIR)
	OTHER_BUILD_DIR := $(RELEASE_BUILD_DIR)
else
	BUILD_DIR := $(RELEASE_BUILD_DIR)
	OTHER_BUILD_DIR := $(DEBUG_BUILD_DIR)
endif

# All of the directories containing code.
SRC_DIRS := $(shell find * -type d -exec bash -c "find {} -maxdepth 1 \
	\( -name '*.cpp' \) | grep -q ." \; -print)

# The target shared library name
LIB_BUILD_DIR := $(BUILD_DIR)/lib
STATIC_NAME := $(LIB_BUILD_DIR)/lib$(PROJECT).a
DYNAMIC_NAME := $(LIB_BUILD_DIR)/lib$(PROJECT).so

##############################
# Get all source files
##############################
# CXX_SRCS are the source files
ifeq ($(BUILD_TEXT), 1)
  CXX_SRCS := $(shell find src/$(PROJECT) -name "*.cpp")
else
  CXX_SRCS := $(filter-out src/$(PROJECT)/text.cpp, $(shell find src/$(PROJECT) -name "*.cpp"))
endif

# TOOL_SRCS are the source files for the tool binaries
ifeq ($(BUILD_TEXT), 1)
  TOOL_SRCS := $(shell find src/tools -name "*.cpp")
else
  TOOL_SRCS := $(filter-out src/tools/aesthetic_feature_extractor.cpp, $(shell find src/tools -name "*.cpp"))
endif

# BUILD_INCLUDE_DIR contains any generated header files we want to include.
BUILD_INCLUDE_DIR := $(BUILD_DIR)/src


##############################
# Derive generated files
##############################
# The generated files for protocol buffers
# The objects corresponding to the source files
# These objects will be linked into the final shared library, so we
# exclude the tool and test objects.
CXX_OBJS := $(addprefix $(BUILD_DIR)/, ${CXX_SRCS:.cpp=.o})
OBJS := $(CXX_OBJS)

# tool objects
TOOL_OBJS := $(addprefix $(BUILD_DIR)/, ${TOOL_SRCS:.cpp=.o})
TOOL_BUILD_DIR := $(BUILD_DIR)/tools

# Output files for automatic dependency generation
DEPS := ${CXX_OBJS:.o=.d}

# tool bins
TOOL_BINS := ${TOOL_OBJS:.o=}


##############################
# Derive compiler warning dump locations
##############################
WARNS_EXT := warnings.txt
CXX_WARNS := $(addprefix $(BUILD_DIR)/, ${CXX_SRCS:.cpp=.o.$(WARNS_EXT)})
TOOL_WARNS := $(addprefix $(BUILD_DIR)/, ${TOOL_SRCS:.cpp=.o.$(WARNS_EXT)})
ALL_CXX_WARNS := $(CXX_WARNS) $(TOOL_WARNS)
ALL_WARNS := $(ALL_CXX_WARNS)
NONEMPTY_WARN_REPORT := $(BUILD_DIR)/$(WARNS_EXT)


##############################
# Derive include and lib directories
##############################
INCLUDE_DIRS += $(BUILD_INCLUDE_DIR) ./include

LIBRARIES += opencv_core opencv_highgui opencv_imgproc
ifeq ($(OPENCV_VERSION), 3)
	LIBRARIES += opencv_imgcodecs opencv_videoio
endif

WARNINGS := -Wall -Wno-sign-compare


##############################
# Set build directories
##############################

DISTRIBUTE_DIR ?= distribute
DISTRIBUTE_SUBDIRS := $(DISTRIBUTE_DIR)/bin $(DISTRIBUTE_DIR)/lib
DIST_ALIASES := dist
ifneq ($(strip $(DISTRIBUTE_DIR)),distribute)
		DIST_ALIASES += distribute
endif

ALL_BUILD_DIRS := $(sort $(BUILD_DIR) $(addprefix $(BUILD_DIR)/, \
$(SRC_DIRS)) $(LIB_BUILD_DIR) $(DISTRIBUTE_SUBDIRS))


##############################
# Configure build
##############################

# Determine platform
UNAME := $(shell uname -s)
ifeq ($(UNAME), Linux)
	LINUX := 1
else ifeq ($(UNAME), Darwin)
	OSX := 1
endif

# Linux
ifeq ($(LINUX), 1)
	CXX ?= /usr/bin/g++
	GCCVERSION := $(shell $(CXX) -dumpversion | cut -f1,2 -d.)
	# older versions of gcc are too dumb to build boost with -Wuninitalized
	ifeq ($(shell echo | awk '{exit $(GCCVERSION) < 4.6;}'), 1)
		WARNINGS += -Wno-uninitialized
	endif
	# We will explicitly add stdc++ to the link target.
	LIBRARIES += stdc++
	COMMON_FLAGS += -fopenmp
endif

# OS X:
# clang++ instead of g++
ifeq ($(OSX), 1)
	CXX := /usr/bin/clang++
	# we need to explicitly ask for the rpath to be obeyed
	DYNAMIC_FLAGS := -install_name @rpath/lib$(PROJECT).so
	ORIGIN := @loader_path
else
	ORIGIN := \$$ORIGIN
endif

# Custom compiler
ifdef CUSTOM_CXX
	CXX := $(CUSTOM_CXX)
endif

# Debugging
ifeq ($(DEBUG), 1)
	COMMON_FLAGS += -std=c++0x -DDEBUG -g -O0
else
	COMMON_FLAGS += -std=c++0x -DNDEBUG -O2
endif

LIBRARY_DIRS += $(LIB_BUILD_DIR)

# Automatic dependency generation
CXXFLAGS += -MMD -MP

# Complete build flags.
COMMON_FLAGS += $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
CXXFLAGS += -pthread -fPIC $(COMMON_FLAGS) $(WARNINGS)
LINKFLAGS += -pthread -fPIC $(COMMON_FLAGS) $(WARNINGS)

LDFLAGS += $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(PKG_CONFIG) \
		$(foreach library,$(LIBRARIES),-l$(library))


##############################
# Define build targets
##############################
.PHONY: all lib tools clean $(DIST_ALIASES)

all: lib tools

lib: $(STATIC_NAME) $(DYNAMIC_NAME)

tools: $(TOOL_BINS)

$(ALL_WARNS): %.o.$(WARNS_EXT) : %.o

$(BUILD_DIR_LINK): $(BUILD_DIR)/.linked

# Create a target ".linked" in this BUILD_DIR to tell Make that the "build" link
# is currently correct, then delete the one in the OTHER_BUILD_DIR in case it
# exists and $(DEBUG) is toggled later.
$(BUILD_DIR)/.linked:
	@ mkdir -p $(BUILD_DIR)
	@ $(RM) $(OTHER_BUILD_DIR)/.linked
	@ $(RM) -r $(BUILD_DIR_LINK)
	@ ln -s $(BUILD_DIR) $(BUILD_DIR_LINK)
	@ touch $@

$(ALL_BUILD_DIRS): | $(BUILD_DIR_LINK)
	@ mkdir -p $@

$(DYNAMIC_NAME): $(OBJS) | $(LIB_BUILD_DIR)
	@ echo LD -o $@
	$(Q)$(CXX) -shared -o $@ $(OBJS) $(LINKFLAGS) $(LDFLAGS) $(DYNAMIC_FLAGS)

$(STATIC_NAME): $(OBJS) | $(LIB_BUILD_DIR)
	@ echo AR -o $@
	$(Q)ar rcs $@ $(OBJS)

$(BUILD_DIR)/%.o: %.cpp | $(ALL_BUILD_DIRS)
	@ echo CXX $<
	$(Q)$(CXX) $< $(CXXFLAGS) -c -o $@ 2> $@.$(WARNS_EXT) \
		|| (cat $@.$(WARNS_EXT); exit 1)
	@ cat $@.$(WARNS_EXT)

# Target for extension-less tool binaries
$(TOOL_BUILD_DIR)/%: $(TOOL_BUILD_DIR)/% | $(TOOL_BUILD_DIR)
	@ $(RM) $@
	@ ln -s $(abspath $<) $@

$(TOOL_BINS): % : %.o | $(DYNAMIC_NAME)
	@ echo CXX/LD -o $@
	$(Q)$(CXX) $< -o $@ $(LINKFLAGS) $(LDFLAGS) -l$(PROJECT)\
		-Wl,-rpath,$(ORIGIN)/../lib

clean:
	@- $(RM) -rf $(ALL_BUILD_DIRS)
	@- $(RM) -rf $(OTHER_BUILD_DIR)
	@- $(RM) -rf $(BUILD_DIR_LINK)
	@- $(RM) -rf $(DISTRIBUTE_DIR)

$(DIST_ALIASES): $(DISTRIBUTE_DIR)

$(DISTRIBUTE_DIR): all | $(DISTRIBUTE_SUBDIRS)
	# add include
	cp -r include $(DISTRIBUTE_DIR)/
	# add tool binaries
	cp $(TOOL_BINS) $(DISTRIBUTE_DIR)/bin
	# add libraries
	cp $(STATIC_NAME) $(DISTRIBUTE_DIR)/lib
	install -m 644 $(DYNAMIC_NAME) $(DISTRIBUTE_DIR)/lib

-include $(DEPS)
