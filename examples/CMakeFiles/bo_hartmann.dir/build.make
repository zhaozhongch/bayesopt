# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/zhaozhong/bayesopt

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/zhaozhong/bayesopt

# Include any dependencies generated for this target.
include examples/CMakeFiles/bo_hartmann.dir/depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/bo_hartmann.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/bo_hartmann.dir/flags.make

examples/CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.o: examples/CMakeFiles/bo_hartmann.dir/flags.make
examples/CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.o: examples/bo_hartmann.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zhaozhong/bayesopt/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.o"
	cd /home/zhaozhong/bayesopt/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.o -c /home/zhaozhong/bayesopt/examples/bo_hartmann.cpp

examples/CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.i"
	cd /home/zhaozhong/bayesopt/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zhaozhong/bayesopt/examples/bo_hartmann.cpp > CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.i

examples/CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.s"
	cd /home/zhaozhong/bayesopt/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zhaozhong/bayesopt/examples/bo_hartmann.cpp -o CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.s

examples/CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.o.requires:

.PHONY : examples/CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.o.requires

examples/CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.o.provides: examples/CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.o.requires
	$(MAKE) -f examples/CMakeFiles/bo_hartmann.dir/build.make examples/CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.o.provides.build
.PHONY : examples/CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.o.provides

examples/CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.o.provides.build: examples/CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.o


# Object files for target bo_hartmann
bo_hartmann_OBJECTS = \
"CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.o"

# External object files for target bo_hartmann
bo_hartmann_EXTERNAL_OBJECTS =

bin/bo_hartmann: examples/CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.o
bin/bo_hartmann: examples/CMakeFiles/bo_hartmann.dir/build.make
bin/bo_hartmann: lib/libbayesopt.a
bin/bo_hartmann: /usr/local/lib/libnlopt.so
bin/bo_hartmann: examples/CMakeFiles/bo_hartmann.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/zhaozhong/bayesopt/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/bo_hartmann"
	cd /home/zhaozhong/bayesopt/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/bo_hartmann.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/bo_hartmann.dir/build: bin/bo_hartmann

.PHONY : examples/CMakeFiles/bo_hartmann.dir/build

examples/CMakeFiles/bo_hartmann.dir/requires: examples/CMakeFiles/bo_hartmann.dir/bo_hartmann.cpp.o.requires

.PHONY : examples/CMakeFiles/bo_hartmann.dir/requires

examples/CMakeFiles/bo_hartmann.dir/clean:
	cd /home/zhaozhong/bayesopt/examples && $(CMAKE_COMMAND) -P CMakeFiles/bo_hartmann.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/bo_hartmann.dir/clean

examples/CMakeFiles/bo_hartmann.dir/depend:
	cd /home/zhaozhong/bayesopt && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zhaozhong/bayesopt /home/zhaozhong/bayesopt/examples /home/zhaozhong/bayesopt /home/zhaozhong/bayesopt/examples /home/zhaozhong/bayesopt/examples/CMakeFiles/bo_hartmann.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/bo_hartmann.dir/depend

