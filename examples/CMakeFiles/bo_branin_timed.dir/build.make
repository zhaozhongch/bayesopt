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
include examples/CMakeFiles/bo_branin_timed.dir/depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/bo_branin_timed.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/bo_branin_timed.dir/flags.make

examples/CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.o: examples/CMakeFiles/bo_branin_timed.dir/flags.make
examples/CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.o: examples/bo_branin_timed.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zhaozhong/bayesopt/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.o"
	cd /home/zhaozhong/bayesopt/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.o -c /home/zhaozhong/bayesopt/examples/bo_branin_timed.cpp

examples/CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.i"
	cd /home/zhaozhong/bayesopt/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zhaozhong/bayesopt/examples/bo_branin_timed.cpp > CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.i

examples/CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.s"
	cd /home/zhaozhong/bayesopt/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zhaozhong/bayesopt/examples/bo_branin_timed.cpp -o CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.s

examples/CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.o.requires:

.PHONY : examples/CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.o.requires

examples/CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.o.provides: examples/CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.o.requires
	$(MAKE) -f examples/CMakeFiles/bo_branin_timed.dir/build.make examples/CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.o.provides.build
.PHONY : examples/CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.o.provides

examples/CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.o.provides.build: examples/CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.o


# Object files for target bo_branin_timed
bo_branin_timed_OBJECTS = \
"CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.o"

# External object files for target bo_branin_timed
bo_branin_timed_EXTERNAL_OBJECTS =

bin/bo_branin_timed: examples/CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.o
bin/bo_branin_timed: examples/CMakeFiles/bo_branin_timed.dir/build.make
bin/bo_branin_timed: lib/libbayesopt.a
bin/bo_branin_timed: /usr/local/lib/libnlopt.so
bin/bo_branin_timed: examples/CMakeFiles/bo_branin_timed.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/zhaozhong/bayesopt/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/bo_branin_timed"
	cd /home/zhaozhong/bayesopt/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/bo_branin_timed.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/bo_branin_timed.dir/build: bin/bo_branin_timed

.PHONY : examples/CMakeFiles/bo_branin_timed.dir/build

examples/CMakeFiles/bo_branin_timed.dir/requires: examples/CMakeFiles/bo_branin_timed.dir/bo_branin_timed.cpp.o.requires

.PHONY : examples/CMakeFiles/bo_branin_timed.dir/requires

examples/CMakeFiles/bo_branin_timed.dir/clean:
	cd /home/zhaozhong/bayesopt/examples && $(CMAKE_COMMAND) -P CMakeFiles/bo_branin_timed.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/bo_branin_timed.dir/clean

examples/CMakeFiles/bo_branin_timed.dir/depend:
	cd /home/zhaozhong/bayesopt && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zhaozhong/bayesopt /home/zhaozhong/bayesopt/examples /home/zhaozhong/bayesopt /home/zhaozhong/bayesopt/examples /home/zhaozhong/bayesopt/examples/CMakeFiles/bo_branin_timed.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/bo_branin_timed.dir/depend

