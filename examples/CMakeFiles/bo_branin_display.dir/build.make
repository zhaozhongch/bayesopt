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
include examples/CMakeFiles/bo_branin_display.dir/depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/bo_branin_display.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/bo_branin_display.dir/flags.make

examples/CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.o: examples/CMakeFiles/bo_branin_display.dir/flags.make
examples/CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.o: examples/bo_branin_display.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zhaozhong/bayesopt/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.o"
	cd /home/zhaozhong/bayesopt/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.o -c /home/zhaozhong/bayesopt/examples/bo_branin_display.cpp

examples/CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.i"
	cd /home/zhaozhong/bayesopt/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zhaozhong/bayesopt/examples/bo_branin_display.cpp > CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.i

examples/CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.s"
	cd /home/zhaozhong/bayesopt/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zhaozhong/bayesopt/examples/bo_branin_display.cpp -o CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.s

examples/CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.o.requires:

.PHONY : examples/CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.o.requires

examples/CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.o.provides: examples/CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.o.requires
	$(MAKE) -f examples/CMakeFiles/bo_branin_display.dir/build.make examples/CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.o.provides.build
.PHONY : examples/CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.o.provides

examples/CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.o.provides.build: examples/CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.o


examples/CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.o: examples/CMakeFiles/bo_branin_display.dir/flags.make
examples/CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.o: utils/displaygp.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zhaozhong/bayesopt/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object examples/CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.o"
	cd /home/zhaozhong/bayesopt/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.o -c /home/zhaozhong/bayesopt/utils/displaygp.cpp

examples/CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.i"
	cd /home/zhaozhong/bayesopt/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zhaozhong/bayesopt/utils/displaygp.cpp > CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.i

examples/CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.s"
	cd /home/zhaozhong/bayesopt/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zhaozhong/bayesopt/utils/displaygp.cpp -o CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.s

examples/CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.o.requires:

.PHONY : examples/CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.o.requires

examples/CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.o.provides: examples/CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.o.requires
	$(MAKE) -f examples/CMakeFiles/bo_branin_display.dir/build.make examples/CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.o.provides.build
.PHONY : examples/CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.o.provides

examples/CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.o.provides.build: examples/CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.o


# Object files for target bo_branin_display
bo_branin_display_OBJECTS = \
"CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.o" \
"CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.o"

# External object files for target bo_branin_display
bo_branin_display_EXTERNAL_OBJECTS =

bin/bo_branin_display: examples/CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.o
bin/bo_branin_display: examples/CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.o
bin/bo_branin_display: examples/CMakeFiles/bo_branin_display.dir/build.make
bin/bo_branin_display: lib/libbayesopt.a
bin/bo_branin_display: lib/libmatplotpp.a
bin/bo_branin_display: /usr/lib/x86_64-linux-gnu/libglut.so
bin/bo_branin_display: /usr/lib/x86_64-linux-gnu/libXmu.so
bin/bo_branin_display: /usr/lib/x86_64-linux-gnu/libXi.so
bin/bo_branin_display: /usr/lib/x86_64-linux-gnu/libGL.so
bin/bo_branin_display: /usr/lib/x86_64-linux-gnu/libGLU.so
bin/bo_branin_display: /usr/local/lib/libnlopt.so
bin/bo_branin_display: examples/CMakeFiles/bo_branin_display.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/zhaozhong/bayesopt/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable ../bin/bo_branin_display"
	cd /home/zhaozhong/bayesopt/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/bo_branin_display.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/bo_branin_display.dir/build: bin/bo_branin_display

.PHONY : examples/CMakeFiles/bo_branin_display.dir/build

examples/CMakeFiles/bo_branin_display.dir/requires: examples/CMakeFiles/bo_branin_display.dir/bo_branin_display.cpp.o.requires
examples/CMakeFiles/bo_branin_display.dir/requires: examples/CMakeFiles/bo_branin_display.dir/__/utils/displaygp.cpp.o.requires

.PHONY : examples/CMakeFiles/bo_branin_display.dir/requires

examples/CMakeFiles/bo_branin_display.dir/clean:
	cd /home/zhaozhong/bayesopt/examples && $(CMAKE_COMMAND) -P CMakeFiles/bo_branin_display.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/bo_branin_display.dir/clean

examples/CMakeFiles/bo_branin_display.dir/depend:
	cd /home/zhaozhong/bayesopt && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zhaozhong/bayesopt /home/zhaozhong/bayesopt/examples /home/zhaozhong/bayesopt /home/zhaozhong/bayesopt/examples /home/zhaozhong/bayesopt/examples/CMakeFiles/bo_branin_display.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/bo_branin_display.dir/depend

