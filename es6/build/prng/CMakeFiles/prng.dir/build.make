# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Produce verbose output by default.
VERBOSE = 1

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
CMAKE_SOURCE_DIR = /home/wahid/Wahid/University/analisi_statistica/esercizi/es6

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/wahid/Wahid/University/analisi_statistica/esercizi/es6/build

# Include any dependencies generated for this target.
include prng/CMakeFiles/prng.dir/depend.make

# Include the progress variables for this target.
include prng/CMakeFiles/prng.dir/progress.make

# Include the compile flags for this target's objects.
include prng/CMakeFiles/prng.dir/flags.make

prng/CMakeFiles/prng.dir/prng.cpp.o: prng/CMakeFiles/prng.dir/flags.make
prng/CMakeFiles/prng.dir/prng.cpp.o: ../prng/prng.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/wahid/Wahid/University/analisi_statistica/esercizi/es6/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object prng/CMakeFiles/prng.dir/prng.cpp.o"
	cd /home/wahid/Wahid/University/analisi_statistica/esercizi/es6/build/prng && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/prng.dir/prng.cpp.o -c /home/wahid/Wahid/University/analisi_statistica/esercizi/es6/prng/prng.cpp

prng/CMakeFiles/prng.dir/prng.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/prng.dir/prng.cpp.i"
	cd /home/wahid/Wahid/University/analisi_statistica/esercizi/es6/build/prng && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/wahid/Wahid/University/analisi_statistica/esercizi/es6/prng/prng.cpp > CMakeFiles/prng.dir/prng.cpp.i

prng/CMakeFiles/prng.dir/prng.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/prng.dir/prng.cpp.s"
	cd /home/wahid/Wahid/University/analisi_statistica/esercizi/es6/build/prng && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/wahid/Wahid/University/analisi_statistica/esercizi/es6/prng/prng.cpp -o CMakeFiles/prng.dir/prng.cpp.s

# Object files for target prng
prng_OBJECTS = \
"CMakeFiles/prng.dir/prng.cpp.o"

# External object files for target prng
prng_EXTERNAL_OBJECTS =

prng/libprng.a: prng/CMakeFiles/prng.dir/prng.cpp.o
prng/libprng.a: prng/CMakeFiles/prng.dir/build.make
prng/libprng.a: prng/CMakeFiles/prng.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/wahid/Wahid/University/analisi_statistica/esercizi/es6/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libprng.a"
	cd /home/wahid/Wahid/University/analisi_statistica/esercizi/es6/build/prng && $(CMAKE_COMMAND) -P CMakeFiles/prng.dir/cmake_clean_target.cmake
	cd /home/wahid/Wahid/University/analisi_statistica/esercizi/es6/build/prng && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/prng.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
prng/CMakeFiles/prng.dir/build: prng/libprng.a

.PHONY : prng/CMakeFiles/prng.dir/build

prng/CMakeFiles/prng.dir/clean:
	cd /home/wahid/Wahid/University/analisi_statistica/esercizi/es6/build/prng && $(CMAKE_COMMAND) -P CMakeFiles/prng.dir/cmake_clean.cmake
.PHONY : prng/CMakeFiles/prng.dir/clean

prng/CMakeFiles/prng.dir/depend:
	cd /home/wahid/Wahid/University/analisi_statistica/esercizi/es6/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/wahid/Wahid/University/analisi_statistica/esercizi/es6 /home/wahid/Wahid/University/analisi_statistica/esercizi/es6/prng /home/wahid/Wahid/University/analisi_statistica/esercizi/es6/build /home/wahid/Wahid/University/analisi_statistica/esercizi/es6/build/prng /home/wahid/Wahid/University/analisi_statistica/esercizi/es6/build/prng/CMakeFiles/prng.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : prng/CMakeFiles/prng.dir/depend
