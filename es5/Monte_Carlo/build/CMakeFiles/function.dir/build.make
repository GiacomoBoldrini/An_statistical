# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.13.1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.13.1/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/boldrinicoder/uni/Analisi-Statistica/code/es5/Monte_Carlo

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/boldrinicoder/uni/Analisi-Statistica/code/es5/Monte_Carlo/build

# Include any dependencies generated for this target.
include CMakeFiles/function.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/function.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/function.dir/flags.make

CMakeFiles/function.dir/function.cpp.o: CMakeFiles/function.dir/flags.make
CMakeFiles/function.dir/function.cpp.o: ../function.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/boldrinicoder/uni/Analisi-Statistica/code/es5/Monte_Carlo/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/function.dir/function.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/function.dir/function.cpp.o -c /Users/boldrinicoder/uni/Analisi-Statistica/code/es5/Monte_Carlo/function.cpp

CMakeFiles/function.dir/function.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/function.dir/function.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/boldrinicoder/uni/Analisi-Statistica/code/es5/Monte_Carlo/function.cpp > CMakeFiles/function.dir/function.cpp.i

CMakeFiles/function.dir/function.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/function.dir/function.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/boldrinicoder/uni/Analisi-Statistica/code/es5/Monte_Carlo/function.cpp -o CMakeFiles/function.dir/function.cpp.s

# Object files for target function
function_OBJECTS = \
"CMakeFiles/function.dir/function.cpp.o"

# External object files for target function
function_EXTERNAL_OBJECTS =

libfunction.a: CMakeFiles/function.dir/function.cpp.o
libfunction.a: CMakeFiles/function.dir/build.make
libfunction.a: CMakeFiles/function.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/boldrinicoder/uni/Analisi-Statistica/code/es5/Monte_Carlo/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libfunction.a"
	$(CMAKE_COMMAND) -P CMakeFiles/function.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/function.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/function.dir/build: libfunction.a

.PHONY : CMakeFiles/function.dir/build

CMakeFiles/function.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/function.dir/cmake_clean.cmake
.PHONY : CMakeFiles/function.dir/clean

CMakeFiles/function.dir/depend:
	cd /Users/boldrinicoder/uni/Analisi-Statistica/code/es5/Monte_Carlo/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/boldrinicoder/uni/Analisi-Statistica/code/es5/Monte_Carlo /Users/boldrinicoder/uni/Analisi-Statistica/code/es5/Monte_Carlo /Users/boldrinicoder/uni/Analisi-Statistica/code/es5/Monte_Carlo/build /Users/boldrinicoder/uni/Analisi-Statistica/code/es5/Monte_Carlo/build /Users/boldrinicoder/uni/Analisi-Statistica/code/es5/Monte_Carlo/build/CMakeFiles/function.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/function.dir/depend

