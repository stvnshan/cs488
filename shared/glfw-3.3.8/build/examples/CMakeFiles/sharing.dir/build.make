# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.27.4/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.27.4/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/build

# Include any dependencies generated for this target.
include examples/CMakeFiles/sharing.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include examples/CMakeFiles/sharing.dir/compiler_depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/sharing.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/sharing.dir/flags.make

examples/sharing.app/Contents/Resources/glfw.icns: /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/examples/glfw.icns
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Copying OS X content examples/sharing.app/Contents/Resources/glfw.icns"
	$(CMAKE_COMMAND) -E copy /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/examples/glfw.icns examples/sharing.app/Contents/Resources/glfw.icns

examples/CMakeFiles/sharing.dir/sharing.c.o: examples/CMakeFiles/sharing.dir/flags.make
examples/CMakeFiles/sharing.dir/sharing.c.o: /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/examples/sharing.c
examples/CMakeFiles/sharing.dir/sharing.c.o: examples/CMakeFiles/sharing.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object examples/CMakeFiles/sharing.dir/sharing.c.o"
	cd /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/build/examples && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT examples/CMakeFiles/sharing.dir/sharing.c.o -MF CMakeFiles/sharing.dir/sharing.c.o.d -o CMakeFiles/sharing.dir/sharing.c.o -c /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/examples/sharing.c

examples/CMakeFiles/sharing.dir/sharing.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/sharing.dir/sharing.c.i"
	cd /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/build/examples && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/examples/sharing.c > CMakeFiles/sharing.dir/sharing.c.i

examples/CMakeFiles/sharing.dir/sharing.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/sharing.dir/sharing.c.s"
	cd /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/build/examples && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/examples/sharing.c -o CMakeFiles/sharing.dir/sharing.c.s

examples/CMakeFiles/sharing.dir/__/deps/glad_gl.c.o: examples/CMakeFiles/sharing.dir/flags.make
examples/CMakeFiles/sharing.dir/__/deps/glad_gl.c.o: /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/deps/glad_gl.c
examples/CMakeFiles/sharing.dir/__/deps/glad_gl.c.o: examples/CMakeFiles/sharing.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object examples/CMakeFiles/sharing.dir/__/deps/glad_gl.c.o"
	cd /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/build/examples && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT examples/CMakeFiles/sharing.dir/__/deps/glad_gl.c.o -MF CMakeFiles/sharing.dir/__/deps/glad_gl.c.o.d -o CMakeFiles/sharing.dir/__/deps/glad_gl.c.o -c /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/deps/glad_gl.c

examples/CMakeFiles/sharing.dir/__/deps/glad_gl.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/sharing.dir/__/deps/glad_gl.c.i"
	cd /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/build/examples && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/deps/glad_gl.c > CMakeFiles/sharing.dir/__/deps/glad_gl.c.i

examples/CMakeFiles/sharing.dir/__/deps/glad_gl.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/sharing.dir/__/deps/glad_gl.c.s"
	cd /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/build/examples && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/deps/glad_gl.c -o CMakeFiles/sharing.dir/__/deps/glad_gl.c.s

# Object files for target sharing
sharing_OBJECTS = \
"CMakeFiles/sharing.dir/sharing.c.o" \
"CMakeFiles/sharing.dir/__/deps/glad_gl.c.o"

# External object files for target sharing
sharing_EXTERNAL_OBJECTS =

examples/sharing.app/Contents/MacOS/sharing: examples/CMakeFiles/sharing.dir/sharing.c.o
examples/sharing.app/Contents/MacOS/sharing: examples/CMakeFiles/sharing.dir/__/deps/glad_gl.c.o
examples/sharing.app/Contents/MacOS/sharing: examples/CMakeFiles/sharing.dir/build.make
examples/sharing.app/Contents/MacOS/sharing: src/libglfw3.a
examples/sharing.app/Contents/MacOS/sharing: examples/CMakeFiles/sharing.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable sharing.app/Contents/MacOS/sharing"
	cd /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/build/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sharing.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/sharing.dir/build: examples/sharing.app/Contents/MacOS/sharing
examples/CMakeFiles/sharing.dir/build: examples/sharing.app/Contents/Resources/glfw.icns
.PHONY : examples/CMakeFiles/sharing.dir/build

examples/CMakeFiles/sharing.dir/clean:
	cd /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/build/examples && $(CMAKE_COMMAND) -P CMakeFiles/sharing.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/sharing.dir/clean

examples/CMakeFiles/sharing.dir/depend:
	cd /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8 /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/examples /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/build /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/build/examples /Users/stevenshan/Documents/cs488/shared/glfw-3.3.8/build/examples/CMakeFiles/sharing.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : examples/CMakeFiles/sharing.dir/depend

