# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.18.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.18.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/yumiao/Documents/Works/LAr_Sim/LAr

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/yumiao/Documents/Works/LAr_Sim/LAr/build

# Include any dependencies generated for this target.
include CMakeFiles/LAr.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/LAr.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/LAr.dir/flags.make

CMakeFiles/LAr.dir/LAr.cc.o: CMakeFiles/LAr.dir/flags.make
CMakeFiles/LAr.dir/LAr.cc.o: ../LAr.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/yumiao/Documents/Works/LAr_Sim/LAr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/LAr.dir/LAr.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LAr.dir/LAr.cc.o -c /Users/yumiao/Documents/Works/LAr_Sim/LAr/LAr.cc

CMakeFiles/LAr.dir/LAr.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LAr.dir/LAr.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/yumiao/Documents/Works/LAr_Sim/LAr/LAr.cc > CMakeFiles/LAr.dir/LAr.cc.i

CMakeFiles/LAr.dir/LAr.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LAr.dir/LAr.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/yumiao/Documents/Works/LAr_Sim/LAr/LAr.cc -o CMakeFiles/LAr.dir/LAr.cc.s

CMakeFiles/LAr.dir/src/LArActionInitialization.cc.o: CMakeFiles/LAr.dir/flags.make
CMakeFiles/LAr.dir/src/LArActionInitialization.cc.o: ../src/LArActionInitialization.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/yumiao/Documents/Works/LAr_Sim/LAr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/LAr.dir/src/LArActionInitialization.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LAr.dir/src/LArActionInitialization.cc.o -c /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArActionInitialization.cc

CMakeFiles/LAr.dir/src/LArActionInitialization.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LAr.dir/src/LArActionInitialization.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArActionInitialization.cc > CMakeFiles/LAr.dir/src/LArActionInitialization.cc.i

CMakeFiles/LAr.dir/src/LArActionInitialization.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LAr.dir/src/LArActionInitialization.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArActionInitialization.cc -o CMakeFiles/LAr.dir/src/LArActionInitialization.cc.s

CMakeFiles/LAr.dir/src/LArAnalysisManager.cc.o: CMakeFiles/LAr.dir/flags.make
CMakeFiles/LAr.dir/src/LArAnalysisManager.cc.o: ../src/LArAnalysisManager.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/yumiao/Documents/Works/LAr_Sim/LAr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/LAr.dir/src/LArAnalysisManager.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LAr.dir/src/LArAnalysisManager.cc.o -c /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArAnalysisManager.cc

CMakeFiles/LAr.dir/src/LArAnalysisManager.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LAr.dir/src/LArAnalysisManager.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArAnalysisManager.cc > CMakeFiles/LAr.dir/src/LArAnalysisManager.cc.i

CMakeFiles/LAr.dir/src/LArAnalysisManager.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LAr.dir/src/LArAnalysisManager.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArAnalysisManager.cc -o CMakeFiles/LAr.dir/src/LArAnalysisManager.cc.s

CMakeFiles/LAr.dir/src/LArAnalysisMessenger.cc.o: CMakeFiles/LAr.dir/flags.make
CMakeFiles/LAr.dir/src/LArAnalysisMessenger.cc.o: ../src/LArAnalysisMessenger.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/yumiao/Documents/Works/LAr_Sim/LAr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/LAr.dir/src/LArAnalysisMessenger.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LAr.dir/src/LArAnalysisMessenger.cc.o -c /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArAnalysisMessenger.cc

CMakeFiles/LAr.dir/src/LArAnalysisMessenger.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LAr.dir/src/LArAnalysisMessenger.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArAnalysisMessenger.cc > CMakeFiles/LAr.dir/src/LArAnalysisMessenger.cc.i

CMakeFiles/LAr.dir/src/LArAnalysisMessenger.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LAr.dir/src/LArAnalysisMessenger.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArAnalysisMessenger.cc -o CMakeFiles/LAr.dir/src/LArAnalysisMessenger.cc.s

CMakeFiles/LAr.dir/src/LArCollectorHit.cc.o: CMakeFiles/LAr.dir/flags.make
CMakeFiles/LAr.dir/src/LArCollectorHit.cc.o: ../src/LArCollectorHit.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/yumiao/Documents/Works/LAr_Sim/LAr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/LAr.dir/src/LArCollectorHit.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LAr.dir/src/LArCollectorHit.cc.o -c /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArCollectorHit.cc

CMakeFiles/LAr.dir/src/LArCollectorHit.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LAr.dir/src/LArCollectorHit.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArCollectorHit.cc > CMakeFiles/LAr.dir/src/LArCollectorHit.cc.i

CMakeFiles/LAr.dir/src/LArCollectorHit.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LAr.dir/src/LArCollectorHit.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArCollectorHit.cc -o CMakeFiles/LAr.dir/src/LArCollectorHit.cc.s

CMakeFiles/LAr.dir/src/LArCollectorSD.cc.o: CMakeFiles/LAr.dir/flags.make
CMakeFiles/LAr.dir/src/LArCollectorSD.cc.o: ../src/LArCollectorSD.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/yumiao/Documents/Works/LAr_Sim/LAr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/LAr.dir/src/LArCollectorSD.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LAr.dir/src/LArCollectorSD.cc.o -c /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArCollectorSD.cc

CMakeFiles/LAr.dir/src/LArCollectorSD.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LAr.dir/src/LArCollectorSD.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArCollectorSD.cc > CMakeFiles/LAr.dir/src/LArCollectorSD.cc.i

CMakeFiles/LAr.dir/src/LArCollectorSD.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LAr.dir/src/LArCollectorSD.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArCollectorSD.cc -o CMakeFiles/LAr.dir/src/LArCollectorSD.cc.s

CMakeFiles/LAr.dir/src/LArDetectorConstruction.cc.o: CMakeFiles/LAr.dir/flags.make
CMakeFiles/LAr.dir/src/LArDetectorConstruction.cc.o: ../src/LArDetectorConstruction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/yumiao/Documents/Works/LAr_Sim/LAr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/LAr.dir/src/LArDetectorConstruction.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LAr.dir/src/LArDetectorConstruction.cc.o -c /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArDetectorConstruction.cc

CMakeFiles/LAr.dir/src/LArDetectorConstruction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LAr.dir/src/LArDetectorConstruction.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArDetectorConstruction.cc > CMakeFiles/LAr.dir/src/LArDetectorConstruction.cc.i

CMakeFiles/LAr.dir/src/LArDetectorConstruction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LAr.dir/src/LArDetectorConstruction.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArDetectorConstruction.cc -o CMakeFiles/LAr.dir/src/LArDetectorConstruction.cc.s

CMakeFiles/LAr.dir/src/LArDetectorConstructionMessenger.cc.o: CMakeFiles/LAr.dir/flags.make
CMakeFiles/LAr.dir/src/LArDetectorConstructionMessenger.cc.o: ../src/LArDetectorConstructionMessenger.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/yumiao/Documents/Works/LAr_Sim/LAr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/LAr.dir/src/LArDetectorConstructionMessenger.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LAr.dir/src/LArDetectorConstructionMessenger.cc.o -c /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArDetectorConstructionMessenger.cc

CMakeFiles/LAr.dir/src/LArDetectorConstructionMessenger.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LAr.dir/src/LArDetectorConstructionMessenger.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArDetectorConstructionMessenger.cc > CMakeFiles/LAr.dir/src/LArDetectorConstructionMessenger.cc.i

CMakeFiles/LAr.dir/src/LArDetectorConstructionMessenger.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LAr.dir/src/LArDetectorConstructionMessenger.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArDetectorConstructionMessenger.cc -o CMakeFiles/LAr.dir/src/LArDetectorConstructionMessenger.cc.s

CMakeFiles/LAr.dir/src/LArEventAction.cc.o: CMakeFiles/LAr.dir/flags.make
CMakeFiles/LAr.dir/src/LArEventAction.cc.o: ../src/LArEventAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/yumiao/Documents/Works/LAr_Sim/LAr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/LAr.dir/src/LArEventAction.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LAr.dir/src/LArEventAction.cc.o -c /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArEventAction.cc

CMakeFiles/LAr.dir/src/LArEventAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LAr.dir/src/LArEventAction.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArEventAction.cc > CMakeFiles/LAr.dir/src/LArEventAction.cc.i

CMakeFiles/LAr.dir/src/LArEventAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LAr.dir/src/LArEventAction.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArEventAction.cc -o CMakeFiles/LAr.dir/src/LArEventAction.cc.s

CMakeFiles/LAr.dir/src/LArParticleSource.cc.o: CMakeFiles/LAr.dir/flags.make
CMakeFiles/LAr.dir/src/LArParticleSource.cc.o: ../src/LArParticleSource.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/yumiao/Documents/Works/LAr_Sim/LAr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/LAr.dir/src/LArParticleSource.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LAr.dir/src/LArParticleSource.cc.o -c /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArParticleSource.cc

CMakeFiles/LAr.dir/src/LArParticleSource.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LAr.dir/src/LArParticleSource.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArParticleSource.cc > CMakeFiles/LAr.dir/src/LArParticleSource.cc.i

CMakeFiles/LAr.dir/src/LArParticleSource.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LAr.dir/src/LArParticleSource.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArParticleSource.cc -o CMakeFiles/LAr.dir/src/LArParticleSource.cc.s

CMakeFiles/LAr.dir/src/LArParticleSourceMessenger.cc.o: CMakeFiles/LAr.dir/flags.make
CMakeFiles/LAr.dir/src/LArParticleSourceMessenger.cc.o: ../src/LArParticleSourceMessenger.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/yumiao/Documents/Works/LAr_Sim/LAr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/LAr.dir/src/LArParticleSourceMessenger.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LAr.dir/src/LArParticleSourceMessenger.cc.o -c /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArParticleSourceMessenger.cc

CMakeFiles/LAr.dir/src/LArParticleSourceMessenger.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LAr.dir/src/LArParticleSourceMessenger.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArParticleSourceMessenger.cc > CMakeFiles/LAr.dir/src/LArParticleSourceMessenger.cc.i

CMakeFiles/LAr.dir/src/LArParticleSourceMessenger.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LAr.dir/src/LArParticleSourceMessenger.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArParticleSourceMessenger.cc -o CMakeFiles/LAr.dir/src/LArParticleSourceMessenger.cc.s

CMakeFiles/LAr.dir/src/LArPhysicsList.cc.o: CMakeFiles/LAr.dir/flags.make
CMakeFiles/LAr.dir/src/LArPhysicsList.cc.o: ../src/LArPhysicsList.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/yumiao/Documents/Works/LAr_Sim/LAr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/LAr.dir/src/LArPhysicsList.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LAr.dir/src/LArPhysicsList.cc.o -c /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArPhysicsList.cc

CMakeFiles/LAr.dir/src/LArPhysicsList.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LAr.dir/src/LArPhysicsList.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArPhysicsList.cc > CMakeFiles/LAr.dir/src/LArPhysicsList.cc.i

CMakeFiles/LAr.dir/src/LArPhysicsList.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LAr.dir/src/LArPhysicsList.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArPhysicsList.cc -o CMakeFiles/LAr.dir/src/LArPhysicsList.cc.s

CMakeFiles/LAr.dir/src/LArPrimaryGeneratorAction.cc.o: CMakeFiles/LAr.dir/flags.make
CMakeFiles/LAr.dir/src/LArPrimaryGeneratorAction.cc.o: ../src/LArPrimaryGeneratorAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/yumiao/Documents/Works/LAr_Sim/LAr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/LAr.dir/src/LArPrimaryGeneratorAction.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LAr.dir/src/LArPrimaryGeneratorAction.cc.o -c /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArPrimaryGeneratorAction.cc

CMakeFiles/LAr.dir/src/LArPrimaryGeneratorAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LAr.dir/src/LArPrimaryGeneratorAction.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArPrimaryGeneratorAction.cc > CMakeFiles/LAr.dir/src/LArPrimaryGeneratorAction.cc.i

CMakeFiles/LAr.dir/src/LArPrimaryGeneratorAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LAr.dir/src/LArPrimaryGeneratorAction.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArPrimaryGeneratorAction.cc -o CMakeFiles/LAr.dir/src/LArPrimaryGeneratorAction.cc.s

CMakeFiles/LAr.dir/src/LArRunAction.cc.o: CMakeFiles/LAr.dir/flags.make
CMakeFiles/LAr.dir/src/LArRunAction.cc.o: ../src/LArRunAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/yumiao/Documents/Works/LAr_Sim/LAr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object CMakeFiles/LAr.dir/src/LArRunAction.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LAr.dir/src/LArRunAction.cc.o -c /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArRunAction.cc

CMakeFiles/LAr.dir/src/LArRunAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LAr.dir/src/LArRunAction.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArRunAction.cc > CMakeFiles/LAr.dir/src/LArRunAction.cc.i

CMakeFiles/LAr.dir/src/LArRunAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LAr.dir/src/LArRunAction.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArRunAction.cc -o CMakeFiles/LAr.dir/src/LArRunAction.cc.s

CMakeFiles/LAr.dir/src/LArSteppingAction.cc.o: CMakeFiles/LAr.dir/flags.make
CMakeFiles/LAr.dir/src/LArSteppingAction.cc.o: ../src/LArSteppingAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/yumiao/Documents/Works/LAr_Sim/LAr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building CXX object CMakeFiles/LAr.dir/src/LArSteppingAction.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LAr.dir/src/LArSteppingAction.cc.o -c /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArSteppingAction.cc

CMakeFiles/LAr.dir/src/LArSteppingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LAr.dir/src/LArSteppingAction.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArSteppingAction.cc > CMakeFiles/LAr.dir/src/LArSteppingAction.cc.i

CMakeFiles/LAr.dir/src/LArSteppingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LAr.dir/src/LArSteppingAction.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/yumiao/Documents/Works/LAr_Sim/LAr/src/LArSteppingAction.cc -o CMakeFiles/LAr.dir/src/LArSteppingAction.cc.s

# Object files for target LAr
LAr_OBJECTS = \
"CMakeFiles/LAr.dir/LAr.cc.o" \
"CMakeFiles/LAr.dir/src/LArActionInitialization.cc.o" \
"CMakeFiles/LAr.dir/src/LArAnalysisManager.cc.o" \
"CMakeFiles/LAr.dir/src/LArAnalysisMessenger.cc.o" \
"CMakeFiles/LAr.dir/src/LArCollectorHit.cc.o" \
"CMakeFiles/LAr.dir/src/LArCollectorSD.cc.o" \
"CMakeFiles/LAr.dir/src/LArDetectorConstruction.cc.o" \
"CMakeFiles/LAr.dir/src/LArDetectorConstructionMessenger.cc.o" \
"CMakeFiles/LAr.dir/src/LArEventAction.cc.o" \
"CMakeFiles/LAr.dir/src/LArParticleSource.cc.o" \
"CMakeFiles/LAr.dir/src/LArParticleSourceMessenger.cc.o" \
"CMakeFiles/LAr.dir/src/LArPhysicsList.cc.o" \
"CMakeFiles/LAr.dir/src/LArPrimaryGeneratorAction.cc.o" \
"CMakeFiles/LAr.dir/src/LArRunAction.cc.o" \
"CMakeFiles/LAr.dir/src/LArSteppingAction.cc.o"

# External object files for target LAr
LAr_EXTERNAL_OBJECTS =

LAr: CMakeFiles/LAr.dir/LAr.cc.o
LAr: CMakeFiles/LAr.dir/src/LArActionInitialization.cc.o
LAr: CMakeFiles/LAr.dir/src/LArAnalysisManager.cc.o
LAr: CMakeFiles/LAr.dir/src/LArAnalysisMessenger.cc.o
LAr: CMakeFiles/LAr.dir/src/LArCollectorHit.cc.o
LAr: CMakeFiles/LAr.dir/src/LArCollectorSD.cc.o
LAr: CMakeFiles/LAr.dir/src/LArDetectorConstruction.cc.o
LAr: CMakeFiles/LAr.dir/src/LArDetectorConstructionMessenger.cc.o
LAr: CMakeFiles/LAr.dir/src/LArEventAction.cc.o
LAr: CMakeFiles/LAr.dir/src/LArParticleSource.cc.o
LAr: CMakeFiles/LAr.dir/src/LArParticleSourceMessenger.cc.o
LAr: CMakeFiles/LAr.dir/src/LArPhysicsList.cc.o
LAr: CMakeFiles/LAr.dir/src/LArPrimaryGeneratorAction.cc.o
LAr: CMakeFiles/LAr.dir/src/LArRunAction.cc.o
LAr: CMakeFiles/LAr.dir/src/LArSteppingAction.cc.o
LAr: CMakeFiles/LAr.dir/build.make
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4Tree.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4GMocren.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4visHepRep.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4RayTracer.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4VRML.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4OpenGL.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4gl2ps.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4interfaces.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4persistency.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4error_propagation.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4readout.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4physicslists.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4parmodels.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4FR.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4vis_management.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4modeling.dylib
LAr: /Library/Developer/CommandLineTools/SDKs/MacOSX10.15.sdk/System/Library/Frameworks/OpenGL.framework/OpenGL.tbd
LAr: /Library/Developer/CommandLineTools/SDKs/MacOSX10.15.sdk/System/Library/Frameworks/OpenGL.framework/OpenGL.tbd
LAr: /opt/X11/lib/libXmu.dylib
LAr: /opt/X11/lib/libXext.dylib
LAr: /opt/X11/lib/libXt.dylib
LAr: /opt/X11/lib/libSM.dylib
LAr: /opt/X11/lib/libICE.dylib
LAr: /opt/X11/lib/libX11.dylib
LAr: /opt/X11/lib/libGLU.dylib
LAr: /opt/X11/lib/libGL.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4run.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4event.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4tracking.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4processes.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4analysis.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4zlib.dylib
LAr: /Library/Developer/CommandLineTools/SDKs/MacOSX10.15.sdk/usr/lib/libexpat.tbd
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4digits_hits.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4track.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4particles.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4geometry.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4materials.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4graphics_reps.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4intercoms.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4global.dylib
LAr: /Users/yumiao/Documents/Geant4/geant4.10.06-install/lib/libG4clhep.dylib
LAr: CMakeFiles/LAr.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/yumiao/Documents/Works/LAr_Sim/LAr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Linking CXX executable LAr"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/LAr.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/LAr.dir/build: LAr

.PHONY : CMakeFiles/LAr.dir/build

CMakeFiles/LAr.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/LAr.dir/cmake_clean.cmake
.PHONY : CMakeFiles/LAr.dir/clean

CMakeFiles/LAr.dir/depend:
	cd /Users/yumiao/Documents/Works/LAr_Sim/LAr/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/yumiao/Documents/Works/LAr_Sim/LAr /Users/yumiao/Documents/Works/LAr_Sim/LAr /Users/yumiao/Documents/Works/LAr_Sim/LAr/build /Users/yumiao/Documents/Works/LAr_Sim/LAr/build /Users/yumiao/Documents/Works/LAr_Sim/LAr/build/CMakeFiles/LAr.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/LAr.dir/depend

