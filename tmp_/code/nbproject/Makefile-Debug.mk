#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/data_encoder/TurtleParser.o \
	${OBJECTDIR}/src/data_encoder/parser.o \
	${OBJECTDIR}/src/data_encoder/partitioner_store.o \
	${OBJECTDIR}/src/data_encoder/profiler.o \
	${OBJECTDIR}/src/data_encoder/utils.o \
	${OBJECTDIR}/src/main.o \
	${OBJECTDIR}/src/query_parser/SPARQLLexer.o \
	${OBJECTDIR}/src/query_parser/SPARQLParser.o \
	${OBJECTDIR}/src/query_parser/node.o \
	${OBJECTDIR}/src/query_parser/query.o \
	${OBJECTDIR}/src/query_parser/query_parser.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-fopenmp
CXXFLAGS=-fopenmp

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-LSS_GraphBLAS_1.10 -Wl,-rpath,'SS_GraphBLAS_1.10'

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/code

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/code: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	g++ -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/code ${OBJECTFILES} ${LDLIBSOPTIONS} -lgraphblas -fopenmp

${OBJECTDIR}/src/data_encoder/TurtleParser.o: src/data_encoder/TurtleParser.cpp
	${MKDIR} -p ${OBJECTDIR}/src/data_encoder
	${RM} "$@.d"
	$(COMPILE.cc) -g -ISS_GraphBLAS_1.10 -Isrc/data_encoder -Isrc/magiq -Isrc/query_parser -std=c++11 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/data_encoder/TurtleParser.o src/data_encoder/TurtleParser.cpp

${OBJECTDIR}/src/data_encoder/parser.o: src/data_encoder/parser.cpp
	${MKDIR} -p ${OBJECTDIR}/src/data_encoder
	${RM} "$@.d"
	$(COMPILE.cc) -g -ISS_GraphBLAS_1.10 -Isrc/data_encoder -Isrc/magiq -Isrc/query_parser -std=c++11 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/data_encoder/parser.o src/data_encoder/parser.cpp

${OBJECTDIR}/src/data_encoder/partitioner_store.o: src/data_encoder/partitioner_store.cpp
	${MKDIR} -p ${OBJECTDIR}/src/data_encoder
	${RM} "$@.d"
	$(COMPILE.cc) -g -ISS_GraphBLAS_1.10 -Isrc/data_encoder -Isrc/magiq -Isrc/query_parser -std=c++11 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/data_encoder/partitioner_store.o src/data_encoder/partitioner_store.cpp

${OBJECTDIR}/src/data_encoder/profiler.o: src/data_encoder/profiler.cpp
	${MKDIR} -p ${OBJECTDIR}/src/data_encoder
	${RM} "$@.d"
	$(COMPILE.cc) -g -ISS_GraphBLAS_1.10 -Isrc/data_encoder -Isrc/magiq -Isrc/query_parser -std=c++11 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/data_encoder/profiler.o src/data_encoder/profiler.cpp

${OBJECTDIR}/src/data_encoder/utils.o: src/data_encoder/utils.cpp
	${MKDIR} -p ${OBJECTDIR}/src/data_encoder
	${RM} "$@.d"
	$(COMPILE.cc) -g -ISS_GraphBLAS_1.10 -Isrc/data_encoder -Isrc/magiq -Isrc/query_parser -std=c++11 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/data_encoder/utils.o src/data_encoder/utils.cpp

${OBJECTDIR}/src/main.o: src/main.cc
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -ISS_GraphBLAS_1.10 -Isrc/data_encoder -Isrc/magiq -Isrc/query_parser -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/main.o src/main.cc

${OBJECTDIR}/src/query_parser/SPARQLLexer.o: src/query_parser/SPARQLLexer.cpp
	${MKDIR} -p ${OBJECTDIR}/src/query_parser
	${RM} "$@.d"
	$(COMPILE.cc) -g -ISS_GraphBLAS_1.10 -Isrc/data_encoder -Isrc/magiq -Isrc/query_parser -std=c++11 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/query_parser/SPARQLLexer.o src/query_parser/SPARQLLexer.cpp

${OBJECTDIR}/src/query_parser/SPARQLParser.o: src/query_parser/SPARQLParser.cpp
	${MKDIR} -p ${OBJECTDIR}/src/query_parser
	${RM} "$@.d"
	$(COMPILE.cc) -g -ISS_GraphBLAS_1.10 -Isrc/data_encoder -Isrc/magiq -Isrc/query_parser -std=c++11 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/query_parser/SPARQLParser.o src/query_parser/SPARQLParser.cpp

${OBJECTDIR}/src/query_parser/node.o: src/query_parser/node.cpp
	${MKDIR} -p ${OBJECTDIR}/src/query_parser
	${RM} "$@.d"
	$(COMPILE.cc) -g -ISS_GraphBLAS_1.10 -Isrc/data_encoder -Isrc/magiq -Isrc/query_parser -std=c++11 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/query_parser/node.o src/query_parser/node.cpp

${OBJECTDIR}/src/query_parser/query.o: src/query_parser/query.cpp
	${MKDIR} -p ${OBJECTDIR}/src/query_parser
	${RM} "$@.d"
	$(COMPILE.cc) -g -ISS_GraphBLAS_1.10 -Isrc/data_encoder -Isrc/magiq -Isrc/query_parser -std=c++11 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/query_parser/query.o src/query_parser/query.cpp

${OBJECTDIR}/src/query_parser/query_parser.o: src/query_parser/query_parser.cpp
	${MKDIR} -p ${OBJECTDIR}/src/query_parser
	${RM} "$@.d"
	$(COMPILE.cc) -g -ISS_GraphBLAS_1.10 -Isrc/data_encoder -Isrc/magiq -Isrc/query_parser -std=c++11 -w -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/query_parser/query_parser.o src/query_parser/query_parser.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
