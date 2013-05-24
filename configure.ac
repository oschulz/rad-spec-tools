# Initialization

AC_INIT([Rad-Spec-Tools], [0], [], [rad-spec-tools])

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror])


# Libtool init and language selection

AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL
AC_PROG_CXX

CXX_PROJECT_INIT
AX_CXX_COMPILE_STDCXX_11([noext])

# Checks for dependencies

ROOTSYS_DEPS([glibs], [-lSpectrum])

CXX_PROJECT_DEPS([froast])


# Output

AC_CONFIG_HEADERS([config.h])

AC_CONFIG_FILES([
	Makefile
	rad-spec-tools.pc
	rad-spec-tools-config
	src/Makefile
	doc/Makefile doc/Doxyfile
])

AC_OUTPUT
