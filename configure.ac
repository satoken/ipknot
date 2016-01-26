#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.15])
AC_INIT([ipknot], [0.0.3], [satoken@bio.keio.ac.jp])
AM_INIT_AUTOMAKE([-Wall foreign])
AM_CONFIG_HEADER(config.h)
AC_CONFIG_SRCDIR([src/ipknot.cpp])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_LANG(C++)

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# debug option
AC_ARG_ENABLE([debug], AS_HELP_STRING([--enable-debug], [turn on debugging (default: off)]))
if test "$enable_debug" == "yes"; then
  CFLAGS="-DDEBUG $CFLAGS"
  CXXFLAGS="-DDEBUG $CXXFLAGS"
else
  CFLAGS="-DNDEBUG $CFLAGS"
  CXXFLAGS="-DNDEBUG $CXXFLAGS"
fi

# Checks for library functions.
AC_CHECK_LIB([pthread], [pthread_create])

# glpk
AC_ARG_WITH([glpk], AS_HELP_STRING([--with-glpk(=PATH)], [support GNU linear programming kit]),
  [
    if test "$with_glpk" != "no"; then
      CPPFLAGS="-DWITH_GLPK $CPPFLAGS"
      if test "$with_glpk" != "yes"; then
        CPPFLAGS="-I${with_glpk}/include $CPPFLAGS"
      	LDFLAGS="-L${with_glpk}/lib $LDFLAGS"
      fi
      AC_CHECK_HEADERS([glpk.h])
      AC_CHECK_LIB([z], [main])
      AC_CHECK_LIB([gmp], [main])
      AC_CHECK_LIB([dl], [main])
      AC_CHECK_LIB([ltdl], [main])
      AC_CHECK_LIB([glpk], [glp_create_prob])
    fi
  ]
)

# ILOG CPLEX
AC_ARG_WITH([cplex-include],
  AC_HELP_STRING([--with-cplex-include=PATH], [location of CPLEX include directory] ),
  [CPPFLAGS="-I$withval ${CPPFLAGS}"] )
AC_ARG_WITH([cplex-lib],
  AC_HELP_STRING([--with-cplex-lib=PATH], [location of CPLEX callable library] ),
  [LDFLAGS="-L$withval ${LDFLAGS}"] )
AC_ARG_WITH([cplex], AS_HELP_STRING([--with-cplex(=PATH)], [support CPLEX optimizer]),
  [
    if test "$with_cplex" != "no"; then
      CPPFLAGS="-DWITH_CPLEX -DIL_STD $CPPFLAGS"
      if test "$with_cplex" != "yes"; then
        CPPFLAGS="-I${with_cplex}/include $CPPFLAGS"
      	LDFLAGS="-L${with_cplex}/lib $LDFLAGS"
      fi
      AC_CHECK_HEADERS([ilcplex/cplex.h])
      AC_CHECK_LIB([cplex], [main])
    fi
  ]
)

# gurobi optimizer
AC_ARG_WITH([gurobi], AS_HELP_STRING([--with-gurobi(=PATH)], [support Gurobi optimizer]),
  [
    if test "$with_gurobi" != "no"; then
      CPPFLAGS="-DWITH_GUROBI $CPPFLAGS"
      if test "$with_gurobi" != "yes"; then
        CPPFLAGS="-I${with_gurobi}/include $CPPFLAGS"
      	LDFLAGS="-L${with_gurobi}/lib $LDFLAGS"
      fi
      AC_CHECK_HEADERS([gurobi_c++.h])
      for i in `echo ${with_gurobi}/lib/libgurobi*.so | sed 's,.*/,,' | sed -e 's;lib\(gurobi.*\).so;\1;'`; do
      	  AC_CHECK_LIB([${i}], [GRBloadenv], LIBS="-lgurobi_c++ -l${i} $LIBS")
      done
    fi
  ]
)

# Vienna RNA package
AC_ARG_WITH([vienna-rna], AS_HELP_STRING([--with-vienna-rna(=PATH)], [support Vienna RNA package]),
  if test "$with_vienna_rna" != "yes" && test "$with_vienna_rna" != "no" ; then
     CPPFLAGS="-I${with_vienna_rna}/include $CPPFLAGS"
     LDFLAGS="-L${with_vienna_rna}/lib $LDFLAGS"
  fi
)
AC_CHECK_LIB([RNA], [pf_fold], [],
[
  AC_MSG_NOTICE([try '-fopenmp'])
  CXXFLAGS_BAK=${CXXFLAGS}
  CXXFLAGS="-fopenmp ${CXXFLAGS}"
  unset ac_cv_lib_RNA_pf_fold
  AC_CHECK_LIB([RNA], [pf_fold], [], [ CXXFLAGS=${CXXFLAGS_BAK} ])
])
AC_CHECK_HEADER([ViennaRNA/ribo.h],
		[AC_DEFINE([HAVE_VIENNA18], [], [have the Vienna RNA package version 1.8 or above])])
AC_CHECK_FUNC(energy_of_structure,
	      [AC_DEFINE([HAVE_VIENNA20], [], [have the Vienna RNA package version 2.0 or above])])

# getopt_long
#AX_FUNC_GETOPT_LONG

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
