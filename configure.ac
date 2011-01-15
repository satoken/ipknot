#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.15])
AC_INIT([ipknot], [0.0.1], [satoken@k.u-tokyo.ac.jp])
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

# Checks for library functions.
AC_CHECK_LIB([pthread], [pthread_create])

# glpk
AC_ARG_WITH([glpk], AS_HELP_STRING([--with-glpk], [support GNU linear programming kit]),
  if test "$with_glpk" != "no"; then
    CPPFLAGS="-DWITH_GLPK $CPPFLAGS"
  fi		  
)
AC_CHECK_HEADERS([glpk.h])
AC_CHECK_LIB([z], [main])
AC_CHECK_LIB([gmp], [main])
AC_CHECK_LIB([dl], [main])
AC_CHECK_LIB([ltdl], [main])
AC_CHECK_LIB([glpk], [glp_create_prob])

# ILOG CPLEX
AC_ARG_WITH([cplex], AS_HELP_STRING([--with-cplex], [support CPLEX optimizer]),
  if test "$with_cplex" != "no"; then
    CPPFLAGS="-DWITH_CPLEX -DIL_STD $CPPFLAGS"
  fi		  
)
AC_CHECK_HEADERS([ilconcert/ilomodel.h])
AC_CHECK_LIB([concert], [main])
AC_CHECK_HEADERS([ilcplex/ilocplex.h])
AC_CHECK_LIB([cplex], [main])
AC_CHECK_LIB([ilocplex], [main])

# gurobi optimizer
AC_ARG_WITH([gurobi], AS_HELP_STRING([--with-gurobi], [support Gurobi optimizer]),
  if test "$with_gurobi" != "no"; then
    CPPFLAGS="-DWITH_GUROBI $CPPFLAGS"
    if test "$with_gurobi" != "yes"; then
      CPPFLAGS="-I${with_gurobi}/include $CPPFLAGS"
      LDFLAGS="-L${with_gurobi}/lib $LDFLAGS"
    fi
  fi
)
AC_CHECK_HEADERS([gurobi_c++.h])
AC_CHECK_LIB([gurobi], [GRBloadenv], LIBS="-lgurobi_c++ -lgurobi $LIBS")
AC_CHECK_LIB([gurobi20], [GRBloadenv], LIBS="-lgurobi_c++ -lgurobi20 $LIBS")
AC_CHECK_LIB([gurobi30], [GRBloadenv], LIBS="-lgurobi_c++ -lgurobi30 $LIBS")
AC_CHECK_LIB([gurobi40], [GRBloadenv], LIBS="-lgurobi_c++ -lgurobi40 $LIBS")

# Vienna RNA package
AC_ARG_WITH([vienna-rna], AS_HELP_STRING([--with-vienna-rna], [support Vienna RNA package]),
  if test "$with_vienna_rna" != "yes" && test "$with_vienna_rna" != "no" ; then
     CPPFLAGS="-I${with_vienna_rna}/include $CPPFLAGS"
     LDFLAGS="-L${with_vienna_rna}/lib $LDFLAGS"
  fi
)
AC_CHECK_LIB(RNA, init_pf_fold)
AC_CHECK_HEADER([ViennaRNA/ribo.h],
		[AC_DEFINE([HAVE_VIENNA18], [], [have the Vienna RNA package version 1.8 or above])])


AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
