TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += C:\armadillo-9.600.6\include
DEPENDPATH += C:\armadillo-9.600.6\include

LIBS += \
        -larmadillo -lblas -llapack

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_
CXX_QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK

SOURCES += \
        gauss_laguerre.cpp \
        gauss_legendre.cpp \
        gauss_quadrature.cpp \
        improved_gauss_quadrature.cpp \
    improved_monte_carlo.cpp \
        main.cpp \
        monte_carlo.cpp \

HEADERS += \
    gauss_laguerre.h \
    gauss_legendre.h \
    gauss_quadrature.h \
    improved_gauss_quadrature.h \ \
    improved_monte_carlo.h \
    monte_carlo.h

DISTFILES += \
    plot_gauss_legendre.py \
    plot_mc_parallellized.py \
    plot_monte_carlo.py
