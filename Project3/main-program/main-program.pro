TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += C:\armadillo-9.600.6\include
DEPENDPATH += C:\armadillo-9.600.6\include

LIBS += \
        -larmadillo -lblas -llapack

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
    plot_gauss_legendre.py
