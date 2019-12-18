TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += C:\armadillo-9.600.6\include
DEPENDPATH += C:\armadillo-9.600.6\include

LIBS += \
        -larmadillo -lblas -llapack

SOURCES += \
        main.cpp \
        transactions.cpp

HEADERS += \
    transactions.h

DISTFILES += \
    plot_5c.py \
    plot_rel_diff.py \
    plot_variance.py
