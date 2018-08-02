TEMPLATE = app
CONFIG += console c++1z -Wall
CONFIG -= app_bundle
QT += widgets

DEPENDPATH += . ../lib
INCLUDEPATH += ../lib
LIBS += -L../lib

SOURCES += main.cpp
