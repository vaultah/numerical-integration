QT -= core
QT += widgets

TARGET = lib
TEMPLATE = lib
CONFIG += staticlib c++1z -Wall

HEADERS += const.h \
           hypercube.h \
           integrationresult.h \
           ellipsoid.h \
           polygon.h \
    normal_distribution.h \
    power_product.h \
    linear.h \
    2d_viewer.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}
