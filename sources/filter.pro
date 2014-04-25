#-------------------------------------------------
#
# Project created by QtCreator 2014-04-25T08:55:43
#
#-------------------------------------------------

QT       -= gui

TARGET = filter
TEMPLATE = lib

INCLUDEPATH += 'J:/Program Files (x86)/Embarcadero/RAD Studio/12.0/include/boost_1_50'



    DESTDIR = dist

DEFINES += FILTER_LIBRARY

SOURCES += filter.cpp

HEADERS += filter.h\
        filter_global.h

    CONFIG += build_all

unix {
    target.path = /usr/lib
    INSTALLS += target
}
