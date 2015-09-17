#!/bin/sh
g++ -O3 -I ./json/jsoncpp-src-0.5.0/include -o s s.cpp json_reader.cpp json_writer.cpp json_value.cpp
