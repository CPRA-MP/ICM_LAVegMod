#!/usr/bin/env python

##\file
# This file defines the Config class. The config class
# is used to read data from the model configuration
# file and present the information in a format that
# that is easy to use.

# STD Python modules
import exceptions
import re
import sys

class Config(dict):
    def __init__(self):
        dict.__init__(self)

    def config(self, argv = sys.argv):
        filename = ''
        if type(argv) == str:
           filename = argv
        elif type(argv) == list:
            if len(argv) != 2:
                errorString = 'Config: Error: There should be exactly one commandline argument, the name of the config file.'
                raise exceptions.RuntimeError(errorString)

            filename = argv[1]
        else:
            errorString = 'Config: Error: unknown type passed to Config.'
            raise exceptions.RuntimeError(errorString)


        try:
            strm = open(filename, 'r')
        except exceptions.IOError as error:
            errorString = 'Config: Error: could not open configuration file named ' + filename
            raise exceptions.RuntimeError(errorString)

        for line in strm:
            line = re.sub(r'//.*', '', line )
            line = re.sub(r'\n'  , '', line )

            while re.search(r'\\\\$', line) != None:
                contLine = strm.readline()
                contLine = re.sub(r'//.*', contLine)
                line     = re.sub(r'\\\\$', contLine, line)

            line = re.sub(r'\t', '', line)
            line = re.sub(r' ', '', line)

            if len(line) == 0:
                continue

            key, value  = line.split('=')
            dict.__setitem__(self, key, value)

        strm.close()

        if dict.has_key(self, 'XFile'):
            print 'Params: Msg: XFile requested. You don\'t have clearance.'

        return



    def __str__(self):
        ret = ''
        for (key,value) in self.iteritems():
            ret += str(key) + ' = ' + str(value) + '\n'
        return ret

