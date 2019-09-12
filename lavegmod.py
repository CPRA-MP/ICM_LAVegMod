#!/usr/bin/env python

##\mainpage
##
##\section Introduction
##     This is the code for the vegetation module for the Louisiana
## Master Plan Modeling 2017 project:  LAVegMod version 3.0.
##
## This model is a direct successor to LAVegMod version 1.0 that
## was used for the Master Plan Modeling 2012.
##
##\section Doc Documentation conventions
## I have used a number of conventions while documenting the
## model.
##
##\subsection CommonClassFeatures Common Class Features
## There are a number of features that are common throughout
## this model. These features are required by the Python
## programming language. These features represent coding elements
## that are needed to make a program function correctly. However,
## the do not really contain any information about how the model
## captures the ecology of wetland plant communities. I have chosen
## to omit a detailed description of these elements. This keeps the
## documentation focused on the important parts of the model without
## getting bogged down in the peculiarities of Python or programming
## in general.
##
## The coding elements that will not receive detailed descriptions are:
## - self
##   The keyword "self" occurs everywhere in python code. This variable
##   name is used as on objects reference to itself. This is a common
##   concept in object oriented programming, and is equivalent to the
##   key work "this" used in C++.
##
##   "self" is always the first argument to any class function. However,
##   only rarely does self have to be explicitly passed as an argument to
##   a member function.
##
## - "__init__"(self)
##   This is the name of a class constructor. Evey Python class has
##   __init__(self) defined, whether you define it explicitly or not.
##   I always define __init__(self) explicitly. The job of the constructor
##   is to bring an object instantiated from a class up to a defined,  initial
##   state as soon as it is created.
##
##   Sometimes you will see the constructor defined with additional arguments
##   such as __init__(self, x, y, z). In these cases the constructor takes
##   three arguments. That's all. Nothing more special than that.
##
##
##\subsection MemFnc Member functions
## For each member function I have documented the arguments
## that are passed to the function as well as the expected type
## , and where appropriate, the units of an argument. By "type"
## I mean the data type used to store values. Typical types
## include int, float, string, class, list, etc ... A variable
## of type "int" is an integer, that is a whole number with no
## fractional part. A "float" is a real number that may have a
## fractional part. That is, values to the right of the decimal
## place. So the number 1 is an integer, while 1.0 and 3.14159265 are
## both considered floats.
##
## Strings are used to contain words, such as "cat" or "Phragmites australis".
##
## Classes are python classes.
##
## List is used to indicate a python list.
##
## For more information about types please see the python documentation.
##
## Units describes the physical units that the values should be converted to
## before the values are passed to the function. For many variables, there is
## no sense of units. For example, a string such as "PHAU" does not have units.
## In these cases, the units are listed as "None".
##
##
##
##\section History
## Verison 1.0 of LAVegMod was originally developed using two programming
## languages. C++ was used for the core model that carried out the actual
## simulation of plant community dynamics. Around this core a set of R
## script was developed to covert input data into a file format that
## the C++ code could handle easily.
##
## As part of the MPM 2017 project, the LAVegMod was extensively modified
## and a number of features added. This formed version 2.0 of the LAVegMod
## A more complete discussion of these
## changes will be available in a report made to TWIG and CPRA.
##
## To facilitate a complete integration of the LAVegMod into an
## automated suite of models the LAVegMod has been translated into
## Python v2.8
##
##\section Authors Authors/Contributers
## __Model design__:
## - Jenneke Visser - University of Louisiana at Lafayette
## - Scott M. Duke-Sylvester - University of Louisiana at Lafayette
##
## __Python code__: Scott Duke-Sylvester - University of Louisiana at Lafayette
##
## __Model documentation__: Scott Duke-Sylvester - University of Louisiana at Lafayette
##
##\section Sponsors Sponsors
## This code is sponsored by the State of Louisiana as part of its on
## going efforts to perserve and manage Louisana's the unique coastal wetland ecosystems.
## Funding is provided by CPRA and the overall Master Plant Modeling project is managed by
## The Water Institute of the Gulf.
##
##
##\section Licence
## Copyright Scott M. Duke-Sylvester, Jenneke Visser, The University of Louisiana at Lafayette, The Water Institute of the Gulf, Louisiana Coastal Protection and Restoration Authority
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.



##\file
# This is the top most file in the run-time heierarchy and is the
# file that should be called from the command line to run the model.
# The command line call should look like:
# > lavegmod.py <config>
# where <config> is the name of a model configuration file.
#
#
#

import cProfile
import sys

#import bimodel
import config
import model_v2 as model
#import model_v1 as model

##\brief Program entry point
##\details This is the first function that is called when you
# run lavegmod.py
def main():

    #configDict = config.Config()
    #configDict.config(sys.argv)

    #if configDict.has_key('BarrierIslandModel'):
    #    mod = bimodel.Model()
    #elif configDict.has_key('WetlandsModel'):
    #    mod = model.Model()
    #else:
    #    mod = model.Model()

    mod = model.Model()
    return(mod.run())

#cProfile.run('main()', filename='cProfile_output.txt')
main()
