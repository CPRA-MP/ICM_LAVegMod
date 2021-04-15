#!/usr/bin/env python

##\file
# This file defined the base Event class and the EventQueue class.
# These classes form the basis of the model update process

import copy
import heapq
import signal
import sys

class Time(object):
    def __init__(self, time, priority):
        self.time     = copy.copy(time)
        self.priority = copy.copy(priority)

    def __lt__(self, other):
        return self.time < other.time or ( self.time == other.time and self.priority < other.priority )

    def __str__(self):
        return str(self.time) + '(' + str(self.priority) + ')'

class Event(object):

    def __init__(self, time, name = 'Event' ):
        self.time = copy.copy(time)
        self.name = copy.copy(name)

    def __lt__(self, other):
        return self.time < other.time

    def __str__(self):
        return str(self.time) + ' : ' + self.name

    def act(self):
        print((self.__str__()))

class GenericEvent(Event):
    def __init__(self, time, name = 'GenericEvent', callable=None):
        Event.__init__(self, time, name)
        self.callable = callable

    def act(self):
        print((self.name))
        self.callable.act()

class PauseEvent(Event):
    def __init__(self, time, name='PauseEvent'):
        Event.__init__(self, time, name)

    def act(self):
        pass

class MsgEvent(Event):

    def __init__(self, time, name = 'MsgEvent', stream=sys.stdout, msg = 'None'):
        Event.__init__(self, time, name)
        self.stream = stream
        self.msg = copy.copy(msg)

    def __str__(self):
        return self.msg

    def act(self):
        print(self.__str__(), file=self.stream)

#class ReadCSVEvent(Event):
#    def __init__(self, time, name = 'ReadCSVEvent', file=0):
#        Event.__init__(self, time, name)
#        self.file = copy.copy(file)
#
#    def __str__(self):
#        return Event.__str__(self)
#
#    def read(self, file):
#        pass
#
#    def act(self):
#        pass

class EventQueue(object):

    def __init__(self):
        self.queue = list()

    def run(self, while_condition=(lambda arg: True)):
        elt = Event(Time(0,0))
        while len(self.queue) != 0 and while_condition(elt):
            elt = heapq.heappop(self.queue)
            elt.act()

    def add_event(self, event):
        heapq.heappush(self.queue, event)

    def clear(self):
        del self.queue
        self.queue = list()

    def __str__(self):
        ret = ''
        for elt in self.queue:
            ret += str(elt) + '\n'
        return ret
