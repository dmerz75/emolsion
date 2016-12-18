#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
print (sys.version)
import time

import numpy as np
my_dir = os.path.abspath(os.path.dirname(__file__))

#  ---------------------------------------------------------  #
#  functions                                                  #
#  ---------------------------------------------------------  #
my_library = os.path.expanduser('~/.pylib')
sys.path.append(my_library)
# mpl_moving_average
# mpl_forcequench
# mpl_worm
from plot.SETTINGS import *

#  ---------------------------------------------------------  #
#  Start matplotlib (1/4)                                     #
#  ---------------------------------------------------------  #
import matplotlib
# default - Qt5Agg
# print matplotlib.rcsetup.all_backends
# matplotlib.use('GTKAgg')
# matplotlib.use('TkAgg')
print 'backend:',matplotlib.get_backend()
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
fig = plt.figure(0)

gs = GridSpec(1,1)
ax1 = plt.subplot(gs[0,:])
# ax2 = plt.subplot(gs[1,:-1])
ax = [ax1]

#  ---------------------------------------------------------  #
#  Import Data! (2/4)                                         #
#  ---------------------------------------------------------  #
result_type = 'emol' # sop | sopnucleo | gsop | namd
plot_type = 'contacts' # fe | tension | rmsd | rdf

#  ---------------------------------------------------------  #
#  mpl_myargs_begin                                           #
#  ---------------------------------------------------------  #
import argparse

def parse_arguments():
    ''' Parse script's arguments.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-o","--option",help="select None,publish,show")
    parser.add_argument("-d","--dataname",help="data name: run280, 76n")
    args = vars(parser.parse_args())
    return args


args = parse_arguments()
option = args['option']
data_name = args['dataname']


#  ---------------------------------------------------------  #
#  Import Data! (3/4)                                         #
#  ---------------------------------------------------------  #

data = np.loadtxt('emol_contacts.dat')
print data.shape


mtdata = np.reshape(data,(data.shape[0],156,data.shape[1]/156))
print mtdata.shape

# for f in mtdata.shape[0]:
# Slice:

# alphas = mtdata[::,::,0]
# betas = mtdata[::,::,1]

# for a in range(alphas.shape[1]):
#     adata = alphas[::,a]
#     plt.plot(adata)

# save_fig(my_dir,0,'fig','%s_%s_%s_alpha' % (result_type,plot_type,data_name),option)
# plt.clf()


# for a in range(betas.shape[1]):
#     adata = betas[::,a]
#     plt.plot(adata)

# save_fig(my_dir,0,'fig','%s_%s_%s_beta' % (result_type,plot_type,data_name),option)
# plt.clf()


# internals = mtdata[::,::,2]

# for a in range(internals.shape[1]):
#     adata = internals[::,a]
#     plt.plot(adata)

# save_fig(my_dir,0,'fig','%s_%s_%s_internals' % (result_type,plot_type,data_name),option)
# plt.clf()


def plot_contacts(tup):
    mdata = mtdata[::,::,tup[1]]

    for a in range(mdata.shape[1]):
        plt.plot(mdata[::,a])

    save_fig(my_dir,0,'fig','%s_%s_%s_%s' % (result_type,plot_type,data_name,tup[0]),option)
    plt.clf()


lst = [("alpha",0),("beta",1),("intra",2),("a-east",3),
       ("a-west",4),("a-south",5),("b-east",6),("b-west",7),("b-north",8)]

for l in lst:
    plot_contacts(l)


# be = mtdata[::,::,1]

# for a in range(alphas.shape[1]):
#     adata = mtdata[::,a]
#     plt.plot(adata)

# save_fig(my_dir,0,'fig','%s_%s_%s_alpha' % (result_type,plot_type,data_name),option)
# plt.clf()



#  ---------------------------------------------------------  #
#  Make final adjustments: (4/4)                              #
#  mpl - available expansions                                 #
#  ---------------------------------------------------------  #
# mpl_rc
# mpl_font
# mpl_label
# mpl_xy
# mpl_ticks
# mpl_tick
# mpl_minorticks
# mpl_legend
# combined_name = '%s_%s_%s' % (result_type, plot_type, data_name)
# save_fig



# mpl_myargs_end
