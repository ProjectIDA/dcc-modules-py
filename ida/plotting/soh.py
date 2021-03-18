#######################################################################################
# Copyright (C) 2016, 2021  Regents of the University of California
#
# This is free software: you can redistribute it and/or modify it under the terms of
# the GNU General Public License (GNU GPL) as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# A copy of the GNU General Public License can be found in LICENSE.TXT in the root of
# the source code repository. It can also be found at http://www.gnu.org/licenses/.
#
# NOTES: Per GNU GPLv3 terms:
#   * This notice must be kept in this source file
#   * Changes to the source must be clearly noted with date & time of change
#
# If you use this software in a product, an explicit acknowledgment in the product
# documentation of the contribution by Project IDA, Institute of Geophysics and
# Planetary Physics, UCSD would be appreciated but is not required.
#######################################################################################
import tempfile
import datetime
import os
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, DateFormatter, DAILY
import matplotlib.dates as mdates
from matplotlib.ticker import AutoMinorLocator
import numpy as np
from obspy import read
from ida.utils import msget


RPM_SOH_PLOT_CHAN_INFO = [
    {'chn': 'rl1', 'loc': '25', 'gain': 1,
     'xlabel': 'Relay 1 State (0 => open; 1 => closed)', 'ylabel': 'RL1-25 (0/1)'},
    {'chn': 'rl2', 'loc': '25', 'gain': 1,
     'xlabel': 'Relay 2 State (0 => open; 1 => closed)', 'ylabel': 'RL2-25 (0/1)'},
    {'chn': 'tpe', 'loc': '25', 'gain': 10.0,
     'xlabel': 'External Temp', 'ylabel': 'TPE-25  (Deg C)'},
    {'chn': 'tpi', 'loc': '25', 'gain': 10.0,
     'xlabel': 'Internal Temp', 'ylabel': 'TPI-25  (Deg C)'},
    {'chn': 'mv1', 'loc': '25', 'gain': 10.0,
     'xlabel': 'Battery O/P Voltage', 'ylabel': 'MV1-25  (V)'},
    {'chn': 'mv2', 'loc': '25', 'gain': 10.0,
     'xlabel': 'Newmar O/P Voltage', 'ylabel': 'MV2-25  (V)'},
    {'chn': 'mv3', 'loc': '25', 'gain': 10.0,
     'xlabel': 'Wilmore O/P Voltage', 'ylabel': 'MV3-25  (V)'},
    {'chn': 'mv4', 'loc': '25', 'gain': 10.0,
     'xlabel': 'AC Indicator (0 => ON; 12 => OFF)', 'ylabel': 'MV4-25'},
    {'chn': 'mc1', 'loc': '25', 'gain': 10.0,
     'xlabel': 'Wilmore O/P Current', 'ylabel': 'MC1-25  (A)'},
    {'chn': 'mc2', 'loc': '25', 'gain': 10.0,
     'xlabel': 'Vault O/P Current', 'ylabel': 'MC2-25  (A)'},
    {'chn': 'mc3', 'loc': '25', 'gain': 10.0,
     'xlabel': 'Shunt Battery Current', 'ylabel': 'MC3-25  (A)'},
]

class plotchan:

    def __init__(self, chn, loc, gain, xlabel, ylabel):

        self.chn = chn
        self.loc = loc
        self.gain = gain
        self.xlabel = xlabel
        self.ylabel = ylabel


def gen_station_soh_fig(fig_title, plot_chans, sta, start_dt, end_dt, plots_per_fig = 4):
    """

    Parameters
    ----------
    fig_title: str
        title text to be integrated at top of fig
    plot_chans: list of plotchan with channel info
    sta: str
        Station code of solar station to generate SOH plots for
    start_dt: datetime.datetime
        start datetime pof period to plot
    end_dt: datetime.datetime
        end datetime of period to plot

    plots_per_fig (optional): int
        number of plots per mnatplotlib figure

    Returns
    -------
    List of matploib.pyplot.Figure
        Up to caller to decide what to do with the Figures. Save to disk or show interactively.


    Produces plots from miniseed archive that are optimized to 7-30 day periods.

"""

    PLOTS_PER_FIG = plots_per_fig
    plot_cnt = len(plot_chans)
    fig_list = []

    for plot_ndx, chan_info in enumerate(plot_chans):

        if plot_ndx % PLOTS_PER_FIG == 0:

            fig = plt.figure(figsize=(8.5, 11))

            fig.suptitle(
                "{} - {}\n{}  -  {}\n(updated {} UTC)".format(
                    sta.upper(), fig_title,
                    start_dt.strftime("%Y-%m-%d %H:%M:%S"), end_dt.strftime("%Y-%m-%d %H:%M:%S"),
                    datetime.datetime.utcnow().strftime("%Y-%m-%d (%j) %H:%M")
                ),
                fontweight='bold'
            )

            sub_plot_ndx = 1

        else:
            sub_plot_ndx += 1

        sta = sta.lower()
        chan = chan_info.chn.lower()
        loc = chan_info.loc.lower()

        arc_dir = os.environ.get('IDA_ARCHIVE_MS_DIR', '/ida/archive/ms')

        # send trimmed output to tempfile for pipping to stdout
        outtf = tempfile.NamedTemporaryFile(delete=False)

        msget(arc_dir, sta, chan, loc, start_dt, end_dt, outfn=outtf.name)

        stream = read(outtf.name)

        os.remove(outtf.name)

        # Note all plots in a single column
        ax = fig.add_subplot("{}1{}".format(str(PLOTS_PER_FIG), str(sub_plot_ndx)))

        ax.set_ylabel(chan_info.ylabel, fontweight='bold')
        ax.set_xlabel(chan_info.xlabel, fontweight='bold')
        ax.set_xlim(start_dt, end_dt)

        # set tick and label properties
        ax.xaxis.set_major_formatter(DateFormatter("%Y-%m-%d"))
        maj_loctr = mdates.AutoDateLocator(minticks=3, maxticks=9, interval_multiples=False)
        maj_loctr.intervald[DAILY] = [2,3,5,7,10,14,21,28]
        ax.xaxis.set_major_locator(maj_loctr)
        ax.xaxis.set_minor_locator(mdates.DayLocator())

        ax.grid(True, linestyle='dotted')

        # set up a second axis so we can show also Y-axis
        # ticks and labels on right side of plot
        ax2 = ax.twinx()
        ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax2.set_ylabel(ax.get_ylabel(), labelpad=12.0,
                       fontweight='bold', rotation=270)
        ax2.tick_params(axis='y', which='both',
                        left=False, right=True,
                        labelleft=False, labelright=True)

        plt.xticks(fontsize=7)

        # loop through traces in stream add data to plot.
        # NOTE: discontinuous traces will appear coninuous on plot
        for tr in stream:
            data = np.require(tr.data, np.float64) / chan_info.gain
            times = ((tr.times() / (60 * 60 * 24)) +
                     date2num(tr.stats.starttime.datetime))

            ax.plot(times, data, 'b-', linewidth=0.5, linestyle='-')

        # set final limits for right-side Y-axis so
        # ticks and lables in the same place as left-side
        ybot, ytop = ax.get_ylim()
        ax2.set_ylim(bottom=ybot, top=ytop)

        # IF last plot on fig do some final adjustments to fig.

        if (sub_plot_ndx == PLOTS_PER_FIG) or (plot_ndx == (plot_cnt- 1)):

            # force Y-lables to align hortizontally
            # and tweak spacing around plots
            fig.align_ylabels()
            plt.subplots_adjust(top=0.92, bottom=0.06, hspace=0.4)

            fig_list.append(fig)

    return fig_list

