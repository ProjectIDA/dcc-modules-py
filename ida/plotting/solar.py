#######################################################################################
# Copyright (C) 2016, 2020  Regents of the University of California
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
from matplotlib.dates import date2num, DateFormatter
import matplotlib.dates as mdates
from matplotlib.ticker import (AutoMinorLocator)
import numpy as np
from obspy import read
from ida.utils import msget


SOLAR_SOH_PLOT_CHANLOCS = [
    {'chn': 'aev', 'loc': '02', 'gain': 100.0,
     'xlabel': 'Load Voltage', 'ylabel': 'AEV-02  (V)'},
    {'chn': 'aec', 'loc': '02', 'gain': 100.0,
     'xlabel': 'Load Current', 'ylabel': 'AEC-02  (A)'},
    {'chn': 'aec', 'loc': '01', 'gain': 100.0,
     'xlabel': 'Charge Current', 'ylabel': 'AEC-01  (A)'},
    {'chn': 'ae2', 'loc': '01', 'gain': 100.0,
     'xlabel': 'Array Voltage', 'ylabel': 'AE2-01  (V)'},
]


def get_solar_station_list():
    """Return list os IDA station codes for solar powered stations

    Returns:
        list of station codes: [str]
    """

    # TODO: replace with call to ida station api endpoint
    sta_list = [
        'aak',
        'abpo',
        'efi',
        'kapi',
        'mbar',
        'sacv',
        'wrab',
    ]
    return sta_list


def gen_station_solar_soh_fig(sta, start_dt, end_dt):
    """

    Parameters
    ----------
    sta: str
        Station code of solar station to generate SOH plots for
    start_dt: datetime.datetime
        start datetime pof period to plot
    end_dt: datetime.datetime
        end datetime of period to plot

    Returns
    -------
    Plot figure: matploib.pyplot.Figure
        Up to caller to decide what to do with the Figure. Save to disk or show interactively.
    """

    fig = plt.figure(figsize=(8.5, 11))

    fig.suptitle(
        "{} - SPS State of Health\n{}  -  {}\n(updated {} UTC)".format(
            sta.upper(),
            start_dt.strftime("%Y-%m-%d %H:%M:%S"), end_dt.strftime("%Y-%m-%d %H:%M:%S"),
            datetime.datetime.utcnow().strftime("%Y-%m-%d (%j) %H:%M")
        ),
        fontweight='bold'
    )

    # this is how many plots in this figure.
    plot_cnt = len(SOLAR_SOH_PLOT_CHANLOCS)

    for ndx, chan_info in enumerate(SOLAR_SOH_PLOT_CHANLOCS):

        sta = sta.lower()
        chan = chan_info['chn'].lower()
        loc = chan_info['loc'].lower()

        arc_dir = os.environ.get('IDA_ARCHIVE_MS_DIR', '/ida/archive/ms')

        # send trimmed output to tempfile for pipping to stdout
        outtf = tempfile.NamedTemporaryFile(delete=False)

        msget(arc_dir, sta, chan, loc, start_dt, end_dt, outfn=outtf.name)

        stream = read(outtf.name)

        os.remove(outtf.name)

        # Note all plots in a single column
        ax = fig.add_subplot("{}1{}".format(str(plot_cnt), str(ndx + 1)))

        ax.set_ylabel(chan_info['ylabel'], fontweight='bold')
        ax.set_xlabel(chan_info['xlabel'], fontweight='bold')
        ax.set_xlim(start_dt, end_dt)

        # set tick and label properties
        ax.xaxis.set_major_formatter(DateFormatter("%Y-%m-%d"))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=7))
        ax.xaxis.set_minor_locator(mdates.DayLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
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

        for tr in stream:
            data = np.require(tr.data, np.float64) / chan_info['gain']
            times = ((tr.times() / (60 * 60 * 24)) +
                     date2num(tr.stats.starttime.datetime))

            ax.plot(times, data, 'b-', linewidth=0.5, linestyle='-')

        # set final limits for right-side Y-axis so
        # ticks and lables in the same place as left-side
        ybot, ytop = ax.get_ylim()
        ax2.set_ylim(bottom=ybot, top=ytop)

    # firce Y-lables to align hortizontally
    # and tweak spacing around plots
    fig.align_ylabels()
    plt.subplots_adjust(top=0.92, bottom=0.06, hspace=0.4)

    return fig

