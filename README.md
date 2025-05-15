Grif-Replay
===============

A data sorting code for processing GRIFFIN-style [MIDAS](https://daq00.triumf.ca/MidasWiki/index.php/Main_Page) files into histograms.

## Architecture

Grif-Replay is composed of two parts; the data-sorting engine and server (written in c), and a set of browser-based web tools (javascript and html) for monitoring/controlling the data sorting and doing data analysis. This repository contains the c code which runs on a computer as a local or remote server. The browser code is in the [SpectrumViewer](https://github.com/GRIFFINCollaboration/SpectrumViewer) repository of the GRIFFINCollaboration GitHub.

## Installation

Download the repository from GitHub into your chosen working directory. Then in the working directory just run `make`.

### Create the Midas shared-object file

(Most users can skip this part and proceed to Running Grif-Replay. This only needs to be done on a TRIUMF DAQ computer.) If this instance of grif-replay is intended to connect to an online Midas experiment then the Midas shared-object file must also be created. This requires two additional files which are deliberately not included in this grif-replay Github repository. The two required files are named libmidas.a and midas.h. They must match the version of Midas installed on the machine that you intend to connect to (which may be a different version than the version of midas that is installed on the local machine where grif-replay is installed.). Copy these two files (libmidas.a and midas.h) to the grif-replay directory.

There is a dependancy in the midas.h file for a file named git-revision.h. This is not required by grif-replay so there are two options; either comment out this line in the midas.h file (#include "git-revision.h"), or alternatively create a blank version of this file in the grif-replay directory ('touch git-revision.h').

On some machines it might be necessary to instal the package libnsl-dev.

Then compile the Midas shared-object file with the following command:

`make midas`

See the Connect to Online instructions for how to attach grif-replay to the Midas experiment.

## Running Grif-Replay

### Launch the server

Once the code has been successfully compiled you can run it. If you will be running this as a remote server then it is a good idea to run Grif-Replay inside a screen session or similar.

Navigate to the working directory where Grif-Replay was compiled. You can launch Grif-Replay with the command:

`./grif-replay`

The output on the terminal does not need to be monitored (you can detach from the screen session now). From this point you will interact with Grif-Replay through the interface in your web browser.

### Connect the web interface

Open a web browser and navigate to the following URL:
`https://griffincollaboration.github.io/SpectrumViewer/analyzerInterface.html?backend=localhost&port=9093`

If you are connecting to a remote instance of Grif-Replay running on another computer than replace the `backend=localhost` part of this URL with the name of the computer hosting the server. For example if you want to connect to an instance of Grif-Replay running of the computer grifstore0 then the URL would be:
`https://griffincollaboration.github.io/SpectrumViewer/analyzerInterface.html?backend=grifstore0&port=9093`

### Set the paths

Set the paths to the directories of the computer running the GRIF-Replay server using the text boxes of the `Sort MIDAS files`:
MIDAS Data File Directory = the directory containing the .mid files to be sorted.
Histogram File Directory = the directory where the .tar histogram output files will be written (can be the same as MIDAS Data File Directory).
Configuration File Directory = the directory to store the configuration file with Grif-Replay settings.

### Sort MIDAS files

Ensure you are on the `Sorting Control` subpage of the Web Interface. In the table at the bottom, click on the title of the MIDAS file to be sorted so that the row will be highlighted. Click the `Submit selected files to the sorting queue`. The sorting progress can be monitored in the box at the top right of the Web Interface.

### Connect to Online

Click on the button `Connect to online MIDAS` for the appropriate server name.

### View histograms

There are a couple of ways to open the .tar histogram files to view spectra. The methods all launch the spectrum viewer with appropriate URL arguments to open directly the file. One can also open the spectrum viewer and open a file by selecting the directory and histogram file name.

On the `Sorting Control` subpage of the Web Interface click on the `Open Histo file` link in the table next to the run.

On the `Spectrum Viewer & Analysis` subpage of the Web Interface click on the title of the histogram file in the table, or click on the `Open in 2D Viewer` link in the same row.

### Convert histogram .tar files to Root format

In a terminal, navigate to the working directory, or any directory containing the tar2root.C script.
Open `grsisort` or `Root`.
run the command

`.x tar2root.C("/tig/grifstore1/grifalt/schedule146/S2232/runXXXX.tar","runXXXX.root)`

The runXXXX.root histogram file can now be opened and viewed in Root or GRSISort.

### Running without an internet connection

An internet connection is required to run the web interface and spectrumViewer codes from the GitHub servers. You can alternatively run completely offline if you launch a second server for the web codes (the first server being the grif-replay one). First download (clone) the spectrumViewer repository from GitHub into a suitable working directory. Navigate to the spectrumViewer working directory and launch a web server with a port number which is different from 9093 because that one is already in use by grif-replay. Port 9000 is used as the example here. Open a web browser and navigate to the following URL:

`localhost:9000/analyzerInterface.html?backend=localhost&port=9093`

In this example the `localhost:9000` is the spectrumViewer server available on port number 9000. The `backend=localhost&port=9093` is the grif-replay server available on port 9093.

There are a number of options for how to launch a local server, for more details (https://developer.mozilla.org/en-US/docs/Learn_web_development/Howto/Tools_and_setup/set_up_a_local_testing_server#running_a_simple_local_http_server). An easy one which is available by default on most linux distributions is to use python:

`python3 -m http.server 9000`

## Development

Grif-Replay was originally developed by Chris Pearson (data-sorting engine and server) and Adam Garnsworthy (browser tools) at TRIUMF Inc. during the 2023-2024 time period. The browser tools use much of the infrastructure originally developed by Bill Mills during the 2013-2016 time period.
