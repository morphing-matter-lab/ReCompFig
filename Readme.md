Date: Apr. 26. 2022
written by Humphrey Yang (hanliny@cs.cmu.edu)
This file is subjected to changes.

1. Ownership

    This code is uploaded and maintained by the Morphing Matter Lab at Carnegie
    Mellon University (CMU). Please cite this work using this DOI:
    https://doi.org/10.1145/3491102.3502065
    This implementation is for academic and noncommerical uses only. Please 
    contact the author and CMU's Center for Technology Transfer and 
    Enterprise Creation if you are interested in using this software for
     commercial uses. See "LICENSE.pdf" for more information.

2. Dependencies:

    2a. System requirements
    The design tool only rans on Windows because it requires Human UI, which
    currently only supports Windows.
    
    2b. For the design tool
    Rhinoceros 3D version 7SR15 with grasshopper. Other versions may be 
    compatible with the software but are not guraranteed.
    Human UI version 0.8.1.3 for grasshopper
    GH_CPython version 0.1-alpha for grasshopper
    Python 3.9.7 (installed through Anaconda3)
        scipy 1.7.1
        numpy 1.20.3
        sympy 1.9

3. Usage instructions:

    3a. parameters.py
    This file contains the parameters used by the design tool.

    3b. Starting to use the software
    Install the dependencies. While installing Python and/or Anaconda, make 
    sure to add the executable to PATH environments. Make sure to choose the 
    correct python interpreter in GH_CPython before loading the design tool.

    3c. Using the design tool
    Run Rhinoceros 3D and grasshopper to lauch the script "UI.gh". When 
    launching the design tool, a user interface will pop up and promt the user 
    to locate the folder contraining the design tool.

    3d. Design tool tutorial
    A brief walkthrough for the design tool provided by the authors can be found
    at (coming soon). We also refer users to the paper for more information.

4. Known issues:

    The system is a work-in-progress prototype to demonstrate the algorithms. 
    The authors will maintain the implementation as much as they are available. 
    If you find any issues that are not listed here while using the tool, 
    please contact the authors.

    4a. Numerical precision
    The Rhino/grasshopper front end and the python backend uses different 
    versions of Python (2 and 3, respectively), which have very different float 
    number implementations and may sometimes leading to unexpected and erroneous 
    results during linear alebraic computations. The linear solvers may also be 
    unable to fiund solutions due to small numeric deviations. This can 
    potentially be resolved by tuning the error thresholds in parameters.py.

    4b. UI malfunctions
    The design tool may sometimes stall and become irresponsive to inputs 
    despite the users still being able to click on the command buttons. To 
    avoid this issue, make sure that after completing an action, the command 
    prompt in Rhino is not asking for further input before clicking on another 
    command. If users find themselves in such situation, press ESC to cancel 
    the ongoing command and proceed as normal. Otherwise, rerunning the UI 
    script in grasshopper may also resolve this problem, though the user may 
    have to redo the design.

    4c. Speed
    The algorithms runs in real time and the authors have optimized the 
    computations as much as possible. Yet, the current implementation is still 
    bottlenecked by the communication overhead between the front- and backend 
    python instances. I.e., the majority of the waiting time comes from the 
    software interface overhead.

    4d. Axis alignment
    The algorithm should work regardless of axis alignment. However, axis-
    aligned designs are preffered as the values and calculations are numerically 
    more stable. If the designs are not axis-aligned, the design tool will 
    generate it's own coordinate systems that best suits the design, which may 
    not be axis aligned. In this case, when interpreting the sensor 
    responsiveness diagrams, the x, y, z axis will be based on the system-
    generated coordinate system, not the model world coordinates.