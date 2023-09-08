To run extract_nirspec_ceers.py [extract_nirspec_ceers_linux.py], place the python program and the file "input_redshift.dat" in the same directory. Then, edit one line of extract_nirspec_ceers.py [extract_nirspec_ceers_linux.py] to specify the path to your 2D/ data directory (line 37 of the code).  

The program is run by typing "python extract_nirspec_ceers.py '<pointing>'" ["python extract_nirspec_ceers_linux.py], where '<pointing>' is, e.g., 'P4'. So, if you're going to look at the data in P4, type: python extract_nirspec_ceers_linux.py 'P4'

Note: instructions in brackets are for Linux, while those not in brackets are for Mac. We will be using Linux.

Keystrokes/actions in the 2D spectrum environment, showing all 3 gratings:
 In the 2D science spectrum panel, click and drag to add pixels in that x-range to the spatial profile.
 Key bindings:
 r -- "reset": Reset the 2D spectrum images.
 y -- "yposition": Toggle marker for the expected y-position of the primary target (orange dotted line).
 e -- "extract": Place the cursor over a panel and press "e" to begin 1D extraction for that grating configuration.
 f -- "fit": Identify a line or rest-frame wavelength at the pixel position of the cursor to estimate the redshift.  Once a redshift has been estimated, the rest-frame wavelength can be displayed for all panels with "w".
 w -- "wavelength": Press "w" to toggle between observed wavelength and rest-frame wavelength for the selected panel.  This feature requires a redshift estimate using "f".
 z -- "zoom": Press "z" twice in the same panel to zoom in to a range with corners defined by the cursor positions when "z" was pressed.
 x -- "axis reset": Resets the plotted range to the original range for the axis in which the cursor is.
 c -- "center": Move the point under the cursor to the center of the panel.
 i -- "in": Zoom the selected axis in, keeping the same center coordinate.
 o -- "out": Zoom the selected axis out, keeping the same center coordinate.
 q -- "quit": Exit and return to the object selection menu.
 arrow keys: Use the arrow keys to shift the plot.  The 2D spectrum panel can only be shifted left and right, while up and down changes the stretch.


Keystrokes/actions in the extraction environment for a single grating:

 In the 2D science spectrum panel, click and drag to add pixels in that x-range to the spatial profile.
 In the extracted 1D spectrum panel, click and drag to add wavelength ranges to the 1D mask.
 Key bindings:
 r -- "reset": Reset the extraction profile, redshift catalog, and 1D extraction.
 p -- "print": Print useful information including the extraction windows and current redshift catalog.
 y -- "yposition": Toggle marker for the expected y-position of the primary target (orange dotted line).
 e -- "extract": With the cursor in the extraction profile panel, press "e" twice to fit a Gaussian profile to the spatial profile.  Hit "e" once on the at the cursor y-position and extract the 1D spectrum.
 b -- "boxcar": Toggle a plot of the boxcar extraction on and off.
 f -- "fit": Fit an emission line with a Gaussian profile.  Press "f" twice, once on each side of the line, to fit the line.  The centroid guess is the highest data point in the selected wavelength range.  Follow the prompt to estimate the redshift and add the line to the redshift catalog.  Line fits currently stored in the redshift catalog will remain plotted in the 1D spectrum panel.
 m -- "modify": Modify the redshift catalog.  This function allows you to remove an entry from the redshift catalog.
 u -- "undo mask": Remove the last masked region from the 1D mask catalog.
 w -- "wavelength": Press "w" in the 2D spectrum panel to switch between pixel and observed wavelength on the x-axis.  If a redshift has been estimated from emission lines, "w" will switch between pixel, observed, and rest-frame wavelength in the 2D spectrum panel.  In the 1D spectrum panel, "w" will toggle rest-frame wavelength tick marks on the x-axis.
 z -- "zoom": Press "z" twice in the same panel to zoom in to a range with corners defined by the cursor positions when "z" was pressed.
 x -- "axis reset": Resets the plotted range to the original range for the axis in which the cursor is.
 c -- "center": Move the point under the cursor to the center of the panel.
 i -- "in": Zoom the selected axis in, keeping the same center coordinate.
 o -- "out": Zoom the selected axis out, keeping the same center coordinate.
 s -- "save": Save the extracted 1D spectrum, extraction information, and redshift catalog.
 q -- "quit": Exit the 1D extraction.  Make sure to save first!
 arrow keys: Use the arrow keys to shift the plot.  The extraction profile and 1D spectrum panels can be shifted up/down/left/right.  The 2D spectrum panel can only be shifted left and right, while up and down changes the stretch.

