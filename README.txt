README 

The surface renewal method is utilized for estimating sensible heat flux (H) and latent heat flux (LE) from ecosystems.
To estimate H flux, measurements are taken at one height together with horizontal wind speeds and temperature.
To estimate LE fluxes, measurements may be taken at one height together with horizontal wind speeds and moisture concentrations.
The Surface renewal process is derived from understandings at an air parcel traveling a a specified height over a canopy.  In the renewal process involving heating; the air parcel moves lower into the canopy and remains for the duration tau; before another parcel replace it ejecting it upwards.

Paw et al, (1995) crated a diagram of this process and likened it to a ramp-like event. The ramp is characterised by both an amplitude (A) and period ( Greek letter :Tau). 
The purpose of this package is to calculate these two characteristics.

The sample dataset (test.csv) utilized in this package consists of measurements (grams per m3) recorded over 30 minutes (half an hour) at 20 Hz per second.
In total we have: 30 * 60 * 20 = 36000 measurements.



R code

The code implemented  begins by converting the H2O concentrations from g per m3 to mmol per m3.
Thereafter, 2nd, 3rd and 5th moments are calculated for each 1/20th of a second.

The package generates a numeric vector that consists of three values:
1) r value - this is the time (seconds) when the ratio between the third moment and observation time (seconds) is a maximum.
2) tau (Greek letter) is the total ramp duration 
3) A is the amplitude of the ramp that models the SR process.





The packages necessary for this package to work include:
NISTunits and polynom.

