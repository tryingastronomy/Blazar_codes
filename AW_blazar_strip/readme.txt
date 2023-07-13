DESCRIPTION OF CODE
The code attempts to vet out all AllWISE sources witin an uncertainty region that share IR colors consistent with blazars (really designed from FAVA flares)

This code takes either a single input (manually) or a dataframe of inputs (naming the ra, dec, and search_radius (i.e r95) columns)

It then pulls all AllWISE sources that are contained in the r95 of that detected source (initially used for g-ray detections)

First, it establishes which sources have been measured (without upper limits) and then checks if, including errors, if the AW source
falls into either the BLL or FSRQ blazar strip in all 3 Dimensions (W1-W2 vs W2-W3, W2-W3 vs W3-W4, and W1-W2 vs. W3-W4). The inputted 
source from your file (and converted to a pandas dataframe) is noted as a blazar if AT LEAST 1 AllWISE source is found to be consistent with 
the blazar strip in all 3 dimensions

Upper limits are also checked (note: the upper limit of a reported magnitude will be the SMALLEST possible value... therefore the real magnitude could be a much larger number (dimmer).)

If two or more FAVA detections share the same association, but only one has AW colors that  make it consistent with a blazar, all are removed.

We also make sure NOT to remove well-known Galactic sources if their association was detected in FAVA


Allows for output as .csv file. 
