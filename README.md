# seeing_calculator
Simple program to determine the seeing of an astronomical image.

## Installation

The easiest way to install the program would be to simply git clone this rep.

```bash
git clone git@github.com:TrystanScottLambert/seeing_calculator.git
```

## Working out seeing of your image.

### Manual determination
At the moment only the manual determination is working. The given fits image is plotted and the user can interact with the matplotlib plot as normal. In particular, the user can zoom into the image and look for non-saturated stars. Once the user finds an appropriate star they can simply double click on it. This will place a red marker showing the user the position that was determined. A postage stamp is then shown in another window. If the user is satisfied that this is indeed a star and it is well centered (meaning the gaussian fit was good) then they need only double click on the postage stamp and then close the postage stamp window. The user can repeat this process until they are satisfied that there are enough stars. After that they can calculate the seeing in pixels.

Using an Ipython terminal:

```python
    run manually_determine_seeing.py    
    INFILE = '../DECAM_analysis/correct_stacks/N964/n964.fits'
    image = FitsPlot(INFILE)
    image.plot_whole() # This is where the user will manually identify stars.
    print(infile, image.calculate_seeing())
```

### Automatic determination
This will come later.