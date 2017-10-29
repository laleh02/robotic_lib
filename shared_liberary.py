"""
MOONS OPS to FPS Shared Library prototype

A prototype module containing functions shared between the MOONS
Observation Preparation Software, Fibre Positioner Control Software
and Fibre Positioner Collision Avoidance Software

This prototype is written in Python. It is expected that the final version
will be translated into C or C++ for greater compatibility with ESO/VLT standard.

08 Jul 2015: First version as fibre positioner simulator.

17 Feb 2016: check_targets function reworked to choose the best parity for
             each positioner. Low-level functions are now check_pair and
             test_pair.
24 Mar 2016: Fail gracefully if the plotting module cannot be imported.
10 May 2016: Improved logging. Make sure angles are reported consistently
             in degrees.
28 Jun 2016: Temporarily added Motor code from path analysis, until the
             path analysis code is successfully merged.

26 Jul 2016: Added set_arm_angles() function, so that the starting angles
             can be defined as parameters. Added weights parameter to
             difficulty_values function. Adjusted avoidance zone sizes and
             removed the circular zone at the end of the datum foot.
27 Jul 2016: target_angles() renamed to get_arm_angles(). Implemented the
             alpha zero-point angle, self.orient, for each positioner.
             Added arm motor tests.


"""
# Python 3 emulation under Python 2.7
from __future__ import absolute_import, unicode_literals, division, print_function

# Import Python logging and select a logging level.
import logging
# logging.basicConfig(level=logging.ERROR)  # Errors only
# logging.basicConfig(level=logging.WARN)   # Warnings and errors only
logging.basicConfig(level=logging.INFO)   # Informational output
# logging.basicConfig(level=logging.DEBUG)  # Some debugging output
logger = logging.getLogger("mocpath.fps_shared")

# Set True to report conflict table after each pass.
EXTRA_CONFLICT_INFO = True

# Use this flag to turn graphical debugging on and off.
# When this flag is True, plots of conflicting positioners will be
# shown in more detail. This can generate a lot of graphics.
GRAPHICAL_DEBUGGING = False

# The beta and datum avoidance zones overlap considerably. They clutter up
# the plots, and and it may be a waste of time including both of them.
# These flags can be used to see the result of turning the zones on and
# off. The flags can also be used to turn plots of the ellipses and padding
# on and off.
INCLUDE_BETA_ZONE = True
INCLUDE_DATUM_FOOT = True
INCLUDE_FIBRE_SNAG = True
PLOT_ELLIPSES = False
PLOT_PADDING = False

# Import standard Python libraries
import math, time

# Import common utility and plotting libraries.
import lib.util as util

try:
    import mocpath.plotting as plotting
except ImportError:
    plotting = None
    logger.warning("Could not import plotting utility. Plotting disabled.")

# Import constants from the utility module.
EPS = util.EPS            # Smallest possible floating point increment
MAXFLOAT = util.MAXFLOAT  # Largest possible floating point value.

# Mathematical constants
ROOT3 = util.ROOT3
ROOT3BY2 = util.ROOT3BY2
PI2 = util.PI2
PIBY2 = util.PIBY2

# Parity flags
PARITY_RIGHT = util.PARITY_RIGHT
PARITY_LEFT = util.PARITY_LEFT


# Global variables that describe the MOONS focal plane as a whole.
#
# Focal plane radius of curvature (mm)
FOCAL_PLANE_CURVATURE = 4212.0
# Focal plane diameter (mm)
FOCAL_PLANE_DIAMETER = 864.0
# Approximate focal plane plate scale (arcsec/mm)
FOCAL_PLANE_PLATE_SCALE = 1.718

# Global variables that define the default lengths and widths
# of the alpha and beta arms in millimetres.
#
# NOTE: In the final software, these parameters will be imported
# from a configuration file, not defined as constants.
#
# TODO: READ FROM CONFIGURATION FILE
INS_POS_LENGTH1 = 8.0     # Length of alpha arm (mm).
INS_POS_WIDTH1 = 8.0      # Width of alpha arm (mm).
ALPHA_LIMITS = [-180.0, 180.0]  # Min and max travel limits of alpha arm (deg)
ALPHA_DEFAULT = 0.0       # Default (starting) location for alpha arm (deg).
ALPHA_RPM_LIMITS = [0.73, 2.78] # Min and max motor speeds for alpha arm

INS_POS_LENGTH2 = 17.0    # Length of beta arm (mm).
INS_POS_WIDTH2 = 4.0      # Width of beta arm (mm).
BETA_LIMITS  = [-180.0, 150.0]  # Min and max travel limits of beta arm (deg)
BETA_DEFAULT = 0.0        # Default (starting) location for beta arm (deg).
BETA_RPM_LIMITS = [0.98, 3.75] # Min and max motor speeds for beta arm

MAX_WAVEFORM_STEPS = 128  # Maximum number of time steps in a motor waveform

# Avoidance zone parameters
INS_POS_B1  = 5.5  # Length of new metrology zone (mm)
INS_POS_B2  = 1.25 # Rear protruding length of metrology zone (mm)
INS_POS_TB2 = 1.5  # Length of small triangular avoidance zone (mm)
INS_POS_B3  = 6.0  # Length of triangular fibre snag zone zone.
INS_POS_B4  = 4.0  # Length of "safe" part of beta arm.
INS_POS_W3 = 2.75  # Width of triangular fibre snag zone zone.
INS_POS_W4 = INS_POS_WIDTH2 + 1.0 # Width of beta metrology (+beta edge) zone
INS_POS_D1  = 8.7   # Datum actuator length
INS_POS_DW  = INS_POS_WIDTH2      # Datum actuator width
INS_POS_D2  = 2.0   # Datum actuator backwards protrusion.
INS_POS_SAFETY = 0.5          # Safety tolerance added to all avoidance zones.

# Fibre holder parameters
INS_POS_MINDIST = 4.0    # Closest fibre approach distance (mm).
INS_POS_TOLER   = 20.0   # Fibre positioning tolerance (microns)
FIBRE_RADIUS    = INS_POS_MINDIST/2.0  # Radius of fibre holder in mm

# Default acquisition camera parameters
AC_FOV_MM = 23.0         # Field of view of acquisition camera in mm on focal plane
AC_FOV_RADIUS = AC_FOV_MM / 2.0
AC_FOV_ARCSEC = AC_FOV_MM * FOCAL_PLANE_PLATE_SCALE # Field of view in arcsec on sky

# Default fiducial parameters
FID_SEMI_MINOR = 12.0     # Semi-minor axis of avoidance zone
FID_SEMI_MAJOR = 14.0    # Semi-major axis of avoidance zone

# Conflict types
CONFLICT_OK = 0
CONFLICT_UNREACHABLE = 1
CONFLICT_LOCKED = 2
CONFLICT_TOOCLOSE = 3
CONFLICT_METROLOGY = 4
CONFLICT_BETA = 5
CONFLICT_DATUM = 6
CONFLICT_FIBRE_SNAG = 7
CONFLICT_AC = 8
CONFLICT_FID = 9

# ***
# ORIGINAL (EXTREMELY INEFFICIENT) CONFLICT FUNCTION
# SUPERCEDED BY THE NEW PositionerGrid.test_pair function.
def check_for_conflict(rcenpos1, thcenpos1, orient1,
                       rtarget1, thtarget1, parity1,
                       rcenpos2, thcenpos2, orient2,
                       rtarget2, thtarget2, parity2):
    """

    A function which implements the conflict check with a non-object-oriented
    interface. This version of the API would be used if the shared library is
    implemented in C.

    :Parameters:

    rcenpos1: float
        The radial distance of the centre of the first FPU,
        in polar coordinates, with respect to the centre of
        the focal plane (in mm).
    thcenpos1: float
        The theta angle of the centre of the first FPU,
        in polar coordinates, with respect to the centre of
        the focal plane (in radians)
    orient1: float.
        The rotation angle between the focal plane coordinate system
        and the local coordinate system of the first positioner
        (in radians). 0.0 means no local rotation.
    rtarget1: float
        The radial distance of the target assigned to the first FPU,
        in polar coordinates, with respect to the centre of
        the focal plane (in mm).
    thtarget1: float
        The theta angle of the target assigned to the first FPU,
        in polar coordinates,  with respect to the centre of
        the focal plane (in radians).
    parity1: int
        The arm elbow orientation to be adopted by the target
        assigned to the first FPU:

        * 1 means elbow right armed
        * -1 means elbow left armed

     rcenpos2: float
        The radial distance of the centre of the second FPU,
        in polar coordinates, with respect to the centre of
        the focal plane (in mm).
    thcenpos2: float
        The theta angle of the centre of the second FPU,
        in polar coordinates, with respect to the centre of
        the focal plane (in radians)
    orient2: float.
        The rotation angle between the focal plane coordinate system
        and the local coordinate system of the second positioner
        (in radians). 0.0 means no local rotation.
    rtarget2: float
        The radial distance of the target assigned to the second FPU,
        in polar coordinates, with respect to the centre of
        the focal plane (in mm).
    thtarget2: float
        The theta angle of the target assigned to the second FPU,
        in polar coordinates,  with respect to the centre of
        the focal plane (in radians).
    parity2: int
        The arm elbow orientation to be adopted by the target
        assigned to the second FPU:

        * 1 means elbow right armed
        * -1 means elbow left armed

    :Returns:

    in_conflict: boolean
        True if the two positioners are in conflict, or if any of the
        targets cannot be reached.
        False if both targets are ok and the positioners are not in conflict.

    """
    # NOTE: very inefficient because it creates, configures an destroys
    # the Positioner objects each time.

    # Create two positioner objects
    positioner1 = Positioner(1, r_centre_focal=rcenpos1,
                             theta_centre_focal=thcenpos1,
                             orient=orient1, column=0, row=0)
    positioner2 = Positioner(2, r_centre_focal=rcenpos2,
                             theta_centre_focal=thcenpos2,
                             orient=orient2, column=0, row=1)

    # Assign the targets
    result1 = positioner1.set_target( rtarget1, thtarget1, parity1 )
    result2 = positioner2.set_target( rtarget2, thtarget2, parity2 )

    if result1 and result2:
        # Check for conflict
        in_conflict = positioner1.in_conflict_with( positioner2 )
        return in_conflict
    else:
        # Targets cannot be reached.
        return True


# ADDITIONAL CLASS FROM PATH ANALYSIS RELEASE - 29 JUL 2016
class Motor(object):
    """

    Class representing the control parameters and path for a MOONS fibre
    positioner motor. Used by the path analysis software.

    """
    def __init__(self, sim_length, position_target, max_speed, min_speed):
        """

        Constructor for Motor class.

        :Parameters:

        sim_length: int
            Number of elements in motor path array.
        position_target: float
            Target motor angle (in radians)
        max_speed: float
            The maximum allowed motor speed (in radians per second)
        min_speed: float
            The minimum allowed motor speed (in radians per second)

        :Returns:

        New Motor object

        """
        self.position_array = [0.0] * sim_length
        self.speed_array = [0.0] * (sim_length - 1)
        self.position_array[sim_length - 1] = position_target
        self.maximum_speed = max_speed
        self.minimum_speed = min_speed


class FiducialMarker(object):
    """

    Class representing a MOONS Fidicual Marker.

    Each Fidicual Marker is identified by a unique identifier and a location
    in MOONS focal plane coordinates.

    NOTE: Fiducial markers are located in (X,Y) Cartesian coordinates.
    They can occupy a cell within the fibre positioner grid or they can
    also be located outside the grid around te edge of the focal plane.

    """

    def __init__(self, ident, x_centre_focal, y_centre_focal,
                 orient, avoid_minor=9.0, avoid_major=14.0,
                 column=None, row=None ):
        """

        Constructor for FiducialMarker class.

        :Parameters:

        ident: int
            Unique identifier for this fidicual marker.
            IDs must be unique amongst fiducial markers. A fiducial marker
            can have the same ID as an acquisition camera or fibre positioner.
        x_centre_focal: float
            The X Cartesian coordinate of the fidicual marker centre
            with respect to the centre of the focal plane (in mm).
        y_centre_focal: float
            The Y Cartesian coordinate of the fidicual marker centre
            with respect to the centre of the focal plane (in mm).
        orient: float
            The rotation angle between the focal plane coordinate system
            and the fiducial marker local coordinate system (in radians).
            Defaults to 0.0 (no local rotation)
        column: int (optional)
            Nearest column number on the positioner grid, if applicable.
            A fiducial marker must not occupy the same grid cell as a
            fibre positioner, although it can share a grid cell with an
            acquisition camera.
        row: int (optional)
            Nearest column number on the positioner grid, if applicable.
            A fiducial marker must not occupy the same grid cell as a
            fibre positioner, although it can share a grid cell with an
            acquisition camera.
        avoid_minor: float (optional)
            Semi-minor axis of the field of view to be kept clear (mm)
            The minor axis is oriented perpendicular to the radial
            direction in the focal plane.
        major_major: float (optional)
            Semi-major axis of the field of view to be kept clear (mm)
            The major axis is oriented along the radial
            direction in the focal plane.

        :Returns:

        New FiducialMarker object

        """
        # Define a unique ID and unique name name for this fiducial marker.
        self.ident = int(ident)
        self.name = "FID%d" % self.ident
        logger.debug("Creating fidicual marker %s" % self.name)

        # Define the location of the centre and the orientation of the
        # fiducial marker on the MOONS focal plane, in polar and Cartesian
        # coordinates.
        # TODO: Is it worth moving these common operations into a GridObject class?
        self.x_centre_focal = float(x_centre_focal)
        self.y_centre_focal = float(y_centre_focal)
        (self.r_centre_focal, self.theta_centre_focal) = \
            util.cartesian_to_polar(self.x_centre_focal, self.y_centre_focal)
        self.orient = float(orient)

        # If the object is part of a grid, record its  location
        # in hexagonal grid coordinates (if defined).
        if column is not None and row is not None:
            self.column = int(column)
            self.row = int(row)
        else:
            self.column = None
            self.row = None

        # Define the avoidance zone. The ellipse is oriented with major axis
        # radially away from the centre (since it is viewed by metrology
        # cameras spaced around the edge of the field).
        self.avoid_minor = avoid_minor
        self.avoid_major = avoid_major
        avoid_tilt = math.radians(90.0) - self.theta_centre_focal
        self.avoidance_ellipse = (self.x_centre_focal, self.y_centre_focal,
                                  self.avoid_major, self.avoid_minor,
                                  avoid_tilt)

        # Initialise the fidicual marker object
        self.initialise()

    def __str__(self):
        """

        Return a readable string describing the fidicual marker.

        :Returns:

        strg: str
            Readable string

        """
        strg = "Fidicual Marker object \'%s\' " % self.name
        if self.column is not None and self.row is not None:
            strg += "[%d,%d] " % (self.column, self.row)
        theta_degrees = math.degrees(self.theta_centre_focal)
        strg += "located on the focal plane at (r=%.3f, theta=%.3f deg) " % \
            (self.r_centre_focal, theta_degrees)
        strg += "or (x=%.3f, y=%.3f).\n" % \
            (self.x_centre_focal, self.y_centre_focal)
        if abs(self.orient) > 0.0:
            orient_degrees = math.degrees(self.orient)
            strg += "  Marker oriented by a0=%.3f (deg)." % orient_degrees
        strg += "  Avoidance ellipse: minor axis=%.3f," % self.avoid_minor
        strg += "  major axis=%.3f." % self.avoid_major
        if self.in_conflict:
            strg += "\n  Fidicual Marker IN CONFLICT: " + self.conflict_reason

        return strg

    # ***
    def initialise(self):
        """

        Initialise the fidicual marker.

        :Parameters:

        None

        :Attributes:

        The function sets the following attributes.

        active, in_conflict, conflict_type, conflict_reason
            Parameters cleared.

        :Returns:

        None

        """
        # A fiducial marker starts inactive.
        self.active = False
        # Start with no conflict and a blank conflict message string
        self.in_conflict = False
        self.conflict_type = CONFLICT_OK
        self.conflict_reason = ''

    def activate(self):
        """

        Activate the fidicual marker. This will enable its avoidance zones.

        :Parameters:

        None

        :Attributes:

        The function sets the following attributes.

        active
            Set to True.

        :Returns:

        None

        """
        self.active = True

    def plot(self, description='', plotlabel=True, plotfig=None, showplot=True):
        """

        Plot the configuration of the fidicual marker.

        NOTE: This function requires access to the plotting
        module, plotting.py. If the module can't be imported
        the function returns with an apologetic message.

        :Parameters:

        description: str, optional
            Optional description to be added to the positioner plot.
        plotlabel: bool, optional
            If True, label the plot with the fidicual marker name.
            Defaults to True,
        plotfig: matplotlib Figure object, optional
            An existing matplotlib figure onto which to add the plot.
            If not given, a new figure is created.
        showplot: boolean, optional
            If True (the default) show the plot when finished.
            Set of False if you are overlaying several plots, and
            this plot is not the last.

        :Returns:

        plotfig: matplotlib Figure object
            A matplotlib figure containing the plot.

        """
        if plotting is None:
            logger.warning("Plotting is disabled")
            return

        # Create a new plotting figure, if necessary.
        if plotfig is None:
            plotfig = plotting.new_figure(1, figsize=(10,9), stitle='')

        # Plot a small diamond marker at the centre of the fiducial and a
        # dotted line showing the avoidance ellipse (green if ok and red
        # if in conflict)
        if self.in_conflict:
            fmt = 'rd'
            colour = 'r '
        else:
            fmt = 'gd'
            colour = 'g '
        if self.active:
            style = '--'
        else:
            style = ':'
        (xcen, ycen, major, minor, tilt) = self.avoidance_ellipse
        plotaxis = plotting.plot_xy( [xcen], [ycen],
                                    plotfig=plotfig, linefmt=fmt, linestyle=' ',
                                    showplot=False )
        plotaxis = plotting.plot_ellipses( [xcen], [ycen], major, minor,
                        tilt, plotfig=plotfig, plotaxis=plotaxis,
                        facecolor='w', ellipsecolor=colour, linefmt=colour,
                        linestyle=style, linewidth=2, grid=False,
                        showplot=False )

        # Add a label giving the fiducial marker ID
        if plotlabel:
            xtxt = self.x_centre_focal - self.avoid_minor/3.0
            ytxt = self.y_centre_focal - self.avoid_minor/3.0
            plotaxis = plotting.plot_text( str(self.name), xtxt, ytxt,
                                       plotfig=plotfig, plotaxis=plotaxis,
                                       showplot=False)

        if showplot:
            plotting.show_plot()
            plotting.close()
        return plotfig

    def plot_status(self, description='', plotfig=None, showplot=True):
        # A status plot shows the same ellipses as a normal plot.
        self.plot(description=description, plotfig=plotfig, plotlabel=False,
                  showplot=showplot)


class AcquisitionCamera(object):
    """

    Class representing a MOONS Acquisition Camera.

    Each acquisition camera is identified by a unique identifier and a location
    in MOONS focal plane coordinates.

    NOTE: Acquisition cameras are located in (R,THETA) polar coordinates.
    They must occupy a cell within the fibre positioner grid. The same cell
    cannot be shared with a fibre positioner, but an acquisition camera and
    fiducial marker can share the same cell.

    """

    def __init__(self, ident, r_centre_focal, theta_centre_focal,
                 orient, fov_radius, column=None, row=None):
        """

        Constructor for AcquisitionCamera class.

        :Parameters:

        ident: int
            Unique identifier for this acquisition camera, AC
            IDs must be unique amongst acquisition cameras. An acquisition
            camera can have the same ID as a fiducial marker or fibre
            positioner.
        r_centre_focal: float
            The radial distance of the AC centre, in polar coordinates,
            with respect to the centre of the focal plane (in mm).
        theta_centre_focal: float
            The theta angle of the AC centre, in polar coordinates,
            with respect to the centre of the focal plane (in radians).
        orient: float
            The rotation angle between the focal plane coordinate system
            and the AC local coordinate system (in radians).
            Defaults to 0.0 (no local rotation)
        fov_radius: float
            The radius of the field of view to be kept clear (mm).
        column: int (optional)
            Column number of this acquisition camera (if part of a grid).
            An acquisition camera cannot occupy the same grid cell as
            a fibre positioner, although it can share a grid cell with
            a fiducial marker.
        row: int (optional)
            Row number of this acquisition camera (if part of a grid)
            An acquisition camera cannot occupy the same grid cell as
            a fibre positioner, although it can share a grid cell with
            a fiducial marker.

        :Returns:

        New AcquisitionCamera object

        """
        # Define a unique ID and unique name for this acquisition camera.
        self.ident = int(ident)
        self.name = "AC%d" % self.ident
        logger.debug("Creating acquisition camera %s" % self.name)

        # Define the location of the centre and the orientation of the
        # acquisition camera on the MOONS focal plane, in polar and
        # Cartesian coordinates.
        # TODO: Is it worth moving these common operations into a GridObject class?
        self.r_centre_focal = float(r_centre_focal)
        self.theta_centre_focal = float(theta_centre_focal)
        (self.x_centre_focal, self.y_centre_focal) = \
            util.polar_to_cartesian(self.r_centre_focal, self.theta_centre_focal)
        self.orient = float(orient)

        # If the object is part of a grid, record its  location
        # in hexagonal grid coordinates (if defined).
        if column is not None and row is not None:
            self.column = int(column)
            self.row = int(row)
        else:
            self.column = None
            self.row = None

        # Define the default avoidance zone
        self.fov_radius_default = fov_radius

        # Initialise the acquisition camera object
        self.initialise()

    def __str__(self):
        """

        Return a readable string describing the acquisition camera.

        :Returns:

        strg: str
            Readable string

        """
        strg = "Acquisition Camera object \'%s\' " % self.name
        if self.column is not None and self.row is not None:
            strg += "[%d,%d] " % (self.column, self.row)
        theta_degrees = math.degrees(self.theta_centre_focal)
        strg += "located on the focal plane at (r=%.3f, theta=%.3f deg) " % \
            (self.r_centre_focal, theta_degrees)
        strg += "or (x=%.3f, y=%.3f).\n" % \
            (self.x_centre_focal, self.y_centre_focal)
        if abs(self.orient) > 0.0:
            orient_degrees = math.degrees(self.orient)
            strg += "  Camera oriented by a0=%.3f (deg)." % orient_degrees
        strg += "  Default field of view radius=%.3f." % self.fov_radius_default
        if self.active:
            theta_degrees = math.degrees(self.theta_star_focal)
            strg += "\n  Guide star defined at (r=%.3f, theta=%.3f deg) " % \
            (self.r_star_focal, theta_degrees)
            strg += "or (x=%.3f, y=%.3f) " % \
                (self.x_star_focal, self.y_star_focal)
            strg += "with avoidance radius=%.3f." % self.fov_radius
        if self.in_conflict:
            strg += "\n  Acquisition Camera IN CONFLICT: " + self.conflict_reason

        return strg

    # ***
    def initialise(self):
        """

        Initialise the acquisition camera.

        :Parameters:

        None

        :Attributes:

        The function sets the following attributes.

        active, x_star_focal, y_star_focal, fov_radius
            Set to default.
        in_conflict, conflict_type, conflict_reason
            Parameters cleared.

        :Returns:

        None

        """
        # An acquisition camera starts inactive with no guide star defined.
        self.active = False
        self.r_star_focal = self.r_centre_focal
        self.theta_star_focal = self.theta_centre_focal
        self.x_star_focal = self.x_centre_focal
        self.y_star_focal = self.y_centre_focal
        self.fov_radius = self.fov_radius_default

        # Start with no conflict and a blank conflict message string
        self.in_conflict = False
        self.conflict_type = CONFLICT_OK
        self.conflict_reason = ''

    # ***
    def set_star(self, r_star_focal=None, theta_star_focal=None,
                 fov_radius=None):
        """

        Define a guide star and activate the acquisition camera.

        :Parameters:

        r_star_focal: float, optional
            R coordinate of guide star on focal plane (in mm)
            Defaults to the centre of the acquisition camera
        theta_star_focal: float, optional
            Theta coordinate of guide star on focal plane (in degrees)
            Defaults to the centre of the acquisition camera.
        fov_radius: float, optional
            Avoidance zone surrounding guide star (in arcsec).
            Cannot be larger than the default field of view of
            the acquisition camera.
            Dafaults to the field of view of the acquisition camera

        :Attributes:

        The function sets the following attributes.



        in_conflict, conflict_type, conflict_reason
            Parameters cleared.

        :Returns:

        None

        """
        if r_star_focal is not None and theta_star_focal is not None:
            theta_radians = math.radians(theta_star_focal)
            self.r_star_focal = r_star_focal
            self.theta_star_focal =  theta_radians
            (x_star_focal, y_star_focal) = \
                util.polar_to_cartesian(r_star_focal, theta_radians)
            self.x_star_focal = x_star_focal
            self.y_star_focal = y_star_focal
            xdiff = self.x_star_focal - self.x_centre_focal
            ydiff = self.y_star_focal - self.y_centre_focal
            rdiff = math.sqrt( xdiff * xdiff + ydiff * ydiff )
            fov_maximum = self.fov_radius_default - rdiff + 1.0
        else:
            self.x_star_focal = self.x_centre_focal
            self.y_star_focal = self.y_centre_focal
            fov_maximum = self.fov_radius_default

        if fov_radius is not None:
            fov_mm = fov_radius / FOCAL_PLANE_PLATE_SCALE
            if fov_mm < fov_maximum:
                self.fov_radius = fov_mm
            else:
                self.fov_radius = fov_maximum
        else:
            self.fov_radius = self.fov_radius_default

        # The acquisition camera is now active
        self.active = True

        # Start with no conflict and a blank conflict message string
        self.in_conflict = False
        self.conflict_type = CONFLICT_OK
        self.conflict_reason = ''

    def plot(self, description='', plotlabel=True, plotfig=None, showplot=True):
        """

        Plot the configuration of the acquisition camera.

        NOTE: This function requires access to the plotting
        module, plotting.py. If the module can't be imported
        the function returns with an apologetic message.

        :Parameters:

        description: str, optional
            Optional description to be added to the positioner plot.
        plotlabel: bool, optional
            If True, label the plot with the acquisition camera name.
            Defaults to True,
        plotfig: matplotlib Figure object, optional
            An existing matplotlib figure onto which to add the plot.
            If not given, a new figure is created.
        showplot: boolean, optional
            If True (the default) show the plot when finished.
            Set of False if you are overlaying several plots, and
            this plot is not the last.

        :Returns:

        plotfig: matplotlib Figure object
            A matplotlib figure containing the plot.

        """
        if plotting is None:
            logger.warning("Plotting is disabled")
            return

        # Create a new plotting figure, if necessary.
        if plotfig is None:
            plotfig = plotting.new_figure(1, figsize=(10,9), stitle='')

        # Plot a solid circle showing the maximum field of view
        # (blue if ok and red if in conflict)
        xcentres = [self.x_centre_focal]
        ycentres = [self.y_centre_focal]
        if self.in_conflict:
            colour = 'r '
        else:
            colour = 'b '
        plotaxis = plotting.plot_circles( xcentres, ycentres,
                        self.fov_radius_default, plotfig=plotfig,
                        facecolor='w', circcolor=colour, linefmt=colour,
                        linestyle='-', linewidth=3, grid=False,
                        showplot=False )

        # If active, plot a * and a dotted circle showing the actual field of view
        # (blue if ok and red if in conflict)
        if self.active:
            xcentres = [self.x_star_focal]
            ycentres = [self.y_star_focal]
            plotaxis = plotting.plot_xy( xcentres, ycentres,
                                        plotfig=plotfig, linefmt='k*', linestyle=' ',
                                        showplot=False )
            plotaxis = plotting.plot_circles( xcentres, ycentres,
                            self.fov_radius, plotfig=plotfig,
                            facecolor='w', circcolor=colour, linefmt=colour,
                            linestyle='--', linewidth=2, grid=False,
                            showplot=False )

        # Add a label giving the acquisition camera ID
        if plotlabel:
            xtxt = self.x_centre_focal - self.fov_radius/3.0
            ytxt = self.y_centre_focal - self.fov_radius/3.0
            plotaxis = plotting.plot_text( str(self.name), xtxt, ytxt,
                                       plotfig=plotfig, plotaxis=plotaxis,
                                       showplot=False)
        if showplot:
            plotting.show_plot()
            plotting.close()
        return plotfig

    def plot_status(self, description='', plotfig=None, showplot=True):
        """

        Plot the status of the acquisition camera.

        NOTE: This function requires access to the plotting
        module, plotting.py. If the module can't be imported
        the function returns with an apologetic message.

        :Parameters:

        description: str, optional
            Optional description to be added to the positioner plot.
        plotfig: matplotlib Figure object, optional
            An existing matplotlib figure onto which to add the plot.
            If not given, a new figure is created.
        showplot: boolean, optional
            If True (the default) show the plot when finished.
            Set of False if you are overlaying several plots, and
            this plot is not the last.

        :Returns:

        plotfig: matplotlib Figure object
            A matplotlib figure containing the plot.

        """
        if plotting is None:
            logger.warning("Plotting is disabled")
            return

        # Create a new plotting figure, if necessary.
        if plotfig is None:
            plotfig = plotting.new_figure(1, figsize=(10,9), stitle='')

        # Plot acquisition cameras as blue filled black circles..
        facecolor = 'blue'
        circcolor = 'b '
        linefmt = 'k '

        # Plot acquisition cameras as filled circles.
        # Add a '*' if the camera has a guide star.
        xcentres = [self.x_centre_focal]
        ycentres = [self.y_centre_focal]
        plotaxis = plotting.plot_circles( xcentres, ycentres,
                        self.fov_radius_default, plotfig=plotfig,
                        facecolor=facecolor, circcolor=circcolor, linefmt=linefmt,
                        linestyle='-', linewidth=3, grid=False, filled=True,
                        showplot=False )
        if self.active:
            plotaxis = plotting.plot_xy( xcentres, ycentres,
                                        plotfig=plotfig, plotaxis=plotaxis,
                                        linefmt='k*', linestyle=' ',
                                        showplot=False )

        if showplot:
            plotting.show_plot()
            plotting.close()
        return plotfig


class Positioner(object):
    """

    Class representing a MOONS fibre positioner unit (FPU).

    Each positioner is identified by a unique identifier and a location
    in MOONS focal plane coordinates. Coordinates within a MOONS fibre
    positioner can be described in 4 different ways:

    1) Focal plane polar coordinates: r,theta relative to the centre
       of the focal plane.

    2) Focal plane Cartesian coordinates: x,y relative to the centre
       of the focal plane.

    3) Local polar coordinates: r,theta relative to the centre of
       the positioner.

    4) Local Cartesian coordinates: x,y relative to the centre of
       the positioner.

    External interfaces use focal plane coordinates (1) or (2) because
    those are universal over the whole focal plane, and can detect
    conflicts between positioners.

    Internal interfaces use local coordinates (3) or (4), which are
    valid within the patrol zone of one positioner only.

    A fibre positioner cannot occupy the same grid cell as an acquisition
    camera or fiducial marker.

    """

    # To start with, common parameter are not defined.
    _params_defined = False

    # ***
    @classmethod
    def define_common_params(cls, length1, length2):
        """

        Class method which defines the parameters which are common to all
        fibre positioners. These parameters only need to be calculated once.

        :Parameters:

        """
        # Define the parameters valid for all positioner arms (saving the
        # square of those parameters for optimisation purposes).
        cls.length1 = float(length1)
        cls.length2 = float(length2)
        cls._length1sq = cls.length1 * cls.length1
        cls._length2sq = cls.length2 * cls.length2
        # Inner and outer limits of the positioner patrol zone.
        cls.outer = cls.length1 + cls.length2
        cls.inner = abs(cls.length2 - cls.length1)
        cls._outersq = cls.outer * cls.outer
        cls._innersq = cls.inner * cls.inner
        # A safe distance beyond which positioners cannot conflict
        cls.safedist = cls.outer + INS_POS_MINDIST
        cls._safedistsq = cls.safedist * cls.safedist
        logger.debug("Target distance beyond which positioners cannot " + \
                     "conflict is %f" % cls.safedist)
        # A limiting distance inside which far neighbours cannot conflict
        cls.limitdist = (cls.outer - INS_POS_MINDIST) * (ROOT3 - 1.0)
#         cls._limitdistsq = cls.limitdist * cls.limitdist
        logger.debug("Reach distance inside which far neighbours cannot " + \
                     "conflict is %f" % cls.limitdist)

        # The first (length2-length1)/2.0 of the beta arm is inside the safe
        # zone. Since the boundary of the safe zone is a circular arc, the
        # following calculation takes into account the finite width of the
        # beta arm and calculates the length at the edges of the arm.
        cls._bsafesq = (((cls.length2-cls.length1)/2.0) * \
                            ((cls.length2-cls.length1)/2.0)) - \
                        (INS_POS_WIDTH2/2.0 * INS_POS_WIDTH2/2.0)
        cls.bsafe = math.sqrt(cls._bsafesq)
        # Allow this parameter to be made stricter with the INS_POS_B4 parameter
        if cls.bsafe > INS_POS_B4:
            cls.bsafe = INS_POS_B4
        logger.debug("Portion of beta arm safe from conflict is %f" % cls.bsafe)

        # Determine the critical extension length beyond which a
        # limited range of travel for the beta arm makes a difference
        # to the choice of parity.
        if BETA_LIMITS and BETA_LIMITS[1] < 180.0:
            upperlimit = math.radians(BETA_LIMITS[1])
            lcritsq = (cls.length1 * cls.length1) + \
                      (cls.length2 * cls.length2) - \
                      (2.0 * cls.length1 * cls.length2 * math.cos(upperlimit))
            cls.lengthcrit = math.sqrt(lcritsq)
            cls.paritycrit = PARITY_RIGHT
        elif BETA_LIMITS and BETA_LIMITS[0] > -180.0:
            upperlimit = math.radians(BETA_LIMITS[0])
            lcritsq = (cls.length1 * cls.length1) + \
                      (cls.length2 * cls.length2) - \
                      (2.0 * cls.length1 * cls.length2 * math.cos(upperlimit))
            cls.lengthcrit = math.sqrt(lcritsq)
            cls.paritycrit = PARITY_LEFT
        else:
            cls.lengthcrit = cls.length2
            cls.paritycrit = PARITY_RIGHT
        logger.debug("Only %s parity is possible beyond %f" % \
                    (util.elbow_parity_str(cls.paritycrit), cls.lengthcrit))
        # The common parameters are now defined.
        cls._params_defined = True

    def __init__(self, ident, r_centre_focal, theta_centre_focal, orient,
                 column=None, row=None, simulated=True, locked=False,
                 rfibre=None, thfibre=None, pfibre=PARITY_RIGHT,
                 length1=INS_POS_LENGTH1, length2=INS_POS_LENGTH2):
        """

        Constructor for Positioner class.

        :Parameters:

        ident: int
            Unique identifier for this fibre positioner
            IDs must be unique amongst fibre positioners.
        r_centre_focal: float
            The radial distance of the FPU centre, in polar coordinates,
            with respect to the centre of the focal plane (in mm)
        theta_centre_focal: float
            The theta angle of the FPU centre, in polar coordinates,
            with respect to the centre of the focal plane (in radians)
        orient: float
            The rotation angle between the focal plane coordinate system
            and the FPU local coordinate system (in radians).
            Defaults to 0.0 (no local rotation)
        column: int (optional)
            Column number of this positioner (if part of a grid)
        row: int (optional)
            Row number of this positioner (if part of a grid)
        simulated: bool (optional)
            True if the positioner is simulated. The default is True.
            (since this code is a simulation).
        locked: bool (optional)
            True if the positioner is locked (broken). The default is False.
        rfibre: float (optional)
            Starting location for the fibre.
            Defaults to the initialised position, near datum.
        thfibre: float (optional)
            Starting location for the fibre.
            Defaults to the initialised position, near datum.
        pfibre: int (optional)
            Starting parity for the fibre.
            Defaults to right-armed parity.
        length1, float (optional)
            Length of alpha arm in mm.
            Defaults to the designed length, INS_POS_LENGTH1.
        length2, float (optional)
            Length of beta arm in mm.
            Defaults to the designed length, INS_POS_LENGTH2.

        :Returns:

        New Positioner object

        """
        # Define a unique ID and unique name for this positioner.
        self.ident = int(ident)
        self.name = "POS%d" % self.ident
        logger.debug("Creating positioner %s" % self.name)

        # Define the location of the centre and the orientation of the
        # fibre positioner on the MOONS focal plane, in polar and Cartesian
        # coordinates.
        # TODO: Is it worth moving these common operations into a GridObject class?
        self.r_centre_focal = float(r_centre_focal)
        self.theta_centre_focal = float(theta_centre_focal)
        (self.x_centre_focal, self.y_centre_focal) = \
            util.polar_to_cartesian(self.r_centre_focal, self.theta_centre_focal)
        self.orient = float(orient)

        # Record the location in the focal plane hexagonal grid (if provided).
        if column is not None and row is not None:
            self.column = int(column)
            self.row = int(row)
        else:
            self.column = None
            self.row = None

        # Initialise the list of neighbouring positioners
        self.near_neighbours = []
        self.far_neighbours = []

        # Initialise the list of nearby acquisition cameras and fiducial markers
        self.near_cameras = []
        self.near_fiducials = []

        # Is the positioner simulated or broken?
        self.simulated = bool(simulated)
        self.locked = bool(locked)

        # If not already defined, define the parameters common to all
        # fibre positioners.
        if not Positioner._params_defined:
            Positioner.define_common_params(length1, length2)

        # Start with no target assigned.
        # Start at the initialised location, near datum, unless a specific
        # starting location has been defined.
        self.starting_rfibre = rfibre
        self.starting_thfibre = thfibre
        self.starting_parity = pfibre
        self.initialise()

        # There are no associated motors until the path analysis software
        # defines them.
        self.motor1 = None
        self.motor2 = None

    def __str__(self):
        """

        Return a readable string describing the positioner.

        :Returns:

        strg: str
            Readable string

        """
        strg = "Positioner object \'%s\' " % self.name
        if self.column is not None and self.row is not None:
            strg += "[%d,%d] " % (self.column, self.row)
        theta_degrees = math.degrees(self.theta_centre_focal)
        strg += "located on the focal plane at (r=%.3f, theta=%.3f deg) " % \
            (self.r_centre_focal, theta_degrees)
        strg += "or (x=%.3f, y=%.3f).\n" % \
            (self.x_centre_focal, self.y_centre_focal)
        if abs(self.orient) > 0.0:
            orient_degrees = math.degrees(self.orient)
            strg += "  Positioner oriented by a0=%.3f (deg)." % orient_degrees
        strg += "  There are %d near and %d far neighbouring positioners" % \
            (len(self.near_neighbours), len(self.far_neighbours))
        nacs = len(self.near_cameras)
        if nacs > 0:
            strg += " + %d close acquisition camera" % nacs
            if nacs > 1:
                strg += "s"
        nfids = len(self.near_fiducials)
        if nfids > 0:
            strg += " + %d close fiducial marker" % nfids
            if nfids > 1:
                strg += "s"
        strg += ".\n"

        if self.locked:
            strg += "  Positioner is LOCKED. No target can be assigned."
        else:
            if self.target_assigned:
                strg += "  Target assigned at "
                strg += "(r=%.3f, theta=%.3f deg, parity=%s).\n  " % \
                    (self.r_target_focal, math.degrees(self.theta_target_focal),
                     util.elbow_parity_str(self.target_parity))
            else:
                strg += "  Target NOT assigned. "
        theta_degrees = math.degrees(self.theta_fibre_focal)
        strg += "Fibre located on the focal plane at (r=%.3f, theta=%.3f deg) " % \
            (self.r_fibre_focal, theta_degrees)
        strg += "or (x=%.3f, y=%.3f),\n" % \
            (self.x_fibre_focal, self.y_fibre_focal)
        theta_local = math.degrees(self.theta_fibre_local)
        strg += "    equivalent to local coordinates (r=%.3f, theta=%.3f deg)" % \
            (self.r_fibre_local, theta_local)
        strg += " or (x=%.3f, y=%.3f).\n" % \
            (self.x_fibre_local, self.y_fibre_local)

        strg += "  Elbow %s parity leads to elbow location" % \
            util.elbow_parity_str(self.parity)
        strg += " of (x=%.3f, y=%.3f) in local coordinates." % \
            (self.x_elbow_local, self.y_elbow_local)

        if self.in_conflict:
            strg += "\n  Positioner IN CONFLICT: " + self.conflict_reason
        return strg

    def add_neighbour(self, neighbour, near=True):
        """

        Add a link to a neighbouring positioner whose patrol zone
        overlaps with this one.

        :Parameters:

        neighbour: Positioner object
            Neighbouring positioner (ignored if None)
        near: bool (optional)
            Set True if the other object is a near neighbour
            or False if the other object is a far neighbour.
            Defaults to True.

        """
#         print("%s: Adding neighbour %s" % (self.name, neighbour.name))
        if neighbour is not None:
            assert isinstance(neighbour, Positioner)
            if near:
                self.near_neighbours.append(neighbour)
            else:
                self.far_neighbours.append(neighbour)

    def get_neighbours(self):
        """

        Return a list of neighbours of the current positioner.

        """
        # Append the lists of near and far neighbours.
        return self.near_neighbours + self.far_neighbours

    def add_acquisition_camera(self, camera):
        """

        Add a link to a fiducial marker which can interfer with the
        patrol zone of this positioner.

        :Parameters:

        camera: AcquisitionCamera object
            Neighbouring object (ignored if None)

        """
#         print("%s: Adding acquisition camera %s" % (self.name, camera.name))
        if camera is not None:
            assert isinstance(camera, AcquisitionCamera)
            self.near_cameras.append(camera)

    def add_fiducial(self, fiducial):
        """

        Add a link to a fiducial marker which can interfer with the
        patrol zone of this positioner.

        :Parameters:

        fiducial: FiducialMarker object
            Neighbouring object (ignored if None)

        """
#         print("%s: Adding fiducial %s" % (self.name, fiducial.name))
        if fiducial is not None:
            assert isinstance(fiducial, FiducialMarker)
            self.near_fiducials.append(fiducial)

    # ***
    def initialise(self):
        """

        Move the positioner to its initial location.

        :Parameters:

        None

        :Attributes:

        The function sets the following attributes.

        target_assigned
            Flag cleared. Any target is unassigned.
        r_target_focal, theta_target_focal
            Requested focal plane polar coordinates of the target.
            Initialised.
        target_parity
            Requested elbow parity of the target. Initialised.
        in_conflict, conflict_type, conflict_reason
            Parameters cleared.
        r_fibre_focal, theta_fibre_focal
            Focal plane polar coordinates of the default location.
        x_fibre_focal, y_fibre_focal
            Focal plane Cartesian coordinates of the default location.
        r_fibre_local, theta_fibre_focal
            Local polar coordinates of the default location.
        x_fibre_local, y_fibre_local
            Local Cartesian coordinates of the default location.
        parity
            Elbow parity of the default positioner orientation.
        x_elbow_local, y_elbow_local
            Local Cartesian coordinates of the beta axis (or elbow)
            at the default location.

        :Returns:

        None

        """
        # Initialise the positioner
        # Start with no conflict and a blank conflict message string
        self.in_conflict = False
        self.conflict_type = CONFLICT_OK
        self.conflict_reason = ''

        # At first there is no target assigned
        self.target_assigned = False
        self.r_target_focal = None
        self.theta_target_focal = None
        self.target_parity = PARITY_RIGHT
        # The positioner starts at its default starting position,
        # unless a different starting location has been specified.
        if self.starting_rfibre is not None and self.starting_thfibre is not None:
            # An alternative start location has been defined.
            result = self.set_target(self.starting_rfibre, self.starting_thfibre,
                            self.starting_parity, starting_location=True)
            if not result:
                strg = "Cannot move positioner %d to its initial location (%f,%f)" \
                    % (self.ident, self.starting_rfibre, self.starting_thfibre)
                raise ValueError(strg)
        else:
            # Move to default location
            self.set_arm_angles( ALPHA_DEFAULT, BETA_DEFAULT )
            strg = "Target initialised at local x=%f,y=%f --> r=%f,theta=%f at %s parity" % \
                (self.x_fibre_local, self.y_fibre_local,
                 self.r_fibre_local, math.degrees(self.theta_fibre_local),
                 util.elbow_parity_str(self.parity))
            if abs(self.orient) > 0.0:
                strg += " (orient=%f (deg))" % math.degrees(self.orient)
            logger.debug(strg)

    # ***
    def can_move_to(self, r_target_focal, theta_target_focal, parity):
        """

        Determine whether this positioner is capable of moving to
        the specified target.

        NOTE: No collisions are detected by this function. It merely tests
        whether a target is within the patrol zone of this positioner.

        :Parameters:

        r_target_focal: float
            The radial distance of the target, in polar coordinates,
            with respect to the centre of the focal plane (in mm)
        theta_target_focal: float
            The theta angle of the target, in polar coordinates,
            clockwise from Y axis,
            with respect to the centre of the focal plane (in radians)
        parity: int
            The arm elbow orientation to be adopted:

            * 1 means elbow right armed
            * -1 means elbow left armed

        :Attributes:

        The function sets the following attributes.

        in_conflict, conflict_type, conflict_reason
            Parameters which indicate if the positioner is in conflict.

        :Returns:

        target_reachable: boolean
            True if the target can be reached
            False if the target cannot be reached

        """
        logger.debug("Can positioner %s move to (r=%f, theta=%f) with %s parity?" %
                      (self.name, r_target_focal,
                       math.degrees(theta_target_focal),
                       util.elbow_parity_str(parity)))
        # Solve the triangle cosine rule to determine the local R
        # for the target
        rls = r_target_focal * r_target_focal + \
              self.r_centre_focal * self.r_centre_focal - \
              2.0 * r_target_focal * self.r_centre_focal * \
              math.cos(theta_target_focal-self.theta_centre_focal)
        r_target_local = math.sqrt(rls)

        # Can the positioner reach this target?
        # TODO: Should INS_POS_TOLER be taken into account?
        if r_target_local >= self.inner and r_target_local <= self.outer:
            # If the range of travel of the beta arm is limited, targets beyond
            # a certain critical distance can only be reached at one parity.
            if r_target_local > self.lengthcrit and parity != self.paritycrit:
                self.conflict_type = CONFLICT_UNREACHABLE
                self.conflict_reason = "%s: Wrong target parity at this reach " % \
                    self.name
                self.conflict_reason += "%.3f > %.3f" % (r_target_local, self.lengthcrit)
                logger.info(self.conflict_reason)
                return False
            else:
                logger.debug("YES.")
                return True
        else:
            self.conflict_type = CONFLICT_UNREACHABLE
            self.conflict_reason = "%s: Target outside patrol zone." % self.name
            logger.info(self.conflict_reason)
            return False

    # ***
    def set_target(self, r_target_focal, theta_target_focal, parity,
                   starting_location=False):
        """

        Propose a target for this positioner. If the target is reachable,
        the simulated positioner is moved immediately to the given
        coordinates and a new target is defined.

        :Parameters:

        r_target_focal: float
            The radial distance of the target, in polar coordinates,
            with respect to the centre of the focal plane (in mm)
        theta_target_focal: float
            The theta angle of the target, in polar coordinates,
            clockwise from Y axis,
            with respect to the centre of the focal plane (in radians)
        parity: int
            The arm elbow orientation to be adopted:

            * 1 means elbow right armed
            * -1 means elbow left armed

        starting_location: bool (optional)
            If True, defines starting location, not a new target.
            The default is False.

        :Attributes:

        The function sets the following attributes.

        target_assigned
            A Boolean flag which indicates if a target is successfully
            assigned to this positioner.
            Only defined when starting_location is False.
        r_target_focal, theta_target_focal
            Requested focal plane polar coordinates of the target.
        target_parity
            Requested elbow parity of the target.
        in_conflict, conflict_type, conflict_reason
            Parameters which indicate if the positioner is in conflict.
        r_fibre_focal, theta_fibre_focal
            Focal plane polar coordinates of the fibre.
        x_fibre_focal, y_fibre_focal
            Focal plane Cartesian coordinates of the fibre.
        r_fibre_local, theta_fibre_focal
            Local polar coordinates of the fibre.
        x_fibre_local, y_fibre_local
            Local Cartesian coordinates of the fibre.
        parity
            Elbow parity of the positioner orientation.
        x_elbow_local, y_elbow_local
            Local Cartesian coordinates of the beta axis (or elbow).

        :Returns:

        target_reachable: boolean
            True if the target can be reached (target location stored)
            False if the target cannot be reached (target ignored)

        """
        # Laleh: It would be easier for PAS to handle targets already presented in local frame of each FPU
        # Can the positioner reach this target?
        if self.can_move_to(r_target_focal, theta_target_focal, parity):

            # Yes, target can be reached. If this is not the starting position,
            # and the positioner is not locked, the previous target is
            # overwritten.
            if not starting_location:
                if not self.locked:
                    self.target_assigned = True
                    self.r_target_focal = r_target_focal
                    self.theta_target_focal = theta_target_focal
                    self.target_parity = parity
                else:
                    # The positioner cannot be moved from its starting location.
                    self.in_conflict = True
                    self.conflict_type = CONFLICT_LOCKED
                    self.conflict_reason = "%s: Positioner locked." % self.name
                    return False

            # Update the fibre location to match the target, in focal plane
            # polar and Cartesian coordinates.
            self.r_fibre_focal = r_target_focal
            self.theta_fibre_focal = theta_target_focal
            (self.x_fibre_focal, self.y_fibre_focal) = \
                util.polar_to_cartesian(self.r_fibre_focal,
                                        self.theta_fibre_focal)
            strg = "Target set to focal coordinates " + \
                "(r=%f,theta=%f) or (x=%f,y=%f)" % \
                (self.r_fibre_focal, math.degrees(self.theta_fibre_focal),
                 self.x_fibre_focal, self.y_fibre_focal)
            logger.debug(strg)

            # Determine the local (wrt positioner centre) coordinates for the
            # target.
            self.x_fibre_local = self.x_fibre_focal - self.x_centre_focal
            self.y_fibre_local = self.y_fibre_focal - self.y_centre_focal
            (self.r_fibre_local, self.theta_fibre_local) = \
                util.cartesian_to_polar(self.x_fibre_local, self.y_fibre_local)
            self.parity = parity
            strg = "Target is at local coordinates " + \
                "(r=%f,theta=%f) or (x=%f,y=%f) at %s parity" % \
                (self.r_fibre_local, math.degrees(self.theta_fibre_local),
                 self.x_fibre_local, self.y_fibre_local,
                 util.elbow_parity_str(self.parity))
            logger.debug(strg)

            # Solve the elbow location for this target and save the result.
            (self.x_elbow_local, self.y_elbow_local) = \
                util.solve_elbow_xy(self.x_fibre_local, self.y_fibre_local,
                             self.parity, self.length1, self._length1sq,
                             self.length2, self._length2sq)
            self.in_conflict = False
            return True
        else:
            # No, target cannot be reached. The previous target is not overwritten.
            return False

    # ***
    def set_arm_angles(self, angle1, angle2):
        """

        Set the alpha and beta motor angles for this positioner.
        The simulated positioner is moved immediately to these angles.
        A new target is NOT defined.

        :Parameters:

        angle1: float
            Alpha motor angle, with respect to its datum location (in degrees)
        angle2: float
            Beta motor angle, with respect to its datum location (in degrees)

        :Attributes:

        The function sets the following attributes.

        parity
            Elbow parity of the positioner orientation.
        x_elbow_local, y_elbow_local
            Local Cartesian coordinates of the beta axis (or elbow).
        x_fibre_local, y_fibre_local
            Local Cartesian coordinates of the fibre.
        x_fibre_focal, y_fibre_focal
            Focal plane Cartesian coordinates of the fibre.
        r_fibre_local, theta_fibre_focal
            Local polar coordinates of the fibre.
        r_fibre_focal, theta_fibre_focal
            Focal plane polar coordinates of the fibre.

        :Returns:

        None

        """
        # Make sure alpha and beta angles are within the defined range
        # Ensure that the requested angles are within their allowed limits.
        while angle1 < ALPHA_LIMITS[0] :
            angle1 += 360.0
        while angle1 > ALPHA_LIMITS[1]:
            angle1 -= 360.0
        while angle2 < BETA_LIMITS[0] :
            angle2 += 360.0
        while angle2 > BETA_LIMITS[1]:
            angle2 -= 360.0
        angle1_rad = math.radians(angle1)
        angle2_rad = math.radians(angle2)

        if angle2 > 0.0:
            self.parity = PARITY_LEFT
        else:
            self.parity = PARITY_RIGHT

        alpha_angle = angle1_rad + self.orient
        self.x_elbow_local = self.length1 * math.cos(alpha_angle)
        self.y_elbow_local = self.length1 * math.sin(alpha_angle)

        elbow_angle = math.pi + alpha_angle + angle2_rad
        self.x_fibre_local = self.x_elbow_local + self.length2 * math.cos(elbow_angle)
        self.y_fibre_local = self.y_elbow_local + self.length2 * math.sin(elbow_angle)

        (self.r_fibre_local, self.theta_fibre_local) = \
            util.cartesian_to_polar(self.x_fibre_local, self.y_fibre_local)
        strg = "Arms moved to (alpha=%f,beta=%f) at orient=%f\n" % \
            (angle1,angle2, math.degrees(self.orient))
        strg += "  which moves fibre to local coordinates " + \
            "(r=%f,theta=%f) or (x=%f,y=%f) at %s parity" % \
            (self.r_fibre_local, math.degrees(self.theta_fibre_local),
                self.x_fibre_local, self.y_fibre_local,
                util.elbow_parity_str(self.parity))
        logger.debug(strg)

        self.x_fibre_focal = self.x_centre_focal + self.x_fibre_local
        self.y_fibre_focal = self.y_centre_focal + self.y_fibre_local
        (self.r_fibre_focal, self.theta_fibre_focal) = \
            util.cartesian_to_polar(self.x_fibre_focal, self.y_fibre_focal)
        strg = "Fibre is at focal coordinates " + \
                "(r=%f,theta=%f) or (x=%f,y=%f)" % \
                (self.r_fibre_focal, math.degrees(self.theta_fibre_focal),
                 self.x_fibre_focal, self.y_fibre_focal)
        logger.debug(strg)

    # ***
    def get_arm_angles(self):
        """

        Calculate the alpha and beta motor angles that would be adopted for
        the current positioner target and parity.

        :Parameters:

        None

        :Returns:

        (angle1,angle2)
            Alpha and beta motor angles in degrees.

        """
        reachsq = self.r_fibre_local * self.r_fibre_local
        # Determine the third length of the triangle that needs
        # to be made by the arm's shoulder and elbow to reach
        # the fibre position.
        # The EPS in the following tests accounts for floating point rounding.
        if self.r_fibre_local > (self.outer - EPS):
            # Special case of arm at full stretch.
            angle1 = math.radians(90.0) - self.theta_fibre_local
            angle2 = -math.radians(180.0)
        elif self.r_fibre_local < (self.inner + EPS):
            # Special case of arm completely doubled back.
            angle1 =  (math.radians(90.0) - self.theta_fibre_local) + \
                        math.radians(180.0)
            angle2 = 0.0
        else:
            # Solve the triangle cosine rule to determine the angle between
            # the shoulder and the fibre target.
            shoulder_fibre = util.solve_triangle(self.r_fibre_local, reachsq,
                                self.length1, self._length1sq, self._length2sq )
            if shoulder_fibre is None:
                # Equation cannot be solved
                logger.error("+++ Shoulder to fibre angle equation cannot be solved")
                return (None, None)
#             print("shoulder_fibre=", math.degrees(shoulder_fibre))

            # Convert shoulder to fibre angle into alpha angle
            angle1 = util.solve_shoulder(self.theta_fibre_local, shoulder_fibre,
                                    self.parity)
#             print("angle1=", math.degrees(angle1))

            # Solve the cosine rule again to determine the angle between
            # the shoulder and elbow.
            shoulder_elbow = util.solve_triangle(self.length1, self._length1sq,
                                            self.length2, self._length2sq,
                                            reachsq)
            if shoulder_elbow is None:
                # Equation cannot be solved
                logger.error("+++ Shoulder to elbow angle equation cannot be solved")
                return (None, None)

            # The beta angle is the negative of this angle. The order of
            # calculation depends on the parity adopted.
            if self.parity == PARITY_RIGHT:
                angle2 = -shoulder_elbow
            else:
                angle2 = shoulder_elbow

        # Subtract the alpha zero-point angle
        angle1 -= self.orient

        # Ensure that the calculated angles are within their allowed limits.
        angle1_deg = math.degrees(angle1)
        angle2_deg = math.degrees(angle2)
        while angle1_deg < ALPHA_LIMITS[0] :
            angle1_deg += 360.0
        while angle1_deg > ALPHA_LIMITS[1]:
            angle1_deg -= 360.0
        while angle2_deg < BETA_LIMITS[0] :
            angle2_deg += 360.0
        while angle2_deg > BETA_LIMITS[1]:
            angle2_deg -= 360.0

        return (angle1_deg, angle2_deg)

    # ***
#     def difficulty_value(self, others=[],
#                          weights=[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]):
#     def difficulty_value(self, others=[],
#                          weights=[0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]):
#     def difficulty_value(self, others=[],
#                          weights=[0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]):
#     def difficulty_value(self, others=[],
#                          weights=[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]):
#     def difficulty_value(self, others=[],
#                          weights=[0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]):
#     def difficulty_value(self, others=[],
#                          weights=[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]):
#     def difficulty_value(self, others=[],
#                          weights=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]):
    def difficulty_value(self, others=[],
                         weights=[0.3, 0.0, 0.1, 0.0, 0.2, 0.2, 0.2]):
        """

        Calculate the difficulty value associated with the current positioner
        target and parity, given the targets assigned to the positioner's
        neighbours.

        :Parameters:

        others: list of Positioner (optional)
            A list of other positioners to be checked. If not specified,
            the positioner's neighbours will be checked.
        weights: list of 6 floats (optional)
            The weights to be applied to the difficulty calculation.

            * Alpha movement weight
            * Beta movement weight
            * Positioner extension/reach weight
            * Closest target-target approach weight
            * Closest metrology-metrology approach weight
            * Closest beta arm-beta arm approach weight
            * Closest beta arm-metrology approach weight

            The weights must add up to 1.

        :Returns:

        difficulty: float
            A number between 0.0 and 1.0, estimating the difficulty in
            achieving the given configuration with real hardware.

        """
        # TODO: INVESTIGATE THE BEST DIFFICULTY DEFINITION.

        # (1) Calculate the amount of alpha and beta axis movement required
        # to reach the target.
        (angle1,angle2) = self.get_arm_angles()
        difficulty_alpha = abs(angle1 - ALPHA_DEFAULT) / 360.0
        difficulty_beta = abs(angle2 - BETA_DEFAULT) / 360.0

        # (2) Calculate the degree to which this positioner's target intrudes
        # into the patrol field of a neighbouring positioner.
        # NOTE: Overlaps somewhat with the beta difficulty.
        reachsq = self.r_fibre_local * self.r_fibre_local
        middlesq = self._outersq / 4.0
        difficulty_reach = max(0.0, (reachsq - middlesq) / self._outersq)

        # (3) Calculate the closest distance between any two fibre
        # holders, belonging to this positioner and its neighbours, taking
        # into account the minimum target approach distance, INS_POS_MINDIST.
        closestsq = self._outersq
        if others:
            check_list = others
        else:
            check_list = self.near_neighbours + self.far_neighbours
        for neighbour in check_list:
            if isinstance(neighbour, Positioner):
                distsq = util.distance_squared(self.x_fibre_focal,
                                               self.y_fibre_focal,
                                               neighbour.x_fibre_focal,
                                               neighbour.y_fibre_focal)
                if distsq < closestsq:
                    closestsq = distsq
        # Base the difficulty on the square of the distance ratio, so the
        # value rises more steeply at close distances.
        ratio = max(0.0, (math.sqrt(closestsq) - INS_POS_MINDIST)) / self.outer
        difficulty_target_target = max(0.0, 1.0 - math.sqrt(ratio))

        # (4) Calculate the closest approach distance between any two
        # metrology avoidance ellipses, any two beta avoidance ellipses
        # or any metrology to beta avoidance ellipse.
        closestsq_met = self._outersq
        closestsq_beta = self._outersq
        closestsq_cross = self._outersq
        if others:
            check_list = others
        else:
            check_list = self.near_neighbours + self.far_neighbours
        my_met_ellipse = self.avoidance_ellipse_metrology()
        my_met_x = my_met_ellipse[0]
        my_met_y = my_met_ellipse[1]
        my_beta_ellipse = self.avoidance_ellipse_beta()
        my_beta_x = my_beta_ellipse[0]
        my_beta_y = my_beta_ellipse[1]
        # Also determine whether any beta arms cross each other.
        beta_arms_cross = False
        my_x_elbow_focal = self.x_centre_focal + self.x_elbow_local
        my_y_elbow_focal = self.y_centre_focal + self.y_elbow_local
        my_beta_arm = [my_x_elbow_focal, my_y_elbow_focal,
                       self.x_fibre_focal, self.y_fibre_focal]
        for neighbour in check_list:
            if isinstance(neighbour, Positioner):
                other_met_ellipse = neighbour.avoidance_ellipse_metrology()
                other_met_x = other_met_ellipse[0]
                other_met_y = other_met_ellipse[1]
                other_beta_ellipse = neighbour.avoidance_ellipse_beta()
                other_beta_x = other_beta_ellipse[0]
                other_beta_y = other_beta_ellipse[1]
                distsq_met = util.distance_squared(my_met_x, my_met_y,
                                                   other_met_x, other_met_y)
                distsq_beta = util.distance_squared(my_beta_x, my_beta_y,
                                                   other_beta_x, other_beta_y)
                if distsq_met < closestsq_met:
                    closestsq_met = distsq_met
                if distsq_beta < closestsq_beta:
                    closestsq_beta = distsq_beta
                distsq_cross = util.distance_squared(my_beta_x, my_beta_y,
                                                   other_met_x, other_met_y)
                if distsq_cross < closestsq_cross:
                    closestsq_cross = distsq_cross
                distsq_cross = util.distance_squared(my_met_x, my_met_y,
                                                   other_beta_x, other_beta_y)
                if distsq_cross < closestsq_cross:
                    closestsq_cross = distsq_cross

                other_x_elbow_focal = neighbour.x_centre_focal + neighbour.x_elbow_local
                other_y_elbow_focal = neighbour.y_centre_focal + neighbour.y_elbow_local
                other_beta_arm = [other_x_elbow_focal, other_y_elbow_focal,
                                  neighbour.x_fibre_focal, neighbour.y_fibre_focal]
                if util.lines_intersect(my_beta_arm, other_beta_arm):
                    beta_arms_cross = True

        # Base the difficulties on the square root of the distance ratio, so the
        # value rises more steeply at close distances.
        ratio = max(0.0, (math.sqrt(closestsq_met) - INS_POS_WIDTH2)) / self.outer
        difficulty_ellipse_met = max(0.0, 1.0 - math.sqrt(ratio))
        ratio = max(0.0, (math.sqrt(closestsq_beta) - INS_POS_MINDIST)) / self.outer
        difficulty_ellipse_beta = max(0.0, 1.0 - math.sqrt(ratio))
        ratio = max(0.0, math.sqrt(closestsq_cross)) / self.outer
        difficulty_ellipse_cross = max(0.0, 1.0 - math.sqrt(ratio))

        # Try varying these weights. The total must add up to 1.0
        if not isinstance(weights, (tuple,list)) or len(weights) < 7:
            weights=[0.2, 0.0, 0.1, 0.1, 0.1, 0.3, 0.2]
        difficulty_alpha_weight = weights[0]
        difficulty_beta_weight = weights[1]
        difficulty_reach_weight = weights[2]
        difficulty_target_target_weight = weights[3]
        difficulty_ellipse_met_weight = weights[4]
        difficulty_ellipse_beta_weight = weights[5]
        difficulty_ellipse_cross_weight = weights[6]

        # Return a combined difficulty estimate.
        strg = "Target/parity difficulties for %d: " % self.ident
        strg += "alpha=%f, beta=%f, reach=%f, target=%f, ellipse_met=%f, ellipse_bet=%f, ellipse_cross=%f" % \
            (difficulty_alpha, difficulty_beta, difficulty_reach,
             difficulty_target_target, difficulty_ellipse_met,
             difficulty_ellipse_beta, difficulty_ellipse_cross)
        logger.debug( strg )
        difficulty = difficulty_alpha_weight * difficulty_alpha + \
                difficulty_beta_weight * difficulty_beta + \
                difficulty_reach_weight * difficulty_reach + \
                difficulty_target_target_weight * difficulty_target_target + \
                difficulty_ellipse_met_weight * difficulty_ellipse_met + \
                difficulty_ellipse_beta_weight * difficulty_ellipse_beta + \
                difficulty_ellipse_cross_weight * difficulty_ellipse_cross

        # Increase the difficulty factor if the beta arms cross.
        if beta_arms_cross:
            difficulty += (1.0 - difficulty) / 2.0

        return difficulty

    def avoidance_polygon_metrology(self, padding=0.0):
        """

        Return, in focal plane coordinates, the vertices of the
        polygon bounding the main metrology target area,
        defined as ((x1, y2), (x2, y2), (x3, y3), (x4, y4), (x5, y5)).

        :Parameters:

        padding: float (optional)
            Safety padding added around the edge of the avoidance zone.
            Defaults to 0.0

        :Returns:

        ((x1, y2), (x2, y2), (x3, y3), (x4, y4), (x5, y5))

        """
        # TODO: These avoidance zone functions need to be optimised.
        # Calculate the vertices of a rectangle which encloses the fibre
        # holder and the metrology targets.
        xdiff = self.x_fibre_local - self.x_elbow_local
        ydiff = self.y_fibre_local - self.y_elbow_local

        half_width = padding + INS_POS_W4/2.0
        xdelta = half_width * ydiff / self.length2
        ydelta = half_width * xdiff / self.length2

        metrology_fraction = INS_POS_B1 / self.length2
        padding_fraction = padding / self.length2
        xmetrology_end = self.x_fibre_focal - metrology_fraction * xdiff
        ymetrology_end = self.y_fibre_focal - metrology_fraction * ydiff

        triangle_fraction = INS_POS_TB2 / self.length2
        xtriangle_end = xmetrology_end - \
            (triangle_fraction + padding_fraction) * xdiff
        ytriangle_end = ymetrology_end - \
            (triangle_fraction + padding_fraction) * ydiff
        xmetrology_start = self.x_fibre_focal + padding_fraction * xdiff
        ymetrology_start = self.y_fibre_focal + padding_fraction * ydiff

        x1 = xmetrology_end - xdelta
        y1 = ymetrology_end + ydelta
        x2 = xtriangle_end
        y2 = ytriangle_end
        x3 = xmetrology_end + xdelta
        y3 = ymetrology_end - ydelta
        x4 = xmetrology_start + xdelta
        y4 = ymetrology_start - ydelta
        x5 = xmetrology_start - xdelta
        y5 = ymetrology_start + ydelta
        return ((x1,y1),(x2,y2),(x3,y3),(x4,y4),(x5,y5))

    def avoidance_rectangle_metrology(self, padding=0.0):
        """

        Return, in focal plane coordinates, the vertices of the
        rectangle bounding the main metrology target area,
        defined as ((x1, y2), (x2, y2), (x3, y3), (x4, y4)).

        :Parameters:

        padding: float (optional)
            Safety padding added around the edge of the avoidance zone.
            Defaults to 0.0

        :Returns:

        ((x1, y2), (x2, y2), (x3, y3), (x4, y4))

        """
        # TODO: These avoidance zone functions need to be optimised.
        # TODO: REPLACE WITH AVOIDANCE_POLYGON_METROLOGY
        # Calculate the vertices of a rectangle which encloses the fibre
        # holder and the metrology targets.
        xdiff = self.x_fibre_local - self.x_elbow_local
        ydiff = self.y_fibre_local - self.y_elbow_local

        half_width = padding + INS_POS_W4/2.0
        xdelta = half_width * ydiff / self.length2
        ydelta = half_width * xdiff / self.length2

        metrology_fraction = INS_POS_B1 / self.length2
        padding_fraction = padding / self.length2
        xmetrology_end = self.x_fibre_focal - \
            (metrology_fraction + padding_fraction) * xdiff
        ymetrology_end = self.y_fibre_focal - \
            (metrology_fraction + padding_fraction) * ydiff
        xmetrology_start = self.x_fibre_focal + padding_fraction * xdiff
        ymetrology_start = self.y_fibre_focal + padding_fraction * ydiff

        x1 = xmetrology_end - xdelta
        y1 = ymetrology_end + ydelta
        x2 = xmetrology_end + xdelta
        y2 = ymetrology_end - ydelta
        x3 = xmetrology_start + xdelta
        y3 = ymetrology_start - ydelta
        x4 = xmetrology_start - xdelta
        y4 = ymetrology_start + ydelta
        return ((x1,y1),(x2,y2),(x3,y3),(x4,y4))

    def avoidance_triangle_metrology(self, padding=0.0):
        """

        Return, in focal plane coordinates, the vertices of the
        triangle bounding the extension of the metrology target area,
        defined as ((x1, y2), (x2, y2), (x3, y3)).

        :Parameters:

        padding: float (optional)
            Safety padding added around the edge of the avoidance zone.
            Defaults to 0.0

        :Returns:

        ((x1, y2), (x2, y2), (x3, y3))

        """
        # TODO: These avoidance zone functions need to be optimised.
        # TODO: REPLACE WITH AVOIDANCE_POLYGON_METROLOGY
        # Calculate the vertices of a triangle bounding the extension
        # of the metrology target area.
        xdiff = self.x_fibre_local - self.x_elbow_local
        ydiff = self.y_fibre_local - self.y_elbow_local

        half_width = padding + INS_POS_W4/2.0
        xdelta = half_width * ydiff / self.length2
        ydelta = half_width * xdiff / self.length2

        metrology_fraction = INS_POS_B1 / self.length2
        padding_fraction = padding / self.length2
        xmetrology_end = self.x_fibre_focal - metrology_fraction * xdiff
        ymetrology_end = self.y_fibre_focal - metrology_fraction * ydiff

        triangle_fraction = (INS_POS_TB2+padding) / self.length2
        xtriangle_end = xmetrology_end - \
            (triangle_fraction+padding_fraction) * xdiff
        ytriangle_end = ymetrology_end - \
            (triangle_fraction+padding_fraction) * ydiff

        x1 = xmetrology_end - xdelta
        y1 = ymetrology_end + ydelta
        x2 = xmetrology_end + xdelta
        y2 = ymetrology_end - ydelta
        x3 = xtriangle_end
        y3 = ytriangle_end
        return ((x1,y1),(x2,y2),(x3,y3))

    def avoidance_ellipse_metrology(self, padding=0.0):
        """

        Return, in focal plane coordinates, the parameters of
        an ellipse covering the metrology avoidance zone.

        :Parameters:

        padding: float (optional)
            Safety padding added around the edge of the avoidance zone.
            Defaults to 0.0

        :Returns:

        (xcen, ycen, major, minor, tilt)
            The centre, major axis, minor axis and tilt angle of the ellipse

        """
        # TODO: These avoidance zone functions need to be optimised.
        # Determine the extent of the metrology target zone.
        # Expand the zone slightly in the direction of the fibre holder.
        xdiff = self.x_fibre_local - self.x_elbow_local
        ydiff = self.y_fibre_local - self.y_elbow_local
        metrology_fraction = (INS_POS_B1+INS_POS_TB2+padding) / self.length2
        fibre_holder_fraction = 1.15 * (FIBRE_RADIUS+padding) / self.length2
        xmetrology_end = self.x_fibre_focal - metrology_fraction * xdiff
        ymetrology_end = self.y_fibre_focal - metrology_fraction * ydiff
        xmetrology_start = self.x_fibre_focal + fibre_holder_fraction * xdiff
        ymetrology_start = self.y_fibre_focal + fibre_holder_fraction * ydiff

        # The ellipse is placed in the middle of this zone
        xcen = (xmetrology_start + xmetrology_end) / 2.0
        ycen = (ymetrology_start + ymetrology_end) / 2.0

        # The major and minor axes are determined by the length and width
        # of the zone. Over-size the width a little to ensure the zone
        # is mostly covered by the ellipse.
        major = (FIBRE_RADIUS+INS_POS_B1+INS_POS_TB2)/2.0
        minor = 1.25 * INS_POS_W4 /  2.0
        tilt = math.atan2(ydiff, xdiff)

        return (xcen, ycen, major, minor, tilt)

    def avoidance_triangle_fibre(self, padding=0.0):
        """

        Return, in focal plane coordinates, the vertices of the
        triangle bounding the fibre avoidance area,
        defined as ((x1, y2), (x2, y2), (x3, y3)).

        :Parameters:

        padding: float (optional)
            Safety padding added around the edge of the avoidance zone.
            Defaults to 0.0

        :Returns:

        ((x1, y2), (x2, y2), (x3, y3))

        """
        # TODO: These avoidance zone functions need to be optimised.
        # Calculate the vertices of a triangle bounding the the fibre
        # avoidance area.
        xdiff = self.x_fibre_local - self.x_elbow_local
        ydiff = self.y_fibre_local - self.y_elbow_local

        half_width = padding + INS_POS_W3/2.0
        xdelta = half_width * ydiff / self.length2
        ydelta = half_width * xdiff / self.length2

        metrology_fraction = INS_POS_B1 / self.length2
        padding_fraction = padding / self.length2
        # The padding is extended at the pointed end of the triangle
        # and reduced at the straight end by this ratio.
        triang_length_ratio = INS_POS_B1 / INS_POS_W3
        xmetrology_end = self.x_fibre_focal - metrology_fraction * xdiff
        ymetrology_end = self.y_fibre_focal - metrology_fraction * ydiff
        xtriangle_start = xmetrology_end + \
            padding_fraction * xdiff/triang_length_ratio
        ytriangle_start = ymetrology_end + \
            padding_fraction * ydiff/triang_length_ratio

        triangle_fraction = INS_POS_B3 / self.length2
        xtriangle_end = xmetrology_end - \
            (triangle_fraction+padding_fraction*triang_length_ratio) * xdiff
        ytriangle_end = ymetrology_end - \
            (triangle_fraction+padding_fraction*triang_length_ratio) * ydiff

        x1 = xtriangle_start - xdelta
        y1 = ytriangle_start + ydelta
        x2 = xtriangle_start + xdelta
        y2 = ytriangle_start - ydelta
        x3 = xtriangle_end
        y3 = ytriangle_end
        return ((x1,y1),(x2,y2),(x3,y3))

    def avoidance_rectangles_beta(self, nzones=3, padding=0.0):
        """

        Return, in focal plane cooordinates, the vertices of the rectangles
        representing successive zones along the beta arm. These zones
        approximate the curved, 3-D design of the beta arm as a "stair case"
        of successive cuboids at different heights.

        Two beta arms will clash if two corresponding zones (at the same
        height) overlap.

        :Parameters:

        nzones: int (optional)
            The number of zones in which to divide the vulnerable
            portion of the beta arm. Default 3.
        padding: float (optional)
            Safety padding added around the edge of the avoidance zone.
            Defaults to 0.0

        :Returns:

        list of rectangles: list of ((x1,y1),(x2,y2),(x3,y3),(x4,y4)))
            The beta arm avoidance rectangles, in focal plane coordinates.

        """
        # TODO: These avoidance zone functions need to be optimised.
        rect_list = []

        # Calculate X and Y increments for the corners of each rectangle
        # from the length and orientation of the beta arm.
        xdiff = self.x_fibre_local - self.x_elbow_local
        ydiff = self.y_fibre_local - self.y_elbow_local
        half_width = padding + INS_POS_WIDTH2/2.0
        xdelta = half_width * ydiff / self.length2
        ydelta = half_width * xdiff / self.length2

        # Find the start of the vulnerable zone on the beta arm.
        padding_fraction = padding / self.length2
        safe_fraction = self.bsafe / self.length2
        xbeta_start = self.x_centre_focal + self.x_elbow_local + \
                        (safe_fraction-padding_fraction) * xdiff
        ybeta_start = self.y_centre_focal + self.y_elbow_local + \
                        (safe_fraction-padding_fraction) * ydiff

        # Divide the vulnerable zone into nzones rectangles
        zone_fraction = (self.length2 + 2*padding - self.bsafe - \
                         INS_POS_B1 - INS_POS_B2) / \
                        ( nzones * self.length2)
        xfraction = xdiff * zone_fraction
        yfraction = ydiff * zone_fraction
        xzone_start = xbeta_start
        yzone_start = ybeta_start
        for zone in range(0, nzones):
            xzone_end = xzone_start + xfraction
            yzone_end = yzone_start + yfraction
            x1 = xzone_end - xdelta
            y1 = yzone_end + ydelta
            x2 = xzone_end + xdelta
            y2 = yzone_end - ydelta
            x3 = xzone_start + xdelta
            y3 = yzone_start - ydelta
            x4 = xzone_start - xdelta
            y4 = yzone_start + ydelta
            rect_list.append( ((x1,y1),(x2,y2),(x3,y3),(x4,y4)) )
            xzone_start = xzone_end
            yzone_start = yzone_end

        return rect_list

    def avoidance_ellipse_beta(self, padding=0.0):
        """

        Return, in focal plane coordinates, the parameters of
        an ellipse covering the beta arm avoidance zone.

        :Parameters:

        padding: float (optional)
            Safety padding added around the edge of the avoidance zone.
            Defaults to 0.0

        :Returns:

        (xcen, ycen, major, minor, tilt)
            The centre, major axis, minor axis and tilt angle of the ellipse

        """
        # TODO: These avoidance zone functions need to be optimised.
        # Determine the extent of the beta arm zone.
        xdiff = self.x_fibre_local - self.x_elbow_local
        ydiff = self.y_fibre_local - self.y_elbow_local

        # Find the start of the vulnerable zone on the beta arm.
        # Expand the zone slightly in the direction of the beta axis.
        safe_fraction = 0.85 * (self.bsafe-padding) / self.length2
        xbeta_start = self.x_centre_focal + self.x_elbow_local + \
                        safe_fraction * xdiff
        ybeta_start = self.y_centre_focal + self.y_elbow_local + \
                        safe_fraction * ydiff

        zone_length = padding + self.length2 - self.bsafe - INS_POS_B1 - INS_POS_B2
        zone_fraction = zone_length / self.length2
        xbeta_end = xbeta_start + xdiff * zone_fraction
        ybeta_end = ybeta_start + ydiff * zone_fraction

        # The ellipse is placed in the middle of this zone
        xcen = (xbeta_start + xbeta_end) / 2.0
        ycen = (ybeta_start + ybeta_end) / 2.0

        # The major and minor axes are determined by the length and width
        # of the zone. Over-size the length a little to ensure the zone
        # is mostly covered by the ellipse.
        major = 1.1 * (self.length2 - self.bsafe - INS_POS_B1 - INS_POS_B2)/2.0
        minor = INS_POS_WIDTH2 /  2.0
        tilt = math.atan2(ydiff, xdiff)

        return (xcen, ycen, major, minor, tilt)

    def avoidance_zone_datum(self, padding=0.0):
        """

        Return, in focal plane cooordinates, the vertices of the rectangle
        connecting the circular "datum foot" to the axis of the beta arm.

        Two beta arms will clash if their datum rectangles overlap.

        :Parameters:

        padding: float (optional)
            Safety padding added around the edge of the avoidance zone.
            Defaults to 0.0

        :Returns:

        (rect, circle)
            rect =((x1, y2), (x2, y2), (x3, y3), (x4, y4))
                The rectangular attachment to the "datum foot".
            circle = (xcen, ycen, radius)
                The centre and radius of the circular "datum foot" zone.

        """
        # TODO: These avoidance zone functions need to be optimised.
        shape_list = []

        # Calculate X and Y increments for the corners of each rectangle
        # from the length and orientation of the beta arm.
        xdiff = self.x_fibre_local - self.x_elbow_local
        ydiff = self.y_fibre_local - self.y_elbow_local
        half_width = padding + INS_POS_DW/2.0
        xdelta = half_width * ydiff / self.length2
        ydelta = half_width * xdiff / self.length2

        # The datum actuator starts at the beta axis (elbow) but
        # extends beyond it (because it contains the beta axis).
        padding_fraction = padding / self.length2
        xdatum_start = self.x_centre_focal + self.x_elbow_local - \
            padding_fraction * xdiff
        ydatum_start = self.y_centre_focal + self.y_elbow_local - \
            padding_fraction * ydiff
        # Extend the start beyond the beta axis
        xdatum_start = xdatum_start - xdiff * INS_POS_D2 / self.length2
        ydatum_start = ydatum_start - ydiff * INS_POS_D2 / self.length2

        # The zone ends at the datum switch above the alpha axis, so it is
        # the alpha arm length from the beta axis.
        datum_fraction = INS_POS_D1 / self.length2
        xdatum_circ = self.x_centre_focal + self.x_elbow_local + datum_fraction * xdiff
        ydatum_circ = self.y_centre_focal + self.y_elbow_local + datum_fraction * ydiff
        xdatum_end = xdatum_circ + padding_fraction * xdiff
        ydatum_end = ydatum_circ + padding_fraction * ydiff

        x1 = xdatum_end-xdelta
        y1 = ydatum_end+ydelta
        x2 = xdatum_end+xdelta
        y2 = ydatum_end-ydelta
        x3 = xdatum_start+xdelta
        y3 = ydatum_start-ydelta
        x4 = xdatum_start-xdelta
        y4 = ydatum_start+ydelta
        rect = ((x1,y1),(x2,y2),(x3,y3),(x4,y4))
        shape_list.append( rect )

        # Assume the radius of the "foot" is the same as the
        # half-width of the beta arm.
        # radius = padding + INS_POS_DW/2.0
        # THE DATUM CIRCLE HAS GONE.  REPLACE FOR NOW WITH A TINY CIRCLE.
        radius = 0.1
        circle = (xdatum_circ, ydatum_circ, radius)
        shape_list.append( circle )
        return shape_list

    # ***
    def in_conflict_with(self, other, safety=INS_POS_SAFETY, symmetrical=True):
        """

        Determine whether this positioner is in conflict with another,
        given their current targets.

        :Parameters:

        other: Positioner, AcquisitionCamera or FiducialMarker object
            The other object with which conflict is to be checked.
        safety: float (optional)
            All conflict avoidance zones are padded with this safety
            tolerance. The default is the design tolerance, INS_POS_SAFETY.
        symmetrical: boolean (optional)
            Set to True to make all conflict checks symmetrical, so that when
            one positioner is in conflict both are in conflict. This only makes
            a difference to the fibre snag check.
            The default is True - when there is a fibre snag, both positioners
            are considered to be in conflict.

        :Attributes:

        The function sets the following attributes.

        in_conflict, conflict_type, conflict_reason
            Parameters which indicate if the positioner is in conflict.

        :Returns:

        in_conflict: boolean
            True if the two positioners are in conflict
            False if the positioners are not in conflict

        """
        if isinstance(other, Positioner):
            # Positioner/Positioner conflict check

            # Determine the distance between the fibre targets associated with
            # the two positioners.
            xdiff = self.x_fibre_focal - other.x_fibre_focal
            ydiff = self.y_fibre_focal - other.y_fibre_focal
            distsq = xdiff * xdiff + ydiff * ydiff

            # Are the targets sufficiently far apart (further than the maximum
            # reach of this positioner plus minimum separation) that a conflict
            # is not possible?
            if distsq > self._safedistsq:
                return False

            # Are the targets closer than the minimum separation distance
            # between fibre holders?
            mindistsq = (INS_POS_MINDIST+2*safety) * (INS_POS_MINDIST+2*safety)
            if distsq <= mindistsq:
                self.in_conflict = True
                self.conflict_type = CONFLICT_TOOCLOSE
                self.conflict_reason = "Conflict between %s and %s: " % \
                    (self.name, other.name)
                self.conflict_reason += "Targets too close: dist^2 %.3f <= %.3f" % \
                    (distsq, mindistsq)
                return True

            # Do the rectangular metrology zones of the two positioners overlap?
            this_metrology_rectangle = \
                self.avoidance_rectangle_metrology(padding=safety)
            other_metrology_rectangle = \
                other.avoidance_rectangle_metrology(padding=safety)
            logger.debug("My metrology rectangle: " + str(this_metrology_rectangle))
            logger.debug("Other metrology rectangle: " + str(other_metrology_rectangle))
#             if util.quadrangles_intersect(this_metrology_rectangle,
#                                           other_metrology_rectangle):
            if util.polygons_intersect(this_metrology_rectangle,
                                          other_metrology_rectangle):
                self.in_conflict = True
                self.conflict_type = CONFLICT_METROLOGY
                self.conflict_reason = "Conflict between %s and %s: " % \
                    (self.name, other.name)
                self.conflict_reason += "Metrology rectangles intersect"
                return True

            # Do the triangular metrology zones of the two positioners overlap?
            this_metrology_triangle = \
                self.avoidance_triangle_metrology(padding=safety)
            other_metrology_triangle = \
                other.avoidance_triangle_metrology(padding=safety)
            logger.debug("My metrology triangle: " + str(this_metrology_triangle))
            logger.debug("Other metrology triangle: " + str(other_metrology_triangle))
#             if util.triangles_intersect(this_metrology_triangle,
#                                         other_metrology_triangle):
            if util.polygons_intersect(this_metrology_triangle,
                                       other_metrology_triangle):
                self.in_conflict = True
                self.conflict_type = CONFLICT_METROLOGY
                self.conflict_reason = "Conflict between %s and %s: " % \
                    (self.name, other.name)
                self.conflict_reason += "Metrology triangles intersect"
                return True
#             if util.quadrangle_intersects_triangle(this_metrology_rectangle,
#                                         other_metrology_triangle):
            if util.polygons_intersect(this_metrology_rectangle,
                                       other_metrology_triangle):
                self.in_conflict = True
                self.conflict_type = CONFLICT_METROLOGY
                self.conflict_reason = "Conflict between %s and %s: " % \
                    (self.name, other.name)
                self.conflict_reason += "Metrology rectangle/triangle intersects"
                return True

            # Do the triangular and rectangular zones intersect?
#             if util.quadrangle_intersects_triangle(other_metrology_rectangle,
#                                         this_metrology_triangle):
            if util.polygons_intersect(other_metrology_rectangle,
                                        this_metrology_triangle):
                self.in_conflict = True
                self.conflict_type = CONFLICT_METROLOGY
                self.conflict_reason = "Conflict between %s and %s: " % \
                    (self.name, other.name)
                self.conflict_reason += "Metrology rectangle/triangle intersects"
                return True

#             # TODO: REPLACE WITH THE FOLLOWING?
#             this_metrology_polygon = self.avoidance_polygon_metrology()
#             other_polygon_polygon = other.avoidance_polygon_metrology()
#             if util.polygons_intersect(this_metrology_polygon, other_polygon_polygon):
#                 self.in_conflict = True
#                 self.conflict_type = CONFLICT_METROLOGY
#                 self.conflict_reason = "Metrology polygons intersect"
#                 return True


            # Does this positioner's fibre holder intersect with the other positioner's
            # metrology rectangle?
#             if util.quadrangle_intersects_circle(other_metrology_rectangle,
#                                                self.x_fibre_focal,
#                                                self.y_fibre_focal,
#                                                (FIBRE_RADIUS+safety)):
            if util.polygon_intersects_circle(other_metrology_rectangle,
                                               self.x_fibre_focal,
                                               self.y_fibre_focal,
                                               (FIBRE_RADIUS+safety)):
                self.in_conflict = True
                self.conflict_type = CONFLICT_METROLOGY
                self.conflict_reason = "Conflict between %s and %s: " % \
                    (self.name, other.name)
                self.conflict_reason += "Fibre holder intersects metrology rectangle"
                logger.debug(self.conflict_reason)
#                 if GRAPHICAL_DEBUGGING and plotting is not None:
#                     (xlist, ylist) = util.polygon_to_points(other_metrology_rectangle)
#                     plotaxis = plotting.plot_circles([self.x_fibre_focal],
#                         [self.y_fibre_focal], FIBRE_RADIUS,
#                         title="Fibre holder intersects metrology rectangle", showplot=False)
#                     plotting.plot_xy(xlist, ylist, plotaxis=plotaxis, showplot=True)
#                     plotting.close()
                return True

            # Does the other's fibre holder intersect with the this positioner's avoidance
            # rectangle (within the given safety tolerance)?
#             if util.quadrangle_intersects_circle(this_metrology_rectangle,
#                                            other.x_fibre_focal,
#                                            other.y_fibre_focal,
#                                            (FIBRE_RADIUS+safety)):
            if util.polygon_intersects_circle(this_metrology_rectangle,
                                           other.x_fibre_focal,
                                           other.y_fibre_focal,
                                           (FIBRE_RADIUS+safety)):
                self.in_conflict = True
                self.conflict_type = CONFLICT_METROLOGY
                self.conflict_reason = "Conflict between %s and %s: " % \
                (self.name, other.name)
                self.conflict_reason += "Other's fibre holder metrology rectangle"
                logger.debug(self.conflict_reason)
#                 if GRAPHICAL_DEBUGGING and plotting is not None:
#                     (xlist, ylist) = util.polygon_to_points(this_metrology_rectangle)
#                     plotaxis = plotting.plot_circles([other.x_fibre_focal],
#                         [other.y_fibre_focal], FIBRE_RADIUS,
#                         title="Other's fibre holder metrology rectangle", showplot=False)
#                     plotting.plot_xy(xlist, ylist, plotaxis=plotaxis, showplot=True)
#                     plotting.close()
                return True

            # Does this fibre holder intersect with the other positioner's avoidance
            # triangle (within the given safety tolerance)?
#             if util.triangle_intersects_circle(other_metrology_triangle,
#                                                self.x_fibre_focal,
#                                                self.y_fibre_focal,
#                                                (FIBRE_RADIUS+safety)):
            if util.polygon_intersects_circle(other_metrology_triangle,
                                               self.x_fibre_focal,
                                               self.y_fibre_focal,
                                               (FIBRE_RADIUS+safety)):
                self.in_conflict = True
                self.conflict_type = CONFLICT_METROLOGY
                self.conflict_reason = "Conflict between %s and %s: " % \
                    (self.name, other.name)
                self.conflict_reason += "Fibre holder intersects metrology triangle"
                logger.debug(self.conflict_reason)
#                 if GRAPHICAL_DEBUGGING and plotting is not None:
#                     (xlist, ylist) = util.polygon_to_points(other_metrology_triangle)
#                     plotaxis = plotting.plot_circles([self.x_fibre_focal],
#                         [self.y_fibre_focal], FIBRE_RADIUS,
#                         title="Fibre holder intersects metrology triangle", showplot=False)
#                     plotting.plot_xy(xlist, ylist, plotaxis=plotaxis, showplot=True)
#                     plotting.close()
                return True

            # Does the other's fibre holder intersect with the this positioner's avoidance
            # triangle (within the given safety tolerance)?
#             if util.triangle_intersects_circle(this_metrology_triangle,
#                                            other.x_fibre_focal,
#                                            other.y_fibre_focal,
#                                            (FIBRE_RADIUS+safety)):
            if util.polygon_intersects_circle(this_metrology_triangle,
                                           other.x_fibre_focal,
                                           other.y_fibre_focal,
                                           (FIBRE_RADIUS+safety)):
                self.in_conflict = True
                self.conflict_type = CONFLICT_FIBRE_SNAG
                self.conflict_reason = "Conflict between %s and %s: " % \
                (self.name, other.name)
                self.conflict_reason += "Other's fibre holder metrology triangle"
                logger.debug(self.conflict_reason)
#                 if GRAPHICAL_DEBUGGING and plotting is not None:
#                     (xlist, ylist) = util.polygon_to_points(this_metrology_triangle)
#                     plotaxis = plotting.plot_circles([other.x_fibre_focal],
#                         [other.y_fibre_focal], FIBRE_RADIUS,
#                         title="Other's fibre holder metrology triangle", showplot=False)
#                     plotting.plot_xy(xlist, ylist, plotaxis=plotaxis, showplot=True)
#                     plotting.close()
                return True


            # Does one of the beta arm zones overlap with the same zone on
            # the other arm?
            if INCLUDE_BETA_ZONE:
                this_rectangles_beta = \
                    self.avoidance_rectangles_beta(padding=safety)
                other_rectangles_beta = \
                    other.avoidance_rectangles_beta(padding=safety)
                # Test the rectangles one at a time.
                zone = 1
                for this_rect, other_rect in zip(this_rectangles_beta,
                                             other_rectangles_beta):
                    logger.debug("My beta rectangle " + str(zone) + ": " + str(this_rect))
                    logger.debug("Other beta rectangle " + str(zone) + ": " + str(other_rect))
#                     if util.quadrangles_intersect(this_rect, other_rect):
                    if util.polygons_intersect(this_rect, other_rect):
                        self.conflict_type = CONFLICT_BETA
                        self.in_conflict = True
                        self.conflict_reason = "Conflict between %s and %s: " % \
                            (self.name, other.name)
                        self.conflict_reason += \
                            "Beta arms intersect at same level (zone %d)" % zone
#                         if GRAPHICAL_DEBUGGING and plotting is not None:
#                             (xlist1, ylist1) = util.polygon_to_points(this_rect)
#                             (xlist2, ylist2) = util.polygon_to_points(other_rect)
#                             plotaxis = plotting.plot_xy(xlist1, ylist1,
#                                         title="Beta zone conflict", showplot=False)
#                             plotting.plot_xy(xlist2, ylist2, plotaxis=plotaxis,
#                                         showplot=True)
                        return True
                    zone = zone + 1

            if INCLUDE_DATUM_FOOT:
                (this_datum_rect, this_datum_circle) = \
                    self.avoidance_zone_datum(padding=safety)
                (other_datum_rect, other_datum_circle) = \
                    other.avoidance_zone_datum(padding=safety)
                logger.debug("My datum circle: " + str(this_datum_circle))
                logger.debug("Other datum circle: " + str(other_datum_circle))
                limitsq = this_datum_circle[2] * this_datum_circle[2] * 4.0
                if util.closer_than(this_datum_circle[0], this_datum_circle[1],
                                    other_datum_circle[0], other_datum_circle[1],
                                    limitsq):
                    self.in_conflict = True
                    self.conflict_type = CONFLICT_DATUM
                    self.conflict_reason = "Conflict between %s and %s: " % \
                            (self.name, other.name)
                    self.conflict_reason += \
                            "Datum actuators intersect"
#                     if GRAPHICAL_DEBUGGING and plotting is not None:
#                         xcen1 = [this_datum_circle[0]]
#                         ycen1 = [this_datum_circle[1]]
#                         rad1 = this_datum_circle[2]
#                         plotaxis = plotting.plot_circles( xcen1, ycen1,
#                             rad1, title="Datum actuators intersect", showplot=False )
#                         xcen2 = [other_datum_circle[0]]
#                         ycen2 = [other_datum_circle[1]]
#                         rad2 = other_datum_circle[2]
#                         plotaxis = plotting.plot_circles( xcen2, ycen2,
#                             rad2, showplot=True )
                    return True

                logger.debug("My datum rect: " + str(this_datum_rect))
                logger.debug("Other datum rect: " + str(other_datum_rect))
#                 if util.quadrangles_intersect(this_datum_rect, other_datum_rect):
                if util.polygons_intersect(this_datum_rect, other_datum_rect):
                    self.in_conflict = True
                    self.conflict_type = CONFLICT_DATUM
                    self.conflict_reason = "Conflict between %s and %s: " % \
                            (self.name, other.name)
                    self.conflict_reason += \
                            "Datum actuator rectangles intersect"
#                     if GRAPHICAL_DEBUGGING and plotting is not None:
#                         (xlist1, ylist1) = util.polygon_to_points(this_datum_rect)
#                         (xlist2, ylist2) = util.polygon_to_points(other_datum_rect)
#                         plotaxis = plotting.plot_xy(xlist1, ylist1,
#                                         title="Datum actuator rectangles conflict",
#                                         showplot=False)
#                         plotting.plot_xy(xlist2, ylist2, plotaxis=plotaxis,
#                                         showplot=True)
                    return True

            # Does the target lie within the other positioner's fibre avoidance triangle?
            if INCLUDE_FIBRE_SNAG:
                other_fibre_triangle = other.avoidance_triangle_fibre(padding=safety)
                logger.debug("Other fibre triangle: " + str(other_fibre_triangle))
                if util.point_inside_triangle(self.x_fibre_focal, self.y_fibre_focal,
                            other_fibre_triangle[0][0], other_fibre_triangle[0][1],
                            other_fibre_triangle[1][0], other_fibre_triangle[1][1],
                            other_fibre_triangle[2][0], other_fibre_triangle[2][1]):
                    self.in_conflict = True
                    self.conflict_type = CONFLICT_FIBRE_SNAG
                    self.conflict_reason = "Conflict between %s and %s: " % \
                        (self.name, other.name)
                    self.conflict_reason += "Target inside fibre snag triangle"
                    logger.debug(self.conflict_reason)
                    return True
                # Does the fibre holder intersect with the other positioner's avoidance
                # triangle (within the given safety tolerance)?
#                 if util.triangle_intersects_circle(other_fibre_triangle,
#                                                    self.x_fibre_focal,
#                                                    self.y_fibre_focal,
#                                                    (FIBRE_RADIUS+safety)):
                if util.polygon_intersects_circle(other_fibre_triangle,
                                                   self.x_fibre_focal,
                                                   self.y_fibre_focal,
                                                   (FIBRE_RADIUS+safety)):
                    self.in_conflict = True
                    self.conflict_type = CONFLICT_FIBRE_SNAG
                    self.conflict_reason = "Conflict between %s and %s: " % \
                        (self.name, other.name)
                    self.conflict_reason += "Fibre holder intersects fibre snag triangle"
                    logger.debug(self.conflict_reason)
#                     if GRAPHICAL_DEBUGGING and plotting is not None:
#                         (xlist, ylist) = util.polygon_to_points(other_fibre_triangle)
#                         plotaxis = plotting.plot_circles([self.x_fibre_focal],
#                             [self.y_fibre_focal], FIBRE_RADIUS,
#                             title="Fibre snag zone conflict", showplot=False)
#                         plotting.plot_xy(xlist, ylist, plotaxis=plotaxis, showplot=True)
#                         plotting.close()
                    return True

                if symmetrical:
                    # If a symmetrical check is needed, also test whether the other
                    # positioner's target lies within this positioner's fibre snag zone.
                    this_fibre_triangle = self.avoidance_triangle_fibre(padding=safety)
                    logger.debug("This fibre triangle: " + str(this_fibre_triangle))
                    if util.point_inside_triangle(other.x_fibre_focal, other.y_fibre_focal,
                            this_fibre_triangle[0][0], this_fibre_triangle[0][1],
                            this_fibre_triangle[1][0], this_fibre_triangle[1][1],
                            this_fibre_triangle[2][0], this_fibre_triangle[2][1]):
                        self.in_conflict = True
                        self.conflict_type = CONFLICT_FIBRE_SNAG
                        self.conflict_reason = "Conflict between %s and %s: " % \
                            (self.name, other.name)
                        self.conflict_reason += "Other's target inside fibre snag triangle"
                        logger.debug(self.conflict_reason)
                        return True
#                     if util.triangle_intersects_circle(this_fibre_triangle,
#                                                    other.x_fibre_focal,
#                                                    other.y_fibre_focal,
#                                                    (FIBRE_RADIUS+safety)):
                    if util.polygon_intersects_circle(this_fibre_triangle,
                                                   other.x_fibre_focal,
                                                   other.y_fibre_focal,
                                                   (FIBRE_RADIUS+safety)):
                        self.in_conflict = True
                        self.conflict_type = CONFLICT_FIBRE_SNAG
                        self.conflict_reason = "Conflict between %s and %s: " % \
                        (self.name, other.name)
                        self.conflict_reason += "Other's fibre holder inside fibre snag triangle"
                        logger.debug(self.conflict_reason)
#                         if GRAPHICAL_DEBUGGING and plotting is not None:
#                             (xlist, ylist) = util.polygon_to_points(this_fibre_triangle)
#                             plotaxis = plotting.plot_circles([other.x_fibre_focal],
#                                 [other.y_fibre_focal], FIBRE_RADIUS,
#                                 title="Fibre snag zone conflict in reverse", showplot=False)
#                             plotting.plot_xy(xlist, ylist, plotaxis=plotaxis, showplot=True)
#                             plotting.close()
                        return True

        elif isinstance(other, AcquisitionCamera):
            # Positioner/AcquisitionCamera conflict check

            # Is the acquisition camera active?
            if other.active:
                # Does the positioner overlap the field of view of the
                # acquisition camera?
                xdiff = self.x_fibre_focal - other.x_star_focal
                ydiff = self.y_fibre_focal - other.y_star_focal
                distsq = xdiff * xdiff + ydiff * ydiff
                mindist = other.fov_radius + FIBRE_RADIUS + safety
                mindistsq = mindist * mindist
                if distsq <= mindistsq:
                    self.in_conflict = True
                    self.conflict_type = CONFLICT_AC
                    self.conflict_reason = "Conflict between %s and %s: " % \
                        (self.name, other.name)
                    self.conflict_reason += "Positioner blocks AC: dist^2 %.3f <= %.3f" % \
                        (distsq, mindistsq)
                    # The acquisition camera is blocked
                    other.in_conflict = True
                    other.conflict_type = CONFLICT_AC
                    other.conflict_reason = "Blocked by Positioner %s." % self.name
                    return True

        elif isinstance(other, FiducialMarker):
            # Positioner/FiducialMarker conflict check

            # Is the fiducial marker active?
            if other.active:
                # Does the positioner overlap the avoidance zone of the
                # fiducial marker?
                xpoint = self.x_fibre_focal
                ypoint = self.y_fibre_focal
                xcen = other.avoidance_ellipse[0]
                ycen = other.avoidance_ellipse[1]
                semimajor = other.avoidance_ellipse[2]
                semiminor = other.avoidance_ellipse[3]
                tilt = other.avoidance_ellipse[4]
                if util.point_inside_ellipse(xpoint, ypoint, xcen, ycen,
                    semimajor, semiminor, tilt):
                    self.in_conflict = True
                    self.conflict_type = CONFLICT_FID
                    self.conflict_reason = "Conflict between %s and %s: " % \
                        (self.name, other.name)
                    self.conflict_reason += "Positioner blocks FID: (%.3f,%.3f) inside ellipse" % \
                        (xpoint, ypoint)
                    # The acquisition camera is blocked
                    other.in_conflict = True
                    other.conflict_type = CONFLICT_FID
                    other.conflict_reason = "Blocked by Positioner %s." % self.name
                    return True

        return False

    # ***
    def in_conflict_with_neighbours(self, test_all=True, safety=INS_POS_SAFETY,
                                    symmetrical=True):
        """

        Determine whether this positioner is in conflict with any of its
        neighbours, given their current targets.

        :Parameters:

        test_all: bool (optional)
            Set to True (the default) to test all the neighbours.
            If False, the function returns after the first conflict
            is detected.
        safety: float (optional)
            All conflict avoidance zones are padded with this safety
            tolerance. The default is the design tolerance, INS_POS_SAFETY.
        symmetrical: boolean (optional)
            Set to True to make all conflict checks symmetrical, so that when
            one positioner is in conflict both are in conflict. This only makes
            a difference to the fibre snag check.
            The default is True - when there is a fibre snag, both positioners
            are considered to be in conflict.

        :Attributes:

        The function sets the following attributes.

        in_conflict, conflict_type, conflict_reason
            Parameters which indicate if the positioner is in conflict.

        :Returns:

        in_conflict: boolean
            True if the positioner is in conflict with any of its neighbours
            False if the positioner is not in conflict with all neighbours.

        """
        in_conflict = False
        for neighbour in self.near_neighbours:
            logger.debug("Check for conflict between %s and near neighbour %s" % \
                          (self.name, neighbour.name))
            if self.in_conflict_with(neighbour, safety=safety,
                                     symmetrical=symmetrical):
                in_conflict = True
                if GRAPHICAL_DEBUGGING:
                    title = "%s in conflict with %s" % (self.name,
                                                        neighbour.name)
                    title += "\n" + self.conflict_reason
                    plotfig = self.plot(showplot=False, description=title)
                    neighbour.plot(plotfig=plotfig)
                if not test_all:
                    break
#             else:
#                 if GRAPHICAL_DEBUGGING:
#                     title = "%s NOT in conflict with %s" % (self.name,
#                                                             neighbour.name)
#                     plotfig = self.plot(showplot=False, description=title)
#                     neighbour.plot(plotfig=plotfig)

        for neighbour in self.far_neighbours:
            logger.debug("Check for conflict between %s and far neighbour %s" % \
                          (self.name, neighbour.name))
            # A far neighbour can only conflict when the positioner's
            # reach is outside a defined limit.
            if self.r_fibre_local >= (self.limitdist - 2.0*safety):
                if self.in_conflict_with(neighbour, safety=safety,
                                         symmetrical=symmetrical):
                    in_conflict = True
                    if GRAPHICAL_DEBUGGING:
                        title = "%s in conflict with %s" % (self.name,
                                                            neighbour.name)
                        title += "\n" + self.conflict_reason
                        plotfig = self.plot(showplot=False, description=title)
                        neighbour.plot(plotfig=plotfig)
                    if not test_all:
                        break
#                 else:
#                     if GRAPHICAL_DEBUGGING:
#                         title = "%s NOT in conflict with %s" % (self.name,
#                                                                 neighbour.name)
#                         plotfig = self.plot(showplot=False, description=title)
#                         neighbour.plot(plotfig=plotfig)

        for camera in self.near_cameras:
            logger.debug("Check for conflict between %s and acquisition camera %s" % \
                          (self.name, camera.name))
            if self.in_conflict_with(camera, safety=safety):
                in_conflict = True
                if GRAPHICAL_DEBUGGING:
                    title = "%s in conflict with %s" % (self.name,
                                                        camera.name)
                    title += "\n" + self.conflict_reason
                    plotfig = self.plot(showplot=False, description=title)
                    camera.plot(plotfig=plotfig)
                if not test_all:
                    break
#             else:
#                 if GRAPHICAL_DEBUGGING:
#                     title = "%s NOT in conflict with %s" % (self.name,
#                                                             camera.name)
#                     plotfig = self.plot(showplot=False, description=title)
#                     camera.plot(plotfig=plotfig)

        for fiducial in self.near_fiducials:
            logger.debug("Check for conflict between %s and fiducial marker %s" % \
                          (self.name, fiducial.name))
            if self.in_conflict_with(fiducial, safety=safety):
                in_conflict = True
                if GRAPHICAL_DEBUGGING:
                    title = "%s in conflict with %s" % (self.name,
                                                        fiducial.name)
                    title += "\n" + self.conflict_reason
                    plotfig = self.plot(showplot=False, description=title)
                    fiducial.plot(plotfig=plotfig)
                if not test_all:
                    break
#             else:
#                 if GRAPHICAL_DEBUGGING:
#                     title = "%s NOT in conflict with %s" % (self.name,
#                                                             fiducial.name)
#                     plotfig = self.plot(showplot=False, description=title)
#                     fiducial.plot(plotfig=plotfig)

        return in_conflict

    def plot(self, description='', plotfig=None, showplot=True):

        """

        Plot the configuration of the positioner.

        NOTE: This function requires access to the plotting
        module, plotting.py. If the module can't be imported
        the function returns with an apologetic message.

        :Parameters:

        description: str, optional
            Optional description to be added to the positioner plot.
        plotfig: matplotlib Figure object, optional
            An existing matplotlib figure onto which to add the plot.
            If not given, a new figure is created.
        showplot: boolean, optional
            If True (the default) show the plot when finished.
            Set of False if you are overlaying several plots, and
            this plot is not the last.

        :Returns:

        plotfig: matplotlib Figure object
            A matplotlib figure containing the plot.

        """
        if plotting is None:
            logger.warning("Plotting is disabled")
            return

        # Create a new plotting figure, if necessary.
        if plotfig is None:
            plotfig = plotting.new_figure(1, figsize=(10,9), stitle='')

        # First plot the centres of the positioners on the grid
        # using "O" symbols.
        xcentres = [self.x_centre_focal]
        ycentres = [self.y_centre_focal]

        # Plot the centre circles showing the patrol field.
        # Inner boundary - black dashed line.
        plotaxis = plotting.plot_circles( xcentres, ycentres,
                        self.inner, plotfig=plotfig,
                        facecolor='w', circcolor='k ', linefmt='ko',
                        linestyle='--', linewidth=0.3,
                        xlabel="X (mm)", ylabel="Y (mm)", title=description,
                        grid=True, showplot=False )
        # Nominally "Safe" boundary - blue dashed line.
        plotaxis = plotting.plot_circles( xcentres, ycentres,
                        self.outer/2.0,
                        plotfig=plotfig, plotaxis=plotaxis,
                        facecolor='w', circcolor='b ', linefmt='b ',
                        linestyle='--', linewidth=0.2, grid=False,
                        showplot=False )
        # Outer boundary - black dashed line.
        plotaxis = plotting.plot_circles( xcentres, ycentres,
                        self.outer, plotfig=plotfig, plotaxis=plotaxis,
                        facecolor='w', circcolor='k ', linefmt='k ',
                        linestyle='--', linewidth=0.3, grid=False,
                        showplot=False )

        # Add a label giving the positioner ID
        xtxt = self.x_centre_focal - self.length1/3.0
        ytxt = self.y_centre_focal - self.length1/3.0
        plotaxis = plotting.plot_text( str(self.ident), xtxt, ytxt,
                                       plotfig=plotfig, plotaxis=plotaxis,
                                       showplot=False)

        # Change the colour of the positioner arms depending on whether there
        # is a target assigned and whether the positioner is in conflict.
        if self.target_assigned:
            target_fmt = 'b+'   # Blue + - assigned
            if self.in_conflict:
                arm_fmt = 'r '  # Red - in conflict
            else:
                arm_fmt = 'k '  # Black - target assigned
        else:
            target_fmt = 'gx'   # Green x - not assigned
            if self.in_conflict:
                arm_fmt = 'r '  # Red - in conflict
            elif self.conflict_type == CONFLICT_UNREACHABLE:
                arm_fmt = 'm '  # Magenta - attempted to assign unreachable target
            else:
                arm_fmt = 'g '  # Green - target not assigned

        # Add to the plot the location of all the fibres using the
        # symbol defined in target_fmt.
        xtargets = [self.x_fibre_focal]
        ytargets = [self.y_fibre_focal]
        plotaxis = plotting.plot_xy( xtargets, ytargets, plotfig=plotfig,
                                     plotaxis=plotaxis,
                                     linefmt=target_fmt, linestyle=' ',
                                     showplot=False )
        # Show the fibre approach tolerance as a thick red circle.
        plotaxis = plotting.plot_circles( xtargets, ytargets,
                        FIBRE_RADIUS,
                        plotfig=plotfig, plotaxis=plotaxis,
                        facecolor='w', circcolor='r ', linefmt='r ',
                        linestyle='-', linewidth=3, grid=False,
                        showplot=False )
        if PLOT_PADDING:
            plotaxis = plotting.plot_circles( xtargets, ytargets,
                            FIBRE_RADIUS+INS_POS_SAFETY,
                            plotfig=plotfig, plotaxis=plotaxis,
                            facecolor='w', circcolor='r ', linefmt='r ',
                            linestyle='--', linewidth=3, grid=False,
                            showplot=False )

        # Show the orientation of the alpha and beta arms as thick lines.
        x_elbow_focal = self.x_centre_focal + self.x_elbow_local
        y_elbow_focal = self.y_centre_focal + self.y_elbow_local
        xlist = [self.x_centre_focal, x_elbow_focal, self.x_fibre_focal]
        ylist = [self.y_centre_focal, y_elbow_focal, self.y_fibre_focal]
        plotaxis = plotting.plot_xy( xlist, ylist, plotfig=plotfig,
                    plotaxis=plotaxis, linefmt=arm_fmt, linestyle='-',
                    linewidth=6, showplot=False )

        # Add thick lined polygons to show the avoidance zones. Note that
        # avoidance zones are plotted without the safety tolerance.
        # (1a) The metrology target zone - red.
        poly1 = self.avoidance_polygon_metrology()
        (xlist, ylist) = util.polygon_to_points(poly1)
        plotaxis = plotting.plot_xy( xlist, ylist, plotfig=plotfig,
                    plotaxis=plotaxis, linefmt='r ', linestyle='-',
                    linewidth=3, showplot=False )
        if PLOT_PADDING:
            poly1_padded = \
                self.avoidance_polygon_metrology(padding=INS_POS_SAFETY)
            (xlist, ylist) = util.polygon_to_points(poly1_padded)
            plotaxis = plotting.plot_xy( xlist, ylist, plotfig=plotfig,
                            plotaxis=plotaxis, linefmt='r ', linestyle='--',
                            linewidth=3, showplot=False )

        # (1b) Cover the metrology zone with an ellipse.
        # FIXME: The ellipse becomes squashed at certain locations and tilts.
        if PLOT_ELLIPSES:
            (xcen, ycen, major, minor, tilt) = \
                self.avoidance_ellipse_metrology()
            plotaxis = plotting.plot_ellipses( [xcen], [ycen], major, minor,
                            tilt, plotfig=plotfig, plotaxis=plotaxis,
                            facecolor='w', ellipsecolor='k ', linefmt='r ',
                            linestyle='-', linewidth=3, grid=False,
                            showplot=False )

        # (2) The fibre snag avoidance zone - yellow.
        if INCLUDE_FIBRE_SNAG:
            triangle2 = self.avoidance_triangle_fibre()
            (xlist, ylist) = util.polygon_to_points(triangle2)
            plotaxis = plotting.plot_xy( xlist, ylist, plotfig=plotfig,
                        plotaxis=plotaxis, linefmt='y ', linestyle='-',
                        linewidth=3, showplot=False )
#         if PLOT_PADDING:
#             triangle2_padded = \
#                 self.avoidance_triangle_fibre(padding=INS_POS_SAFETY)
#             (xlist, ylist) = util.polygon_to_points(triangle2_padded)
#             plotaxis = plotting.plot_xy( xlist, ylist, plotfig=plotfig,
#                             plotaxis=plotaxis, linefmt='y ', linestyle='--',
#                             linewidth=3, showplot=False )

        # (3a) Add a thick dotted rectangles showing the beta arm
        # collision zones - magenta.
        if INCLUDE_BETA_ZONE:
            rectangles_beta = self.avoidance_rectangles_beta()
            for rect in rectangles_beta:

                (xlist, ylist) = util.polygon_to_points(rect)
                plotaxis = plotting.plot_xy( xlist, ylist, plotfig=plotfig,
                                plotaxis=plotaxis, linefmt='m ', linestyle='-.',
                                linewidth=3, showplot=False )
            if PLOT_PADDING:
                rectangles_beta_padded = \
                    self.avoidance_rectangles_beta(padding=INS_POS_SAFETY)
                for rect in rectangles_beta_padded:

                    (xlist, ylist) = util.polygon_to_points(rect)
                    plotaxis = plotting.plot_xy( xlist, ylist, plotfig=plotfig,
                                    plotaxis=plotaxis, linefmt='m ',
                                    linestyle=':', linewidth=3,
                                    showplot=False )

        # (3b) Cover the beta zone with an ellipse.
        # FIXME: The ellipse becomes squashed at certain locations and tilts.
        if PLOT_ELLIPSES:
            (xcen, ycen, major, minor, tilt) = \
                self.avoidance_ellipse_beta()
            plotaxis = plotting.plot_ellipses( [xcen], [ycen], major, minor,
                            tilt, plotfig=plotfig, plotaxis=plotaxis,
                            facecolor='w', ellipsecolor='b ', linefmt='r ',
                            linestyle='-', linewidth=3, grid=False,
                            showplot=False )

        # (4) Add a thick dotted circle showing the location of the
        # "datum foot" of the beta arm - cyan.
        if INCLUDE_DATUM_FOOT:
            (rect4, circle4) = self.avoidance_zone_datum()
            plotaxis = plotting.plot_circles( [circle4[0]], [circle4[1]],
                            circle4[2], plotfig=plotfig, plotaxis=plotaxis,
                            facecolor='w', circcolor='c ', linefmt='co',
                            linestyle=':', linewidth=3,
                            xlabel="X (mm)", ylabel="Y (mm)", title=description,
                            grid=True, showplot=False )
            (xlist, ylist) = util.polygon_to_points(rect4)
            plotaxis = plotting.plot_xy( xlist, ylist, plotfig=plotfig,
                            plotaxis=plotaxis, linefmt='c ', linestyle=':',
                            linewidth=3, showplot=False )
            if PLOT_PADDING:
                (rect4_pad, circle4_pad) = \
                    self.avoidance_zone_datum(padding=INS_POS_SAFETY)
                plotaxis = plotting.plot_circles( [circle4_pad[0]],
                            [circle4_pad[1]], circle4_pad[2], plotfig=plotfig,
                            plotaxis=plotaxis,
                            facecolor='w', circcolor='c ', linefmt='co',
                            linestyle=':', linewidth=2,
                            xlabel="X (mm)", ylabel="Y (mm)", title=description,
                            grid=True, showplot=False )
                (xlist, ylist) = util.polygon_to_points(rect4_pad)
                plotaxis = plotting.plot_xy( xlist, ylist, plotfig=plotfig,
                            plotaxis=plotaxis, linefmt='c ', linestyle=':',
                            linewidth=2, showplot=False )
        if showplot:
            plotting.show_plot()
            plotting.close()
        return plotfig

    def plot_status(self, description='', plotfig=None, showplot=True):
        """

        Plot the status of the positioner.

        NOTE: This function requires access to the plotting
        module, plotting.py. If the module can't be imported
        the function returns with an apologetic message.

        :Parameters:

        description: str, optional
            Optional description to be added to the positioner plot.
        plotfig: matplotlib Figure object, optional
            An existing matplotlib figure onto which to add the plot.
            If not given, a new figure is created.
        showplot: boolean, optional
            If True (the default) show the plot when finished.
            Set of False if you are overlaying several plots, and
            this plot is not the last.

        :Returns:

        plotfig: matplotlib Figure object
            A matplotlib figure containing the plot.

        """
        if plotting is None:
            logger.warning("Plotting is disabled")
            return

        # Create a new plotting figure, if necessary.
        if plotfig is None:
            plotfig = plotting.new_figure(1, figsize=(10,9), stitle='')

        # Change the colour of the filled circle, depending on the status of
        # the fibre positioner. The ISFS suggests the following colour scheme.
        # Flashing white - Fibre positioner is moving
        # Green - Fibre positioner in position on science target.
        # Cyan - Fibre positioner in position on blank sky.
        # Magenta - Positioner at default location. (Target not assigned or out of reach)
        # Yellow - Fibre positioner in deadlock (no collision free path).
        # Red - Collision or conflict detected.
        # Blue - Cell occupied by acquisition camera
        # Black - Fibre positioner broken or unavailable.
        if self.simulated or self.locked:
                facecolor = 'black'
                circcolor = 'k '
                linefmt = 'k '
        elif self.target_assigned:
            if self.in_conflict:
                facecolor = 'red'
                circcolor = 'r '
                linefmt = 'r '
            else:
                facecolor = 'green'
                circcolor = 'g '
                linefmt = 'g '
        else:
            if self.in_conflict:
                facecolor = 'red'
                circcolor = 'r '
                linefmt = 'r '
            elif self.conflict_type == CONFLICT_UNREACHABLE:
                facecolor = 'magenta'
                circcolor = 'm '
                linefmt = 'm '
            else:
                facecolor = 'cyan'
                circcolor = 'c '
                linefmt = 'c '

        # Plot fibre positioners as filled circles.
        xcentres = [self.x_centre_focal]
        ycentres = [self.y_centre_focal]
        plotaxis = plotting.plot_circles( xcentres, ycentres,
                        (self.length1+self.length2)/2.0, plotfig=plotfig,
                        facecolor=facecolor, circcolor=circcolor, linefmt=linefmt,
                        linestyle='-', linewidth=0.3, grid=False, filled=True,
                        showplot=False )

        if showplot:
            plotting.show_plot()
            plotting.close()
        return plotfig

    # ADDITIONAL METHODS FROM PATH ANALYSIS RELEASE - 29 JUL 2016
    def motor_position_to_focal_cartesian(self, angle1, angle2):
        """

        Function used by the path analysis software to calculate the
        (X,Y) coordinate of the the fibre when the arm motors are
        moved to the given angles. It is used for determining a path
        without defining the target coordinates.

        The simulated positioner is NOT moved to these angles.

        c.f. the set_arm_angles function makes the same calculation but
        moves the simulated positioner to the defined angles.

        :Parameters:

        angle1: float
            Alpha motor angle, with respect to its datum location (in radians)
        angle2: float
            Beta motor angle, with respect to its datum location (in radians)

        :Attributes:

        The function sets NO attributes.

        :Returns:

        (x_fibre_focal, y_fibre_focal)
            Focal plane Cartesian coordinates of the fibre

        """
        # TODO: THIS FUNCTION OVERLAPS WITH SET_ARM_ANGLES. FIND A COMMON SOLUTION.
        x_fibre_focal = self.x_centre_focal + self.length1 * math.cos(angle1) + \
                        self.length2 * math.cos(angle1 + angle2 + math.pi)
        y_fibre_focal = self.y_centre_focal + self.length1 * math.sin(angle1) + \
                        self.length2 * math.sin(angle1 + angle2 + math.pi)

        return (x_fibre_focal, y_fibre_focal)

    # Add motors
    # TODO: REWORK. Read the maximum and minimum speed from an interface with hardware
    # max_speed_alpha = 2.0*math.pi*ALPHA_RPM_LIMITS[1]/60.0
    # min_speed_alpha = 2.0*math.pi*ALPHA_RPM_LIMITS[0]/60.0
    # max_speed_beta = 2.0*math.pi*BETA_RPM_LIMITS[1]/60.0
    # min_speed_beta = 2.0*math.pi*BETA_RPM_LIMITS[0]/60.0
    def add_motors(self, sim_length,
                   alpha_position_target, beta_position_target,
                   max_speed_alpha, min_speed_alpha,
                   max_speed_beta, min_speed_beta):
        """

        Function used by the path analysis software to associate two
        Motor objects with a Positioner object.

        :Parameters:

        sim_length: int
            Number of elements in motor path array.
        alpha_position_target: float
            Target alpha motor angle (in radians)
        beta_position_target: float
            Target alpha motor angle (in radians)
        max_speed_alpha: float
            The maximum allowed alpha motor speed (in radians per second)
        min_speed_alpha: float
            The minimum allowed alpha motor speed (in radians per second)
        max_speed_beta: float
            The maximum allowed alpha motor speed (in radians per second)
        min_speed_beta: float
            The minimum allowed alpha motor speed (in radians per second)

        :Attributes:

        motor1: Motor
            Alpha motor object.
        motor2: Motor
            Beta motor object.

        :Returns:

        None

        """
        self.motor1 = Motor( sim_length, alpha_position_target,
                             max_speed_alpha, min_speed_alpha )
        self.motor2 = Motor( sim_length, beta_position_target,
                             max_speed_beta, min_speed_beta )


class PositionerGrid(object):
    """

    Class representing a hexagonal grid of fibre positioners.

    Each grid is built from a collection of fibre positioners at
    different MOONS focal plane locations and orientations. Column and
    row numbers in the grid are also given explicitly to disambiguate
    the cells.

    The defined locations and orientations are used to acquire targets
    and check for conflicts between fibre positioners.

    The defined column and row numbers are used to determine which fibre
    positioners are neighbours.

    """
    # Class variables defining the location of information within the
    # fibre positioner configuration list.
    C_IDENT = 0
    C_RCEN = 1
    C_THCEN = 2
    C_ORIENT = 3
    C_COLUMN = 4
    C_ROW = 5
    C_SIMULATED = 6
    C_LOCKED = 7
    C_RFIBRE = 8
    C_THFIBRE = 9
    C_PFIBRE = 10

    # Class variables defining the location of information within the
    # acquisition camera configuration list.
    A_IDENT = 0
    A_RCEN = 1
    A_THCEN = 2
    A_ORIENT = 3
    A_FOV = 4
    A_COLUMN = 5
    A_ROW = 6

    # Class variables defining the location of information within the
    # fiducial marker configuration list.
    F_IDENT = 0
    F_XCEN = 1
    F_YCEN = 2
    F_ORIENT = 3
    F_MINOR = 4
    F_MAJOR = 5
    F_COLUMN = 6
    F_ROW = 7

    def __init__(self, config_list, camera_list=[], fiducial_list=[],
                 length1=INS_POS_LENGTH1, length2=INS_POS_LENGTH2):
        """

        Constructor for PositionerGrid class.

        :Parameters:

        config_list: list of (ident, rcen, thcen, orient, column, row,
            [simulated, locked, rfibre, thfibre, pfibre])
            List of configuration parameters describing the IDs,
            locations, orientations, column numbers and row numbers
            of the fibre positioners. Each list may also contain
            a simulation flag, positioner locked (broken) flag and
            a default location for the fibre.
        camera_list: list of (ident, rcen, thcen, orient, fov, column, row) (optional)
            List of configuration parameters describing the IDs,
            locations, orientations, field of view radii, column numbers and
            row numbers of the acquisition cameras. Defaults to empty list.
        fiducial_list: list of (ident, xcen, ycen, orient, sminor, smajor, [column, row]) (optional)
            List of configuration parameters describing the IDs, locations,
            orientations, semi-minor and semi-major axes. column numbers and
            row numbers of the fiducial markers. Defaults to empty list.
        length1: float (optional)
            Length of all alpha arms.
            Defaults to the designed length, INS_POS_LENGTH1
        length2: float (optional)
            Length of all beta arms.
            Defaults to the designed length, INS_POS_LENGTH2

        """
        assert isinstance(config_list, (tuple,list))
        assert isinstance(config_list[0], (tuple,list))
        assert len(config_list) > 1
        assert len(config_list[0]) >= 6

        # A list in which to store the positioners and a dictionary
        # in which to translate ident number into row and column.
        self.grid = []
        self.positioner_count = 0
        self.ident_dict = {}

        # The pitch size of the grid is assumed to be the same as the reach of
        # the positioners.
        self.pitch = length1 + length2

        # Determine the number of columns and rows in the grid.
        cmin = config_list[0][self.C_COLUMN]
        rmin = config_list[0][self.C_ROW]
        cmax = config_list[0][self.C_COLUMN]
        rmax = config_list[0][self.C_ROW]
        idmax = config_list[0][self.C_IDENT]
        for config in config_list:
            if config[self.C_COLUMN] < cmin:
                cmin = config[self.C_COLUMN]
            if config[self.C_ROW] < rmin:
                rmin = config[self.C_ROW]
            if config[self.C_COLUMN] > cmax:
                cmax = config[self.C_COLUMN]
            if config[self.C_ROW] > rmax:
                rmax = config[self.C_ROW]
            if config[self.C_IDENT] > idmax:
                idmax = config[self.C_IDENT]
        # Use the range of column and row numbers given
        self.max_columns = cmax + 1
        self.max_rows = rmax + 1
        self.columns = cmax - cmin + 1
        self.rows = rmax - rmin + 1
        self.idmax = idmax
        pmax = idmax + 1
        strg = "Fibre positioner grid contains %d grid locations, " % idmax
        strg += "columns %d-%d and rows %d-%d (or %d x %d)." % \
            (cmin, cmax, rmin, rmax, self.columns, self.rows)
        logger.info(strg)

        # Create an empty list in which to store the positioners
        self.positioners = [None] * pmax

        # Create an empty grid of the correct size.
        # TODO: Creating a rectangular grid to cover a circular focal plane
        # leaves several blank locations. Is there a better way of addressing
        # positioners but still allowing neighbours to be quickly identified?
        for row in range(0, self.max_rows):
            grid_rows = []
            for column in range(0, self.max_columns):
                grid_rows.append(None)
            self.grid.append(grid_rows)

        # Populate the grid with positioners. column and row are assumed
        # to have been already derived correctly from the centre of each
        # positioner.
        for config in config_list:
            ident = config[self.C_IDENT]
            rcen = config[self.C_RCEN]
            thcen = config[self.C_THCEN]
            orient = config[self.C_ORIENT]
            column = config[self.C_COLUMN]
            row = config[self.C_ROW]
            if self.positioners[ident] is None:
                if self.grid[row][column] is None:
                    # TODO: Locked (broken) positioners not taken into account.
                    positioner = Positioner(ident, rcen, thcen, orient,
                                            column=column, row=row,
                                            simulated=False, locked=False,
                                            length1=length1, length2=length2)
                    self.positioners[ident] = positioner
#                     print("Adding positioner %d at (%d,%d)" % (ident, row, column))
                    self.grid[row][column] = ident
                    self.ident_dict[str(ident)] = (column,row) # NEEDED?
                    self.positioner_count += 1
                    logger.debug("Positioner %s added at row=%d, column=%d" % \
                             (positioner.name, row, column) )
                else:
                    # Two or more positioners defined at the same column and row
                    strg = "More than one positioner found at row=%d, column=%d\n" % \
                        (row, column)
                    strg += "\tAttempt to add positioner "
                    strg += "ident=%d, rcen=%f, thcen=%f\n" % (ident, rcen, thcen)
                    strg += "\tSpace occupied by "
                    ident = self.grid[row][column]
                    other = self.positioners[ident]
                    strg += "ident=%d, rcen=%f, thcen=%f" % \
                        (other.ident, other.r_centre_focal, other.theta_centre_focal)
                    raise AttributeError(strg)
            else:
                # Two or more positioners defined with the same identifier
                strg = "More than one positioner with identifier %d\n" % ident
                strg += "\tAttempt to add positioner "
                strg += "rcen=%f, thcen=%f\n" % (rcen, thcen)
                strg += "\tSpace occupied by "
                other = self.positioners[ident]
                strg += "rcen=%f, thcen=%f" % \
                    (other.r_centre_focal, other.theta_centre_focal)
                raise AttributeError(strg)

        # Create a list of acquisition cameras. If the column and row are
        # provided, check that no positioner is occupying the same cell
        # in the grid.
        self.acquisition_cameras = []
        for config in camera_list:
            ident = config[self.A_IDENT]
            rcen = config[self.A_RCEN]
            thcen = config[self.A_THCEN]
            orient = config[self.A_ORIENT]
            fov_radius = config[self.A_FOV]
            if len(config) > self.A_ROW:
                column = config[self.A_COLUMN]
                row = config[self.A_ROW]
            else:
                column = None
                row = None
            if row is not None and column is not None:
                # Only test acquisition cameras located inside the bounds of
                # the fibre positioner grid.
                if row >= rmin and row <= rmax and \
                   column >= cmin and column <= cmax and \
                   self.grid[row][column] is not None:
                    # Two or more objects defined at the same column and row
                    strg = "More than one object found at row=%d, column=%d\n" % \
                        (row, column)
                    strg += "\tAttempt to add acquisition camera "
                    strg += "ident=%d, rcen=%f, thcen=%f\n" % (ident, rcen, thcen)
                    strg += "\tSpace occupied by positioner "
                    ident = self.grid[row][column]
                    other = self.positioners[ident]
                    strg += "ident=%d, rcen=%f, thcen=%f" % \
                        (other.ident, other.r_centre_focal, other.theta_centre_focal)
                    raise AttributeError(strg)

            ac = AcquisitionCamera( ident, rcen, thcen, orient, fov_radius,
                                    column=column, row=row )
            self.acquisition_cameras.append(ac)

        # Create a list of fiducial markers. If the column and row are
        # provided, check that no positioner is occupying the same cell
        # in the grid.
        self.fiducials = []
        for fconfig in fiducial_list:
            ident = fconfig[self.F_IDENT]
            xcen = fconfig[self.F_XCEN]
            ycen = fconfig[self.F_YCEN]
            orient = fconfig[self.F_ORIENT]
            sminor = config[self.F_MINOR]
            smajor = config[self.F_MAJOR]
            if len(fconfig) > self.F_ROW:
                column = fconfig[self.F_COLUMN]
                row = fconfig[self.F_ROW]
            else:
                column = None
                row = None
            if row is not None and column is not None:
                # Only test fiducial markers located inside the bounds of
                # the fibre positioner grid.
                if row >= rmin and row <= rmax and \
                   column >= cmin and column <= cmax and \
                   self.grid[row][column] is not None:
                    # Two or more objects defined at the same column and row
                    strg = "More than one object found at row=%d, column=%d\n" % \
                        (row, column)
                    strg += "\tAttempt to add fiducial marker "
                    strg += "ident=%d, xcen=%f, ycen=%f\n" % (ident, xcen, ycen)
                    strg += "\tSpace occupied by positioner "
                    ident = self.grid[row][column]
                    other = self.positioners[ident]
                    strg += "ident=%d, xcen=%f, ycen=%f" % \
                        (other.ident, other.x_centre_focal, other.y_centre_focal)
                    raise AttributeError(strg)

            fid = FiducialMarker( ident, xcen, ycen, orient, sminor, smajor,
                                  column=column, row=row )
            self.fiducials.append(fid)

        # Set up the references between each positioner and its neighbours.
        self._define_neighbours()

        # Initialise the conflict counters.
        self.init_counters()

    def __del__(self):
        """

        Destructor

        """
        try:
            for positioner in self.positioners:
                if positioner is not None:
                    del positioner
            del self.positioners
            for camera in self.acquisition_cameras:
                if camera is not None:
                    del camera
            del self.acquisition_cameras
            for fiducial in self.fiducials:
                if fiducial is not None:
                    del fiducial
            del self.fiducials
        except:
            pass

    def __str__(self):
        """

        Return a readable string describing the positioner grid.

        """
        strg = '=' * 75 + '\n'
        strg += "Grid of %d positioners, %d ACs and %d fiducials" % \
            (self.positioner_count, len(self.acquisition_cameras), len(self.fiducials))
        strg += " in %d cols x %d rows with pitch=%.3f:\n" % \
            (self.columns, self.rows, self.pitch)
        for positioner in self.positioners:
            if positioner is not None:
                strg += '-' * 75 + '\n'
                strg += str(positioner) + '\n'

        count = 0
        for ac in self.acquisition_cameras:
            if count < 1:
                strg += '=' * 75 + '\n'
            else:
                strg += '-' * 75 + '\n'
            strg += str(ac) + '\n'
            count += 1

        count = 0
        for fid in self.fiducials:
            if count < 1:
                strg += '=' * 75 + '\n'
            else:
                strg += '-' * 75 + '\n'
            strg += str(fid) + '\n'
            count += 1

        return strg

    def _define_neighbours(self):
        """

        Helper function which defines the neighbours surrounding
        each positioner in the grid. It also associates positioners
        with fiducials which may affect their patrol zone.
        Called on initialisation.

        :Parameters:

        None.

        """
        # Near neighbours are within an inner ring adjacent to each
        # positioner. Far neighbours are within the second ring adjacent
        # to the near neighbours. Near and far neighbours both overlap
        # each positioner's patrol field, but far neighbours only need to
        # be considered when a positioner is reaching beyond a certain limit.
        for row in range(0, self.max_rows):
            # The column offsets are slightly different for odd and even
            # numbered rows in the grid.
            if row % 2 == 0:
                # Even numbered row
                for column in range(0, self.max_columns):
                    positioner = self.get_row_column(row, column)
                    if positioner is not None:
                        neighbour_offsets = [ \
                            (row-2, column-1, False),
                            (row-2, column, False),
                            (row-2, column+1, False),
                            (row-1, column-1, False),
                            (row-1, column, True),
                            (row-1, column+1, True),
                            (row-1, column+2, False),
                            (row, column-2, False),
                            (row, column-1, True),
                            (row, column+1, True),
                            (row, column+2, False),
                            (row+1, column-1, False),
                            (row+1, column, True),
                            (row+1, column+1, True),
                            (row+1, column+2, False),
                            (row+2, column-1, False),
                            (row+2, column, False),
                            (row+2, column+1, False),
                            ]
                        # Add neighbours at the specified row and column offsets.
                        for (roff, coff, near) in neighbour_offsets:
                            neighbour = self.get_row_column(roff, coff)
                            if neighbour is not None:
                                positioner.add_neighbour(neighbour, near=near)

                        if logger.getEffectiveLevel() == logging.DEBUG:
                            strg = "Even row. Positioner %s has the following neighbours." % \
                                positioner.name
                            strg += "\nNEAR neighbours: "
                            for neighbour in positioner.near_neighbours:
                                strg += "%s " % neighbour.name
                            strg += "\nFAR neighbours: "
                            for neighbour in positioner.far_neighbours:
                                strg += "%s " % neighbour.name
                            logger.debug(strg)
            else:
                # Odd numbered row
                for column in range(0, self.max_columns):
                    positioner = self.get_row_column(row, column)
                    if positioner is not None:
                        neighbour_offsets = [ \
                            (row-2, column-1, False),
                            (row-2, column, False),
                            (row-2, column+1, False),
                            (row-1, column-2, False),
                            (row-1, column-1, True),
                            (row-1, column, True),
                            (row-1, column+1, False),
                            (row, column-2, False),
                            (row, column-1, True),
                            (row, column+1, True),
                            (row, column+2, False),
                            (row+1, column-2, False),
                            (row+1, column-1, True),
                            (row+1, column, True),
                            (row+1, column+1, False),
                            (row+2, column-1, False),
                            (row+2, column, False),
                            (row+2, column+1, False),
                            ]
                        # Add neighbours at the specified row and column offsets.
                        for (roff, coff, near) in neighbour_offsets:
                            neighbour = self.get_row_column(roff, coff)
                            if neighbour is not None:
                                positioner.add_neighbour(neighbour, near=near)
                        if logger.getEffectiveLevel() == logging.DEBUG:
                            strg = "Odd row. Positioner %s has the following neighbours." % \
                                positioner.name
                            strg += "\nNEAR neighbours: "
                            for neighbour in positioner.near_neighbours:
                                strg += "%s " % neighbour.name
                            strg += "\nFAR neighbours: "
                            for neighbour in positioner.far_neighbours:
                                strg += "%s " % neighbour.name
                            logger.debug(strg)

        # Now link each positioner with acquisition cameras and fiducials
        # likely to interfere with its patrol zone.
        for positioner in self.positioners:
            if positioner is not None and isinstance(positioner, Positioner):

                for camera in self.acquisition_cameras:
                    xdiff = (positioner.x_centre_focal - camera.x_centre_focal)
                    ydiff = (positioner.y_centre_focal - camera.y_centre_focal)
                    dist = math.sqrt(xdiff * xdiff + ydiff * ydiff)
                    limitdist = positioner.outer + camera.fov_radius
                    if dist <= limitdist:
                        positioner.add_acquisition_camera(camera)

                for fiducial in self.fiducials:
                    xdiff = (positioner.x_centre_focal - fiducial.x_centre_focal)
                    ydiff = (positioner.y_centre_focal - fiducial.y_centre_focal)
                    dist = math.sqrt(xdiff * xdiff + ydiff * ydiff)
                    limitdist = positioner.outer + fiducial.avoid_major
                    if dist <= limitdist:
                        positioner.add_fiducial(fiducial)

    def get_positioner(self, ident):
        """

        Return the positioner in the grid with the given ID.

        :Parameters:

        ident: int
            Positioner ID.

        :Returns:

        positioner: Positioner
            The positioner object with the given ID, or None if the positioner
            cannot be found.

        """
        if ident >= 0 and ident < len(self.positioners):
            return self.positioners[ident]
        else:
            return None

    def get_row_column(self, row, column):
        """

        Return the positioner in the grid at the given row and column location.

        :Parameters:

        row: int
            Row number in the grid.
        column: int
            Column number in the grid.

        :Returns:

        positioner: Positioner
            The positioner object with the given location, or None if the
            positioner cannot be found.

        """
        if row >= 0 and row < self.max_rows and column >= 0 and column < self.max_columns:
            ident = self.grid[row][column]
            if ident is not None:
                return self.positioners[ident]
            else:
                return None
        else:
            return None

    def define_guide_stars(self, ac_stars):
        """

        Activates the given list of acquisition cameras and define the
        locations of the guide stars.

        :Parameters:

        ac_stars: list of (ID, [r_star, theta_star, fov_star])
            or list of ID.
            The IDs of the acquisition cameras to be activated and
            (optionally) the location and surrounding avoidance fov
            of the guide star. If a guide star coordinate is not
            given, the centre of the acquisition camera is used.

            Acquisition cameras not included in the list are deactivated.

        :Returns:

        None

        """
        assert isinstance(ac_stars, (tuple,list))

        # Start by deactivating all acquisition cameras
        for camera in self.acquisition_cameras:
            if camera is not None:
                camera.initialise()

        if isinstance(ac_stars[0], (tuple,list)):
            # A list of acquisition camera parameters has been provided
            # Activate each camera and define its guide stars.
            for camera in self.acquisition_cameras:
                if camera is not None:
                    for tobeactivated in ac_stars:
                        ident = tobeactivated[0]
                        if camera.ident == ident:
                            if len(tobeactivated) > 2:
                                r_star = tobeactivated[1]
                                theta_star = tobeactivated[2]
                            else:
                                r_star = None
                                theta_star = None
                            if len(tobeactivated) > 3:
                                fov_star = tobeactivated[3]
                            else:
                                fov_star = None
                            camera.set_star(r_star_focal=r_star,
                                            theta_star_focal=theta_star,
                                            fov_radius=fov_star)
                            break
        else:
            # A list of acquisition camera IDs has been provided (no stars).
            # Activate each camera.
            for camera in self.acquisition_cameras:
                if camera is not None:
                    for tobeactivated in ac_stars:
                        if camera.ident == tobeactivated:
                            camera.set_star()
                            break

    def activate_fiducials(self):
        """

        Activates all the fiducial markers so their avoidance zones
        are taken into account.

        :Parameters:

        None

        :Returns:

        None

        """
        for fid in self.fiducials:
            if fid is not None:
                fid.activate()

    def deactivate_fiducials(self):
        """

        Deactivates all the fiducial markers so their avoidance zones
        are no longer taken into account.

        :Parameters:

        None

        :Returns:

        None

        """
        for fid in self.fiducials:
            if fid is not None:
                fid.initalise()

    def _get_figure_size(self):
        """

        Helper function which returns a matplotlib figure size appropriate
        to the size of fibre positioner grid.

        """
        # Calculate a figure size based on 4 inches per hexagonal cell.
        width_per_cell = 6.0
        figwidth = width_per_cell * self.columns
        if self.rows > 1:
            figheight = width_per_cell * self.rows * ROOT3BY2
        else:
            figheight = width_per_cell * self.rows

        # Truncate the figure size to a reasonable screen size.
        max_width = 18.0
        max_height = 12.0
        if figwidth > max_width:
            figheight = figheight * max_width/figwidth
            figwidth = max_width
        if figheight > max_height:
            figwidth = figwidth * max_height/figheight
            figheight = max_height

        figsize = (figwidth, figheight)
        #print("Figure size %.1f x %.1f inches" % figsize)
        return figsize

    def plot(self, description='', targetlist=[], withcircles=True, simple=True,
             trivial=False, showplot=True):
        """

        Plot all the positioners in the grid.

        :Parameters:

        description: str (optional)
            Optional description to be added to the plot title.
        targetlist: list of Target objects (optional)
            Optional list of targets whose positions are to be
            shown on the plot.
        withcircles: boolean (optional)
            Set to True to include circles on the plot showing
            the overlapping ranges of the positioners.
        simple: boolean (optional)
            If True, show the arms as thick lines.
            If False, show the arms as realistic rectangles.
        trivial: boolean (optional)
            If True, show all positioners as a simple +.
            If False, plot all positioners.
        showplot: boolean, optional
            If True (the default) show the plot when finished.
            Set of False if you are overlaying several plots, and
            this plot is not the last.

        :Returns:

        plotfig: matplotlib Figure object
            A matplotlib figure containing the plot.

        """
        if plotting is None:
            logger.warning("Plotting is disabled")
            return

        # Determine a figure size appropriate to the grid size,
        # and some sensible line widths.
        figsize = self._get_figure_size()

        title = "Grid of %d positioners in %d cols x %d rows with pitch=%.3f:" % \
            (self.positioner_count, self.max_columns-1, self.max_rows-1, self.pitch)
        if description:
            title += "\n%s" % description
        plotfig = plotting.new_figure(1, figsize=figsize, stitle=title)

        # Plot the positioners.
        for positioner in self.positioners:
            if positioner is not None:
                if trivial:
                    # For a trivial plot, simply mark the location of each
                    # positioner with a +.
                    plotting.plot_xy( positioner.x_centre_focal,
                                      positioner.y_centre_focal,
                                      plotfig=plotfig, linefmt='b+',
                                      linestyle=' ', showplot=False )
                else:
                    # Plot the actual configuration of each positioner.
                    plotfig = positioner.plot(plotfig=plotfig,
                                              showplot=False)

        # If there is a list of fiducials, show these on the plot as well
        if len(self.acquisition_cameras) > 0:
            if trivial:
                xcameras = []
                ycameras = []
                for camera in self.acquisition_cameras:
                    xcameras.append(camera.x_centre_focal)
                    ycameras.append(camera.y_centre_focal)
                plotting.plot_xy( xcameras, ycameras,
                                    plotfig=plotfig, linefmt='ko', linestyle=' ',
                                    showplot=False )
            else:
                for camera in self.acquisition_cameras:
                    # Plot the actual configuration of each cameraucial.
                    plotfig = camera.plot(plotfig=plotfig, showplot=False)

        # If there is a list of fiducials, show these on the plot as well
        if len(self.fiducials) > 0:
            if trivial:
                xfids = []
                yfids = []
                for fid in self.fiducials:
                    xfids.append(fid.x_centre_focal)
                    yfids.append(fid.y_centre_focal)
                plotting.plot_xy( xfids, yfids,
                                    plotfig=plotfig, linefmt='gd', linestyle=' ',
                                    showplot=False )
            else:
                for fid in self.fiducials:
                    # Plot the actual configuration of each fiducial.
                    plotfig = fid.plot(plotfig=plotfig, showplot=False)

        # If there is a list of targets, show these on the plot as well,
        # with x symbols.
        if targetlist:
            xtargets = []
            ytargets = []
            for target in targetlist:
                xtargets.append(target[1])
                ytargets.append(target[2])
                plotting.plot_xy( xtargets, ytargets,
                                  plotfig=plotfig, linefmt='cx', linestyle=' ',
                                  showplot=False )

        if showplot:
            plotting.show_plot()
            plotting.close()
        return plotfig

    def plot_subset(self, ident_list, description=''):
        """

        Plot a subset the positioners in the grid.

        :Parameters:

        ident_list: list of int
            A list of all the positioners to be plotted.
        description: str (optional)
            Optional description to be added to the plot title.

        """
        if plotting is None:
            logger.warning("Plotting is disabled")
            return

        last_positioner = None
        for ident in ident_list:
            positioner = self.get_positioner(ident)
            if positioner is not None:
                last_positioner = positioner
                positioner.plot(showplot=False)
        if last_positioner is not None:
            last_positioner.plot(showplot=True, description=description)

    def plot_status(self, description='', showplot=True):
        """

        Plot the overall status of the positioners in the grid,
        in a similar manner to that expected in the ESO status GUI.

        :Parameters:

        description: str (optional)
            Optional description to be added to the plot title.
        showplot: boolean, optional
            If True (the default) show the plot when finished.
            Set of False if you are overlaying several plots, and
            this plot is not the last.

        :Returns:

        plotfig: matplotlib Figure object
            A matplotlib figure containing the plot.

        """
        if plotting is None:
            logger.warning("Plotting is disabled")
            return

        # Determine a figure size appropriate to the grid size,
        # and some sensible line widths.
        figsize = self._get_figure_size()

        title = "Grid of %d positioners in %d cols x %d rows with pitch=%.3f:" % \
            (self.positioner_count, self.max_columns-1, self.max_rows-1, self.pitch)
        if description:
            title += "\n%s" % description
        plotfig = plotting.new_figure(1, figsize=figsize, stitle=title)

        # Plot the positioners.
        for positioner in self.positioners:
            if positioner is not None:
                # Plot the status of each positioner object.
                plotfig = positioner.plot_status(plotfig=plotfig,
                                                 showplot=False)

        for camera in self.acquisition_cameras:
            if camera is not None:
                # Plot the status of each acquisition camera object.
                plotfig = camera.plot_status(plotfig=plotfig,
                                             showplot=False)
        for fid in self.fiducials:
            if fid is not None:
                # Plot the status of each acquisition camera object.
                plotfig = fid.plot_status(plotfig=plotfig,
                                          showplot=False)

        if showplot:
            plotting.show_plot()
            plotting.close()
        return plotfig

    def count(self):
        """

        Count the number of positioners in the grid.

        :Returns:

        count: int
            The total number of positioners.

        """
        count = 0
        for positioner in self.positioners:
            if positioner is not None:
                count += 1
        return count

    def reset(self):
        """

        Reset all objects in the grid to their initial state.

        """
        for positioner in self.positioners:
            if positioner is not None:
                positioner.initialise()
        for camera in self.acquisition_cameras:
            if camera is not None:
                camera.initialise()
        for fid in self.fiducials:
            if fid is not None:
                fid.initialise()

    def init_counters(self):
        """

        Initialise the conflict counters

        """
        self.counters = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    def counters_str(self):
        """

        Return a summary of the conflict counters in a string.

        """
        strg = "Summary of conflicts:\n"
        strg += "  %d positioners ok\n" % self.counters[CONFLICT_OK]
        strg += "  %d positioners with unreachable targets\n" % self.counters[CONFLICT_UNREACHABLE]
        strg += "  %d locked (broken) positioners\n" % self.counters[CONFLICT_LOCKED]
        strg += "  %d targets too close\n" % self.counters[CONFLICT_TOOCLOSE]
        strg += "  %d conflicts between metrology zones\n" % self.counters[CONFLICT_METROLOGY]
        strg += "  %d conflicts between beta arms\n" % self.counters[CONFLICT_BETA]
        strg += "  %d conflicts between datum actuators\n" % self.counters[CONFLICT_DATUM]
        strg += "  %d conflicts between fibre and beta arm\n" % self.counters[CONFLICT_FIBRE_SNAG]
        strg += "  %d conflicts between fibre and acquisition camera\n" % self.counters[CONFLICT_AC]
        strg += "  %d conflicts between fibre and fiducial marker" % self.counters[CONFLICT_FID]
        return strg

    #+++
    # Class variables defining the location of information within the
    # conflict list.
    CL_ID = 0
    CL_R_OK = 1
    CL_L_OK = 2
    CL_R_DIFF = 3
    CL_L_DIFF = 4
    # Class variables defining the location of information within the
    # reply list.
    R_IDENT = 0
    R_OK = 1
    R_PARITY = 2
    R_DIFFICULTY = 3
    R_ALT_DIFFICULTY = 4
    def check_targets(self, fibretargets, all_neighbours=False,
                      setup_positioners=False):
        """

        Check a collection of targets, with nominated positioners, to
        determine if any will result in conflict between positioners
        in the grid when in their final configuration. The best
        parity is suggested for each positioner.

        Target fibre coordinates are expressed in polar coordinates
        on the MOONS focal plane.

        :Parameters:

        fibretargets: list of (ident, rfibre, thfibre)
            List of target fibre positions required for each
            positioner.
        all_neighbours: boolean (optional)
            If True each positioner is checked for a conflict with all
            its known neighbours, regardless of whether they are included
            in the fibretargets list.
            If False (the default) each positioner is checked for a
            conflict with other positioners in the fibretarget list.
            True is more efficient when checking large collections of
            targets, for example an entire grid configuration.
        setup_positioners: boolean (optional)
            If True, two extra passes through the data are used to setup
            the fibre positioners into the recommended configuration.
            Set to True if the final configuration of the fibre
            positioners is important, or False (the default) if only the
            reply is important and performance is the priority.

        :Returns:

        reply: list of (ident, ok, parity, difficulty, alt_difficulty)
            ident is the fibre positioner identifier
            ok is True if the target can be reached without conflict
            parity is the recommended parity
            difficulty is the degree of difficulty at the recommended parity
            alt_difficulty is the degree of dificulty at the other parity.

        """
        assert isinstance(fibretargets, (tuple,list))
        assert isinstance(fibretargets[0], (tuple,list))
        assert len(fibretargets) > 1
        assert len(fibretargets[0]) == 3

        # First assign all the targets and record whether they can be reached
        # at each parity setting. Make a list of all the positioners included
        # in the list.
        conflict_table = []
        positioner_ids = []
        # For each target
        for (ident, rfibre, thfibre) in fibretargets:
            logger.debug("Checking target %d at r=%f, theta=%f" % \
                         (ident, rfibre, math.degrees(thfibre)))
            # Locate the positioner with the given identifier.
            positioner = self.get_positioner(ident)
            if positioner is not None and isinstance(positioner, Positioner):
                positioner_ids.append(ident)
                # Attempt to assign the target to the positioner at RIGHT parity
                # and record the result
                if positioner.set_target(rfibre, thfibre, parity=PARITY_RIGHT):
                    # Target can be reached
                    ok_right = True
                    diff_right = 0.0
                else:
                    # Target cannot be reached.
                    ok_right = False
                    diff_right = 1.0
                # Attempt to assign the target to the positioner at LEFT parity
                # and record the result
                if positioner.set_target(rfibre, thfibre, parity=PARITY_LEFT):
                    # Target can be reached
                    ok_left = True
                    diff_left = 0.0
                else:
                    # Target cannot be reached.
                    ok_left = False
                    diff_left = 1.0
                conflict_table.append( [ident, ok_right, ok_left,
                                        diff_right, diff_left] )
            elif positioner is None:
                # Positioner is not known
                strg = "Positioner %d is not known." % ident
                logger.error(strg)
                #raise ValueError(strg)
                conflict_table.append( [ident, False, False, 1.0, 1.0] )
            else:
                # Positioner is not known
                strg = "Object %d (%s) is not a fibre positioner." % \
                    (ident, positioner.name)
                logger.error(strg)
                #raise ValueError(strg)
                conflict_table.append( [ident, False, False, 1.0, 1.0] )
        if EXTRA_CONFLICT_INFO:
            logger.info("After first pass: ID R_OK L_OK R_DIFF L_DIFF")
            for entry in conflict_table:
                logger.info(str(entry))

        # Now make a second and third pass through the target list and check
        # each positioner for a potential conflict between neighbours.
        for passtr in ['Second', 'Third']:
            r_record = 0
            for (ident, rfibre, thfibre) in fibretargets:
                logger.debug("%s pass for positioner %d" % (passtr, ident))
                assert( conflict_table[r_record][self.CL_ID] == ident)
                positioner = self.get_positioner(ident)
                # Check for conflict at the RIGHT parity
                if positioner.set_target(rfibre, thfibre, parity=PARITY_RIGHT):
                    if all_neighbours:
                        # Check all neighbours.
                        if positioner.in_conflict_with_neighbours():
                            # Change the status of this parity to not ok.
                            conflict_table[r_record][self.CL_R_OK] = False
                            conflict_table[r_record][self.CL_R_DIFF] = 1.0
                        else:
                            # Confirm this parity is ok.
                            conflict_table[r_record][self.CL_R_OK] = True
                            # Calculate the difficulty for this parity
                            right_diff = positioner.difficulty_value()
                            conflict_table[r_record][self.CL_R_DIFF] = right_diff
                    else:
                        # Only check neighbours which are part of this subset.
                        subset = []
                        neighbours = positioner.get_neighbours()
                        conflict_table[r_record][self.CL_R_OK] = True
                        for neighbour in neighbours:
                            logging.debug("Checking neighbour %s" % neighbour.name)
                            if neighbour.ident in positioner_ids:
                                subset.append(neighbour)
                                if positioner.in_conflict_with(neighbour):
                                    # Change the status of this parity to not ok.
                                    conflict_table[r_record][self.CL_R_OK] = False
                                    conflict_table[r_record][self.CL_R_DIFF] = 1.0
                                    logger.info("RIGHT: Positioner %s in conflict with %s (%s)" % \
                                          (positioner.name, neighbour.name, positioner.conflict_reason))
#                                     if GRAPHICAL_DEBUGGING:
#                                         title = "%s in conflict with %s" % \
#                                             (positioner.name, neighbour.name)
#                                         title += "\n" + positioner.conflict_reason
#                                         plotfig = positioner.plot(showplot=False,
#                                                     description=title)
#                                         neighbour.plot(plotfig=plotfig,
#                                                     showplot=True)
                        if conflict_table[r_record][self.CL_R_OK]:
                            # Calculate the difficulty for this parity
                            right_diff = positioner.difficulty_value(others=subset)
                            logger.debug("RIGHT: Positioner %s ok. Difficulty %f" % \
                                    (positioner.name, right_diff))
                            conflict_table[r_record][self.CL_R_DIFF] = right_diff
                # Check for conflict at the LEFT parity
                if positioner.set_target(rfibre, thfibre, parity=PARITY_LEFT):
                    if all_neighbours:
                        # Check all neighbours.
                        if positioner.in_conflict_with_neighbours():
                            # Change the status of this parity to not ok.
                            conflict_table[r_record][self.CL_L_OK] = False
                            conflict_table[r_record][self.CL_L_DIFF] = 1.0
                        else:
                            # Confirm this parity is ok.
                            conflict_table[r_record][self.CL_L_OK] = True
                            # Calculate the difficulty for this parity
                            left_diff = positioner.difficulty_value()
                            conflict_table[r_record][self.CL_L_DIFF] = left_diff
                    else:
                        # Only check neighbours which are part of this subset.
                        neighbours = positioner.get_neighbours()
                        subset = []
                        conflict_table[r_record][self.CL_L_OK] = True
                        for neighbour in neighbours:
                            logging.debug("Checking neighbour %s" % neighbour.name)
                            if neighbour.ident in positioner_ids:
                                subset.append(neighbour)
                                if positioner.in_conflict_with(neighbour):
                                    # Change the status of this parity to not ok.
                                    conflict_table[r_record][self.CL_L_OK] = False
                                    conflict_table[r_record][self.CL_L_DIFF] = 1.0
                                    logger.info("LEFT: Positioner %s in conflict with %s" % \
                                          (positioner.name, neighbour.name))
#                                     if GRAPHICAL_DEBUGGING:
#                                         title = "%s in conflict with %s" % \
#                                             (positioner.name, neighbour.name)
#                                         title += "\n" + positioner.conflict_reason
#                                         plotfig = positioner.plot(showplot=False,
#                                                     description=title)
#                                         neighbour.plot(plotfig=plotfig,
#                                                     showplot=True)
                        if conflict_table[r_record][self.CL_L_OK]:
                            # Calculate the difficulty for this parity
                            left_diff = positioner.difficulty_value(others=subset)
                            logger.debug("LEFT: Positioner %s ok. Difficulty %f" % \
                                    (positioner.name, left_diff))
                            conflict_table[r_record][self.CL_L_DIFF] = left_diff
                # Move the positioner back its optimum parity
                if conflict_table[r_record][self.CL_R_OK] and \
                   conflict_table[r_record][self.CL_L_OK]:
                    if conflict_table[r_record][self.CL_R_DIFF] <= \
                       conflict_table[r_record][self.CL_L_DIFF]:
                        positioner.set_target(rfibre, thfibre, parity=PARITY_RIGHT)
                    else:
                        positioner.set_target(rfibre, thfibre, parity=PARITY_LEFT)
                elif conflict_table[r_record][self.CL_R_OK]:
                    positioner.set_target(rfibre, thfibre, parity=PARITY_RIGHT)
                elif conflict_table[r_record][self.CL_L_OK]:
                    positioner.set_target(rfibre, thfibre, parity=PARITY_LEFT)

                r_record = r_record + 1

            if EXTRA_CONFLICT_INFO:
                logger.info("After %s pass: ID R_OK L_OK R_DIFF L_DIFF" % passtr)
                for entry in conflict_table:
                    logger.info(str(entry))

        # Now make a fourth pass through the target list, choose the best
        # configuration, and construct the reply.
        reply = []
        r_record = 0
        for (ident, rfibre, thfibre) in fibretargets:
            logger.debug("Fourth pass for positioner %d" % ident)
            assert( conflict_table[r_record][0] == ident)
#             positioner = self.get_positioner(ident)
            if conflict_table[r_record][self.CL_R_OK] and \
               conflict_table[r_record][self.CL_L_OK]:
                # No conflict at either parity. Choose the one with the
                # smallest difficulty.
                if conflict_table[r_record][self.CL_R_DIFF] < \
                   conflict_table[r_record][self.CL_L_DIFF]:
                    # Best choice RIGHT, but LEFT also possible.
                    ok = True
                    parity_choice = PARITY_RIGHT
                    difficulty = conflict_table[r_record][self.CL_R_DIFF]
                    alt_difficulty = conflict_table[r_record][self.CL_L_DIFF]
                else:
                    # Best choice LEFT, but RIGHT also possible.
                    ok = True
                    parity_choice = PARITY_LEFT
                    difficulty = conflict_table[r_record][self.CL_L_DIFF]
                    alt_difficulty = conflict_table[r_record][self.CL_R_DIFF]
            elif conflict_table[r_record][self.CL_R_OK]:
                # No conflict at RIGHT parity only
                ok = True
                parity_choice = PARITY_RIGHT
                difficulty = conflict_table[r_record][self.CL_R_DIFF]
                alt_difficulty = 1.0
            elif conflict_table[r_record][self.CL_L_OK]:
                # No conflict at LEFT parity only
                ok = True
                parity_choice = PARITY_RIGHT
                difficulty = conflict_table[r_record][self.CL_L_DIFF]
                alt_difficulty = 1.0
            else:
                # Conflict
                ok = False
                parity_choice = PARITY_RIGHT
                difficulty = 1.0
                alt_difficulty = 1.0
            reply.append( [ident, ok, parity_choice, difficulty, alt_difficulty] )
            r_record = r_record + 1

        if EXTRA_CONFLICT_INFO:
            logger.info("After fourth pass: ID OK PARITY DIFF ALT_DIFF")
            for entry in reply:
                logger.info(str(entry))

        # Now make a fifth pass through the target list, in the
        # fifth pass moving the positioners to their final configurations
        # and in the sixth pass checking for any unresolved conflict.
        # These last two passes are mainly for simulation purposes, since
        # the reply has already been constructed.
        if setup_positioners:
            r_record = 0
            for (ident, rfibre, thfibre) in fibretargets:
                assert( conflict_table[r_record][self.R_IDENT] == ident)
                positioner = self.get_positioner(ident)
                if positioner is not None and isinstance(positioner, Positioner):
                    positioner.initialise()
                    parity_choice = reply[r_record][self.R_PARITY]
                    logger.debug("Fourth pass for positioner %d at parity %s" % \
                          (ident, util.elbow_parity_str(parity_choice)))
                    positioner.in_conflict = False
                    positioner.set_target(rfibre, thfibre, parity=parity_choice)
                r_record = r_record + 1
            r_record = 0
            for (ident, rfibre, thfibre) in fibretargets:
                assert( conflict_table[r_record][self.R_IDENT] == ident)
                positioner = self.get_positioner(ident)
                if positioner is not None and isinstance(positioner, Positioner):
                    parity_choice = reply[r_record][self.R_PARITY]
                    logger.debug("Fifth pass for positioner %d at parity %s" % \
                          (ident, util.elbow_parity_str(parity_choice)))
                    # Go to the chosen configuration and log any unresolved conflict.
                    if positioner.set_target(rfibre, thfibre, parity=parity_choice):
                        if all_neighbours:
                            if positioner.in_conflict_with_neighbours():
                                logger.warning(positioner.conflict_reason)
                        else:
                            # Only check neighbours which are part of this subset.
                            neighbours = positioner.get_neighbours()
                            for neighbour in neighbours:
                                if neighbour.ident in positioner_ids:
                                    if positioner.in_conflict_with(neighbour):
                                        logger.debug("Positioner %s in conflict with %s" % \
                                                (positioner.name, neighbour.name))
                                        logger.warning(positioner.conflict_reason)
                                        if GRAPHICAL_DEBUGGING:
                                            title = "%s in conflict with %s" % \
                                                (positioner.name, neighbour.name)
                                            title += "\n" + positioner.conflict_reason
                                            plotfig = positioner.plot(showplot=False,
                                                        description=title)
                                            neighbour.plot(plotfig=plotfig, showplot=True)
                    else:
                        logger.warning(positioner.conflict_reason)
                    self.counters[positioner.conflict_type] += 1
                r_record = r_record + 1
        return reply

    def test_pair(self, target1, parity1, target2, parity2):
        """

        Test for a conflict between two targets assigned to specific
        positioners for a particular combination of parities.

        Target fibre coordinates are expressed in polar coordinates
        on the MOONS focal plane.

        :Parameters:

        target1: tuple of (ident, rfibre, thfibre)
            The positioner identifier and proposed target R and theta
            coordinates for the first target.
        parity1: int
            The proposed elbow parity for the first target.
        target2: tuple of (ident, rfibre, thfibre)
            The positioner identifier and proposed target R and theta
            coordinates for the second target.
        parity2: int
            The proposed elbow parity for the second target.

        :Returns:

        (in_conflict, difficulty1, difficulty2)
            in_conflict is True if the two positioners are in conflict,
            or if any of the targets cannot be reached, and False if both
            targets are ok and the positioners are not in conflict.
            difficulty1 and difficulty2 give the degree of difficulty
            estimates for each positioner in the proposed configuration
            (in the range 0.0-1.0).

        """
        assert isinstance(target1, (tuple,list))
        assert isinstance(target2, (tuple,list))
        assert len(target1) >= 3
        assert len(target2) >= 3
        # Unpack the target parameters
        (ident1, rfibre1, thfibre1) = target1
        (ident2, rfibre2, thfibre2) = target2

        # Verify that both positioners exist.
        positioner1 = self.positioners[ident1]
        if positioner1 is None:
            # Positioner is not known
            strg = "Positioner %d is not known." % ident1
            logger.error(strg)
            return (True, 1.0, 1.0)
        positioner2 = self.positioners[ident2]
        if positioner2 is None:
            # Positioner is not known
            strg = "Positioner %d is not known." % ident2
            logger.error(strg)
            return (True, 1.0, 1.0)

        # Initialise both positioners
        positioner1.initialise()
        positioner2.initialise()

        # Attempt to assign the given targets to the positioners.
        if not positioner1.set_target(rfibre1, thfibre1, parity=parity1):
            # Target cannot be reached.
            logger.warning(positioner1.conflict_reason)
            return (True, 1.0, 1.0)
        if not positioner2.set_target(rfibre2, thfibre2, parity=parity2):
            # Target cannot be reached.
            logger.warning(positioner2.conflict_reason)
            return (True, 1.0, 1.0)

        # Test for conflict.
        in_conflict = positioner1.in_conflict_with( positioner2 )
        if in_conflict:
            difficulty1 = 1.0
            difficulty2 = 1.0
        else:
            difficulty1 = positioner1.difficulty_value(others=[positioner2])
            difficulty2 = positioner2.difficulty_value(others=[positioner1])
        return (in_conflict, difficulty1, difficulty2)

    CONFLICT_RR = 0
    CONFLICT_RL = 1
    CONFLICT_LR = 2
    CONFLICT_LL = 3
    def check_pair(self, target1, target2):
        """

        Check for a conflict between two targets assigned to specific
        positioners, trying all four combinations of parities (RR, RL,
        LR, LL) and ranking the ones that don't conflict.

        Target fibre coordinates are expressed in polar coordinates
        on the MOONS focal plane.

        :Paramete,
        rs:

        target1: tuple of (ident, rfibre, thfibre)
            The positioner identifier and proposed target R and theta
            coordinates for the first target.
        target2: tuple of (ident, rfibre, thfibre)
            The positioner identifier and proposed target R and theta
            coordinates for the second target.

        :Returns:

        (rr, rl, lr, ll)
            A conflict list giving the result of the test at each
            parity combination (e.g. rl means target1 at RIGHT parity
            and target2 at LEFT parity). Each entry in the list
            contains a difficulty value, in the range 0.0 to 1.0,
            giving the likelihood of conflict. Values greater than or
            equal to 1.0 indicate conflict.

        """
        # Try each combination of target parities.
        conflict_list = []
        parity_combinations = [(PARITY_RIGHT, PARITY_RIGHT),
                               (PARITY_RIGHT, PARITY_LEFT),
                               (PARITY_LEFT, PARITY_RIGHT),
                               (PARITY_LEFT, PARITY_LEFT),
                               ]
        for (parity1, parity2) in parity_combinations:
            (in_conflict, diff1, diff2) = self.test_pair(target1, parity1,
                                                         target2, parity2)
            if in_conflict:
                conflict_list.append( 1.0 )
            else:
                combined_diff = (diff1 + diff2)/2.0
                conflict_list.append( combined_diff )

        return conflict_list


def count_conflicts( conflict_reply ):
    """

    Count the number of conflicts in a conflict check reply list.

    :Parameters:

    conflict_reply: list of (ident, ok, parity, difficulty, alt_difficulty)
        The reply list returned by the PositionerGrid  "check_targets"
        method.

    :Returns:

    (ngood, nconflicts): list of int
        The number of good replies and the number of conflicts.

    """
    ngood = 0
    nconflicts = 0
    for entry in conflict_reply:
        if entry[PositionerGrid.R_OK]:
            ngood += 1
        else:
            nconflicts += 1
    return (ngood, nconflicts)


if __name__ == "__main__":
    print("\nMOONS OPS/FPS shared library test...")

    PLOTTING = True
    if plotting is None:
        PLOTTING = False

    # Decide which tests to run.
    POSITIONER_TESTS = True
    MOTOR_TESTS = True
    POSITIONER_GRID_TESTS = True
    POSITIONER_PAIR_TESTS = True
    TEST_CASES = True

    if POSITIONER_TESTS:
        def set_target_angles( positioner, r, theta, parity, description):
            """

            Test function to set a positioner to a given target, display
            and then plot the configuration, including the given description
            in the title.

            """
            if r is not None and theta is not None and parity is not None:
                if not positioner.set_target(r, math.radians(theta), parity):
                    print("***Failed to set target (%f,%f)" % (r,theta))
            print("")
            print(positioner)
            strg = "Arm angles are alpha=%.3f, beta=%.3f" % positioner.get_arm_angles()
            strg += " (orient=%.3f)" % positioner.orient
            print(strg + " (" + description + ").")
            if PLOTTING:
                positioner.plot(showplot=True,
                                description=description + "\n" + strg)

        def set_motor_angles( positioner, alpha, beta, description):
            """

            Test function to set a positioner's arms to given orientations,
            display and then plot the configuration, including the given
            description in the title.

            """
            if alpha is not None and beta is not None:
                positioner.set_arm_angles(alpha, beta)
            print("")
            print(positioner)
            strg = "Arm angles are alpha=%.3f, beta=%.3f" % positioner.get_arm_angles()
            strg += " (orient=%.3f)" % math.degrees(positioner.orient)
            print(strg + " (" + description + ").")
            if PLOTTING:
                positioner.plot(showplot=True,
                                description=description + "\n" + strg)

        # Create a positioner and move its arms to different orientations.
        if MOTOR_TESTS:
            print("\n----------------------------------------------------------------------")
            print("Positioner motor tests...")
            import numpy as np
            positioner = Positioner(1, 0.0, 0.0, math.radians(0.0), 17, 19)
            beta = BETA_DEFAULT
            for alpha in np.arange(ALPHA_LIMITS[0], ALPHA_LIMITS[1], 60.0):
                set_motor_angles( positioner, alpha, beta, "Moving alpha")
            set_motor_angles( positioner, ALPHA_LIMITS[1], beta, "Moving alpha")
            positioner.initialise()
            alpha = ALPHA_DEFAULT
            for beta in np.arange(BETA_LIMITS[0], BETA_LIMITS[1], 60.0):
                set_motor_angles( positioner, alpha, beta, "Moving beta")
            set_motor_angles( positioner, alpha, BETA_LIMITS[1], "Moving beta")
            del positioner
            for orient in [0.0, 45.0, 225.0]:
                positioner = Positioner(1, 0.0, 0.0, math.radians(orient), 17, 19)
                set_motor_angles( positioner, 42.0, -120.0, "Changing orientation")
                del positioner

        # Create a positioner, assign a target to it and then plot it,
        # showing the avoidance zones used by the conflict detection
        # function.
        # ident, r_centre_focal, theta_centre_focal, orient, column, row
        print("\n----------------------------------------------------------------------")
        print("Positioner tests...")
        positioner = Positioner(1, 0.0, 0.0, math.radians(0.0), 17, 19)
        strg = "Inner value=%f, " % positioner.inner
        strg += "outer value=%f " % positioner.outer
        strg += "and critical length=%f" % positioner.lengthcrit
        print(strg)
        set_target_angles( positioner, None, None, None, "Starting position")
        set_target_angles( positioner, positioner.outer, 90.0, PARITY_RIGHT,
                           "Fully stretched position")
        set_target_angles( positioner, positioner.lengthcrit, 90.0, PARITY_LEFT,
                           "Opposite fully stretched position")
        set_target_angles( positioner, 20.0, 15.0, PARITY_RIGHT,
                           "r=20.0 theta=15.0 with right-armed parity")
        set_target_angles( positioner, 20.0, 15.0, PARITY_LEFT,
                           "r=20.0 theta=15.0 with left-armed parity")
        set_target_angles( positioner, 12.0, 180.0, PARITY_RIGHT,
                           "r=12.0 theta=180.0 with right-armed parity")
        set_target_angles( positioner, 12.0, 180.0, PARITY_LEFT,
                           "r=12.0 theta=180.0 with left-armed parity")
        set_target_angles( positioner, 16.0, -65.0, PARITY_RIGHT,
                           "r=16.0 theta=-65.0 with right-armed parity")
        set_target_angles( positioner, 16.0, -65.0, PARITY_LEFT,
                           "r=16.0 theta=-65.0 with left-armed parity")

        print("\nAdd a second positioner and check for conflict")
        positioner2 = Positioner(2, 20, math.radians(45.0), 0.0, 16, 21)
        set_target_angles(positioner2, 30.0, 30.0, PARITY_RIGHT,
                        "r=30.0 theta=30.0 with right-armed parity");

        # Test for conflict
        print("\nin_conflict_with function test.")
        if (positioner2.in_conflict_with( positioner ) ):
            print("Positioners are in conflict (" + \
                  positioner.conflict_reason + ")")
        else:
            print("Positioners are NOT in conflict")

        print("\ncheck_for_conflict function test.")
        in_conflict = check_for_conflict(0.0, 0.0, 0.0,
                            16.0, math.radians(-65.0), PARITY_LEFT,
                            20.0, math.radians(45.0), 0.0,
                            30.0, math.radians(30.0), PARITY_RIGHT)
        if ( in_conflict ):
            print("Positioners are in conflict (" + positioner.conflict_reason + ")")
        else:
            print("Positioners are NOT in conflict")

    if POSITIONER_GRID_TESTS:
        print("\n----------------------------------------------------------------------")
        print("PositionerGrid test...\n")
        # Define a fibre positioner grid in Polar coordinates
        positioner_configs = [
                # id,   rcen,     thcen,        orient,      col, row
#                 ( 1,  0.00000, math.radians( 90.0), math.radians(0.0), 17, 19),
                ( 2, 25.25870, math.radians( 30.0), math.radians(0.0), 17, 20),
                ( 3, 25.25870, math.radians( 90.0), math.radians(0.0), 18, 19),
                ( 4, 25.25872, math.radians(150.0), math.radians(0.0), 17, 18),
                ( 5, 25.25875, math.radians(210.0), math.radians(0.0), 16, 18),
                ( 6, 25.25870, math.radians(270.0), math.radians(0.0), 16, 19),
                ( 7, 25.25873, math.radians(330.0), math.radians(0.0), 16, 20),
                ( 8, 43.74782, math.radians(  0.0), math.radians(0.0), 17, 21),
                ( 9, 50.51474, math.radians( 30.0), math.radians(0.0), 18, 21),
                (10, 43.74782, math.radians( 60.0), math.radians(0.0), 18, 20),
                (11, 50.51474, math.radians( 90.0), math.radians(0.0), 19, 19),
                (12, 43.74781, math.radians(120.0), math.radians(0.0), 18, 18),
                (13, 50.51470, math.radians(150.0), math.radians(0.0), 18, 17),
                (14, 43.74780, math.radians(180.0), math.radians(0.0), 17, 17),
                (15, 50.51471, math.radians(210.0), math.radians(0.0), 16, 17),
                (16, 43.74779, math.radians(240.0), math.radians(0.0), 15, 18),
                (17, 50.51470, math.radians(270.0), math.radians(0.0), 15, 19),
                (18, 43.74780, math.radians(300.0), math.radians(0.0), 15, 20),
                (19, 50.51476, math.radians(330.0), math.radians(0.0), 16, 21)]

        # Define some acquisition cameras
        ac_configs = [
                # id,   rcen,     thcen,              orient,          fov, col, row
                ( 1,  0.00000, math.radians( 90.0), math.radians(0.0), AC_FOV_RADIUS, 17, 19),
            ]

        # Define a guide star
        ac_guide_stars = [
                (1, 3.0, 136.0, 12.0)
            ]

        fid_configs = [
                # id,  xcen, ycen,    orient,        sminor, smajor
                ( 1,   0.00, -75.00, math.radians(0.0), FID_SEMI_MINOR, FID_SEMI_MAJOR),
                ( 2,  64.95, -37.50, math.radians(0.0), FID_SEMI_MINOR, FID_SEMI_MAJOR),
                ( 3,  64.95,  37.50, math.radians(0.0), FID_SEMI_MINOR, FID_SEMI_MAJOR),
                ( 4,   0.00,  75.00, math.radians(0.0), FID_SEMI_MINOR, FID_SEMI_MAJOR),
                ( 5, -64.95,  37.50, math.radians(0.0), FID_SEMI_MINOR, FID_SEMI_MAJOR),
                ( 6, -64.95, -37.50, math.radians(0.0), FID_SEMI_MINOR, FID_SEMI_MAJOR)
            ]

        # Create a fibre positioner grid object.
        positioner_grid = PositionerGrid( positioner_configs,
                                          camera_list=ac_configs,
                                          fiducial_list=fid_configs )
        positioner_grid.define_guide_stars(ac_guide_stars)
        positioner_grid.activate_fiducials()
        print( positioner_grid )

        print("\nConflict check test 1 (check_targets with 19 positioners)...\n")
        # Define a set of targets in Cartesian coordinates (for convenience),
        # then convert to polar coordinates to match the agreed interface.
        fibrecoords = [(2,   26.50, 42.0),
                       (3,   22.94, -16.0),
                       (4,    8.80, -33.0),
                       (5,    6.86, -32.0),
                       (6,  -11.57, 14.0),
                       (7,  -16.00,  0.0),
                       (8,  -20.0, 60.0),
                       (9,   16.2, 43.75),
                       (10,  48.0, 19.0),
                       (11,  45.0, 16.0),
                       (12,  40.87, -11.87),
                       (13,  10.0, -38.5),
                       (14, -24.9, -43.75),
                       (15, -25.0, -60.0 ),
                       (16, -50.7, -10.0),
                       (17, -54.0, -18.0),
#                        (18, -37.0,  40.0),   # Original location
                       (18, -55.0,  35.0),   # Conflicts with fiducial
                       (19, -20.0,  35.0)]
        fibretargets = []
        for (ident, x, y) in fibrecoords:
            (r, theta) = util.cartesian_to_polar(x, y)
            fibretargets.append( (ident, r, theta) )

        time.sleep(1) # Ensure any logging text appears before the reply string.
        strg = "Target suggestions [ident, r, theta] (focal plane coordinates):"
        for entry in fibretargets:
            theta_degrees = math.degrees(entry[2])
            strg += "\n  [%d, %.5f, %.5f (%.5f deg)]" % \
                (entry[0], entry[1], entry[2], theta_degrees)
        print(strg)

        # Now assign these targets to the fibre positioner grid and
        # check for conflicts.
        replies = positioner_grid.check_targets(fibretargets,
                                                all_neighbours=True,
                                                setup_positioners=True)
        print(positioner_grid.counters_str())

        time.sleep(1) # Ensure any logging text appears before the reply string.
        strg = "Reply string 1 [ident, ok, parity, difficulty, alt_difficulty]: "
        for entry in replies:
            strg += "\n  [%d, %s, %s, %.4f, %.4f]" % \
                (entry[0], str(entry[1]), util.elbow_parity_str(entry[2]),
                 entry[3], entry[4])
        print(strg)

        targetlist = []
        for ident, rtarget, thtarget in fibretargets:
            xtarget, ytarget = util.polar_to_cartesian(rtarget, thtarget)
            targetlist.append( [ident, xtarget, ytarget] )
        if PLOTTING:
            plotfig = positioner_grid.plot_status()
            plotfig = positioner_grid.plot(targetlist=targetlist,
                            description="(Green=no target; " + \
                            "Magenta=target unreachable; Black=target assigned; " + \
                            "Red=In conflict; x=targets)")

        print("\nConflict check test 2 (check_targets with just 2 positioners)...\n")
        # Define two simple targets for 16 and 17
        (r1, theta1) = util.cartesian_to_polar(-50.7, -10.0)
        (r2, theta2) = util.cartesian_to_polar(-54.0, -18.0)
        fibretargets = [(16, r1, theta1),
                        (17, r2, theta2)]

        # Now assign these targets to the fibre positioner grid and
        # check for conflicts.
        replies = positioner_grid.check_targets(fibretargets)

        time.sleep(1) # Ensure any logging text appears before the reply string.
        strg = "Reply string 2 [ident, ok, parity, difficulty, alt_difficulty]: "
        for entry in replies:
            strg += "\n  [%d, %s, %s, %.4f, %.4f]" % \
                (entry[0], str(entry[1]), util.elbow_parity_str(entry[2]),
                 entry[3], entry[4])
        print(strg)
    elif POSITIONER_PAIR_TESTS:
        # Define the minimum needed by POSITIONER_PAIR_TESTS
        positioner_configs_xy = [
                # id,   x,         y,           orient,      col, row
                (16, -37.88670, -21.87390, math.radians(0.0), 15, 18),
                (17, -50.51470,   0.00000, math.radians(0.0), 15, 19)
            ]
        # Convert to a grid configuration in polar coordinates.
        positioner_configs = []
        for config in positioner_configs_xy:
            newconfig = []
            newconfig.append(config[0])
            (r,theta) = util.cartesian_to_polar(config[1], config[2])
            newconfig.append(r)
            newconfig.append(theta)
            newconfig.append(config[3])
            newconfig.append(config[4])
            newconfig.append(config[5])
            positioner_configs.append(newconfig)
        # Create a fibre positioner grid object.
        positioner_grid = PositionerGrid( positioner_configs )
        print( positioner_grid )

    if POSITIONER_PAIR_TESTS:
        print("\nConflict check test 3 (test_pair)...\n")
        # Try the same pair of targets in 4 different parity configurations.
        (r1, theta1) = util.cartesian_to_polar(-50.7, -10.0)
        (r2, theta2) = util.cartesian_to_polar(-54.0, -18.0)
        target1 = (16, r1, theta1)
        target2 = (17, r2, theta2)
        parity_checks = [[PARITY_RIGHT, PARITY_RIGHT],
                         [PARITY_RIGHT, PARITY_LEFT],
                         [PARITY_LEFT,  PARITY_RIGHT],
                         [PARITY_LEFT,  PARITY_LEFT]
                         ]
        for (parity1, parity2) in parity_checks:
            # Construct a string from the first letter of each parity name.
            pstrg1 = util.elbow_parity_str(parity1)
            pstrg2 = util.elbow_parity_str(parity2)
            paritystrg = "%.1s%.1s" % (pstrg1, pstrg2)
            (conflict, diff1, diff2) = positioner_grid.test_pair( target1, parity1,
                                                                  target2, parity2 )
            if conflict:
                strg = "Positioners 16 and 17 are in conflict (%s)." % paritystrg
            else:
                strg = "Positioners 16 and 17 are not in conflict (%s)." % paritystrg
                strg += " Difficulty estimates: %.3f, %.3f" % (diff1, diff2)
                strg += " (mean %.3f)." % ((diff1 + diff2)/2.0)
            print( strg )
            if PLOTTING:
                positioner_grid.plot_subset( [16, 17], description=strg )

        print("\nConflict check test 4 (check_pair)...\n")
        conflict_list = positioner_grid.check_pair(target1, target2)
        print("Conflict list: rr=%.3f, rl=%.3f, lr=%.3f, ll=%.3f" % tuple(conflict_list))

    # Fibre positioner and path analysis test cases.
    # Awkward situations that the fibre positioner and path analysis
    # software must be able to cope with.
    if TEST_CASES:
        positioner_configs = [
            # id,     r,       theta,          orient,           col, row
            (1,  43.7478, math.radians(180.0), math.radians(0.0), 1, 1),
            (2,  50.5147, math.radians(210.0), math.radians(0.0), 2, 1)
            ]
        positioner_grid2 = PositionerGrid( positioner_configs )

        # Fibre positioner test case 1 - demonstrated at RR
        (r1, theta1) = util.cartesian_to_polar(-9.0, -44.0) # 44.911, 191.560
        (r2, theta2) = util.cartesian_to_polar(-1.1, -49.0) # 49.012, 181.286
        fp_test_case_1 = [(1, r1, theta1),
                          (2, r2, theta2)]
        # Fibre positioner test case 2 - demonstrated at RR
        (r1, theta1) = util.cartesian_to_polar(-9.0, -44.0)  # 44.911, 191.560
        (r2, theta2) = util.cartesian_to_polar(-12.0, -47.0) # 48.508, 194.323
        fp_test_case_2 = [(1, r1, theta1),
                          (2, r2, theta2)]
        # Fibre positioner test case 3 - demonstrated at RL
        (r1, theta1) = util.cartesian_to_polar(-9.0, -43.5)  # 44.421, 191.689
        (r2, theta2) = util.cartesian_to_polar(-2.0, -38.0)  # 38.053, 183.013
        fp_test_case_3 = [(1, r1, theta1),
                          (2, r2, theta2)]
        # Fibre positioner test case 4 - demonstrated at RR
        (r1, theta1) = util.cartesian_to_polar(-23.5, -50.994) # 56.148, 204.742
        (r2, theta2) = util.cartesian_to_polar(-34.3, -43.994) # 55.785, 217.942
        fp_test_case_4 = [(1, r1, theta1),
                          (2, r2, theta2)]
        # Fibre positioner test case 5 - demonstrated at RR
        (r1, theta1) = util.cartesian_to_polar(-10.0, -38.0) # 39.294, 194.744
        (r2, theta2) = util.cartesian_to_polar(-0.3, -43.7)  # 43.701, 180.393
        fp_test_case_5 = [(1, r1, theta1),
                          (2, r2, theta2)]
        # Fibre positioner test case 6 - demonstrated at RR
        (r1, theta1) = util.cartesian_to_polar(-22.0, -55.0) # 59.237, 201.801
        (r2, theta2) = util.cartesian_to_polar(-15.8, -47.0) # 49.585, 198.581
        fp_test_case_6 = [(1, r1, theta1),
                          (2, r2, theta2)]
        # Fibre positioner test case 7 - demonstrated at RR (LR avoids conflict)
        (r1, theta1) = util.cartesian_to_polar(5.0, -54.5)  # 54.729, 174.758
        (r2, theta2) = util.cartesian_to_polar(-2.5, -50.0) # 50.062, 182.862
        fp_test_case_7 = [(1, r1, theta1),
                          (2, r2, theta2)]
        # Fibre positioner test case 8 - demonstrated at RR (LL avoids conflict)
        (r1, theta1) = util.cartesian_to_polar(-0.5, -59.5) # 59.502, 180.481
        (r2, theta2) = util.cartesian_to_polar(-5.0, -37.5) # 37.832, 187.595
        fp_test_case_8 = [(1, r1, theta1),
                          (2, r2, theta2)]
        # Fibre positioner test case 9 - demonstrated at RR (LL has smaller difficulty)
        (r1, theta1) = util.cartesian_to_polar(-0.5, -59.5) # 59.502, 180.481
        (r2, theta2) = util.cartesian_to_polar(-2.0, -49.0) #
        fp_test_case_9 = [(1, r1, theta1),
                          (2, r2, theta2)]
        # Fibre positioner test case 10 - demonstrated at RR (LL avoids conflict)
        (r1, theta1) = util.cartesian_to_polar(-0.5, -59.5) # 59.502, 180.481
        (r2, theta2) = util.cartesian_to_polar(-3.0, -43.0) # 43.105, 183.991
        fp_test_case_10 = [(1, r1, theta1),
                           (2, r2, theta2)]
        # Fibre positioner test case 11 - demonstrated at RR
        (r1, theta1) = util.cartesian_to_polar(-20.0, -57.0) #
        (r2, theta2) = util.cartesian_to_polar(-21.0, -52.0) #
        fp_test_case_11 = [(1, r1, theta1),
                           (2, r2, theta2)]

        # Path analysis test case 1
        pa_test_case_1 = [(1, 40.3434, math.radians(212.668)),
                          (2, 28.0970, math.radians(219.449))]
        # Path analysis test case 2
        pa_test_case_2 = [(1, 50.3396, math.radians(209.646)),
                          (2, 65.0000, math.radians(202.620))]

        test_cases = [
                    fp_test_case_1,
                    fp_test_case_2,
                    fp_test_case_3,
                    fp_test_case_4,
                    fp_test_case_5,
                    fp_test_case_6,
                    fp_test_case_7,
                    fp_test_case_8,
                    fp_test_case_9,
                    fp_test_case_10,
                    fp_test_case_11,
                    pa_test_case_1,
                    pa_test_case_2
                    ]

        parity_checks = ([PARITY_RIGHT, PARITY_RIGHT],
                         [PARITY_RIGHT, PARITY_LEFT],
                         [PARITY_LEFT,  PARITY_RIGHT],
                         [PARITY_LEFT,  PARITY_LEFT]
                         )
        #parity_checks = [[PARITY_RIGHT, PARITY_RIGHT]]

        ncase = 1
        for (target1, target2) in test_cases:
            test_title = "Test case %d: " % ncase
            print("\n" + test_title + "...\n")

            for (parity1, parity2) in parity_checks:
                # Construct a string from the first letter of each parity name.
                pstrg1 = util.elbow_parity_str(parity1)
                pstrg2 = util.elbow_parity_str(parity2)
                paritystrg = "%.1s%.1s" % (pstrg1, pstrg2)
                (conflict, diff1, diff2) = positioner_grid2.test_pair(
                                                        target1, parity1,
                                                        target2, parity2 )
                if conflict:
                    strg = "Positioners 1 and 2 in conflict (%s).\n" % paritystrg
                    pos1 = positioner_grid2.get_positioner(1)
                    if pos1.conflict_reason:
                        strg += pos1.conflict_reason
                    pos2 = positioner_grid2.get_positioner(2)
                    if pos2.conflict_reason:
                        strg += " " + pos2.conflict_reason
                else:
                    strg = "Positioners 1 and 2 not in conflict (%s)." % paritystrg
                    strg += " Difficulty estimates: %.3f, %.3f" % (diff1, diff2)
                    strg += " (mean %.3f)." % ((diff1 + diff2)/2.0)
                print( strg )
                print( positioner_grid2 )
                if PLOTTING:
                    positioner_grid2.plot(description=test_title+strg)
            ncase += 1

    print("Test finished.")
