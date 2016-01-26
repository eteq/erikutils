from __future__ import division, print_function
"""
Velocity transformation functions
"""

__all__ = ['vhelio_to_vlsr', 'vlsr_to_vgsr',
            'vhelio_to_lsr_term', 'vlsr_to_gsr_term'] #these deprecated

import warnings

import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord, CIRS, AltAz
from astropy.coordinates.builtin_frames.utils import get_polar_motion
from astropy import _erfa as erfa

@u.quantity_input(vhelio=u.km/u.s, uvwlsr=u.km/u.s)
def vhelio_to_vlsr(coord, vhelio=0*u.km/u.s, lsr_definition='dynamical'):
    """
    Convert heliocentric radial velocity to Local Standard of Rest radial velocity

    Parameters
    ----------
    coord : SkyCoord
        The direction at which to compute the conversion factor.
    vhelio : Quantity with velocity units
        The heliocentric radial velocity.  Should be a scalar quantity or match the shape of `coord`.
    lsr_definition : str or Quantity with velocity units
        The definition of LSR to assume.  This can be one of three options:
        * 'kinematic' : 20 km/s towards 18h,+30 deg (at 1900)
        * 'dynamical' : IAU definition of (9, 12, 7) km/sec in Galactic cartesian coordinates
        * A length-3 vector with the U, V, and W components (velocity along the galactic x, y, and z axes)
          to assume for the LSR.

    Returns
    -------
    vlsr : Quantity
        The velocity in the Local Standard of Rest frame.
    """
    if lsr_definition == 'kinematic':
        direction = SkyCoord('18h', '30d', frame='fk5', equinox='J1900')
        velocity = 20*u.km/u.s
        uvw_lsr = direction.galactic.cartesian.xyz * velocity
    elif lsr_definition == 'dynamical':
        uvw_lsr = (9, 12, 7)*u.km/u.s
    else:
        uvw_lsr = lsr_definition

    # the unitspherical conversion ensures that the resulting cartesian is a *unit* vector
    usr = coord.galactic.represent_as(UnitSphericalRepresentation)
    cart = usr.to_cartesian()

    vlsr_vector = (cart.xyz.T*uvw_lsr).T

    return vhelio + np.sum(vlsr_vector)


@u.quantity_input(vlsr=u.km/u.s, vrot=u.km/u.s)
def vlsr_to_vgsr(coord, vlsr=0*u.km/u.s, vrot=220*u.km/u.s):
    """
    Convert Local Standard of Rest radial velocity to galactocentric velocity

    Parameters
    ----------
    coord : SkyCoord
        The direction at which to compute the conversion factor.
    vlsr : Quantity with velocity units
        The heliocentric radial velocity.  Should be a scalar quantity or match the shape of `coord`.
    vrot : Quantity with velocity units
        The value of the galactic rotation curve at the solar cycle. The default is the IAU convention of 220 km/s.

    Returns
    -------
    vlsr : Quantity
        The velocity in the Local Standard of Rest frame.
    """
    # this adds the projection along the l=90, b=0 axis.  Could do the whole vector approach as above
    # but along the axis makes it trivial to do this simpler version
    l = coord.galactic.l
    b = coord.galactic.b
    return vlsr + vrot * np.sin(l) * np.cos(b)


# old versions

def vhelio_to_lsr_term(coo):
    """
    From http://www.atnf.csiro.au/people/Tobias.Westmeier/tools_hihelpers.php#restframes

    *Add* this to vhelio to get lsr
    """
    warnings.warn(DeprecationWarning('vhelio_to_lsr_term was replaced by vhelio_to_vlsr'))
    gcoo = coo.transform_to('galactic')
    sb = np.sin(gcoo.b)
    cb = np.cos(gcoo.b)
    sl = np.sin(gcoo.l)
    cl = np.cos(gcoo.l)
    return (9*cl*cb + 12*sl*cb + 7*sb)*u.km/u.s


def vlsr_to_gsr_term(coo, vrot=220*u.km/u.s):
    """
    From http://www.atnf.csiro.au/people/Tobias.Westmeier/tools_hihelpers.php#restframes

    *Add* this to vlsr to get gsr
    """
    warnings.warn(DeprecationWarning('vlsr_to_gsr_term was replaced by vlsr_to_vgsr'))
    gcoo = coo.transform_to('galactic')
    return vrot * np.sin(gcoo.l) * np.cos(gcoo.b)
