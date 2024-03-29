# -*- coding: utf-8 -*-
"""
/***************************************************************************
 DirSlope
                                 A QGIS plugin
 This plugin calculates directional slope along a line. A DEM is required.
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                             -------------------
        begin                : 2024-03-26
        copyright            : (C) 2024 by Joe Tayabji, Dara Seidl
        email                : jtayabji@gmail.com
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load DirSlope class from file DirSlope.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .dir_slope import DirSlope
    return DirSlope(iface)
