# coding: utf-8
"""
Source Name:   generatepointsfromlines.py
Version:       ArcGIS Pro 3.1
Author:        Environmental Systems Research Institute Inc.
Description:   Source for Generate Points From Line geoprocessing tool.
"""

import arcpy
import os
from collections import namedtuple

point_placement = dict(DISTANCE=False, PERCENTAGE=True)


def create_points_from_lines(input_fc, output_fc, spatial_ref, percent=False,
                             dist=True, add_end_points=False, add_chainage=False):
    """Convert line features to feature class of points

    :param input_fc: Input line features
    :param output_fc: Output point feature class
    :param spatial_ref: The spatial reference of the input
    :param percent: If creating points by percentage (or distance)
    :param dist: The distance used to create points (if percentage == False).
    The distance should be in units of the input (see convert_units)
    :param add_end_points: If True, an extra point will be added from start
    and end point of each line feature
    :param add_chainage: If True, fields with sequence and distance are
    added to the output
    :return: None
    """

    if percent:
        is_percentage = True
    else:
        is_percentage = False

    create_feature_class(input_fc, output_fc, spatial_ref)

    fid_name = 'ORIG_FID'
    len_name = 'ORIG_LEN'
    seq_name = 'ORIG_SEQ'

    # Add necessary fields
    if add_chainage:
        arcpy.management.AddFields(
            output_fc,
            [[fid_name, 'LONG'], [len_name, 'DOUBLE'], [seq_name, 'LONG']])
    else:
        arcpy.management.AddField(output_fc, fid_name, 'LONG')

    # Create new points based on input lines
    in_fields = ['SHAPE@', 'OID@']
    out_fields = ['SHAPE@', fid_name]
    if add_chainage:
        out_fields += [len_name, seq_name]

    out_fc_is_empty = True
    with arcpy.da.SearchCursor(input_fc, in_fields) as search_cursor:
        with arcpy.da.InsertCursor(output_fc, out_fields) as insert_cursor:
            for row in search_cursor:
                line = row[0]

                if line:  # if null geometry--skip
                    i = 1
                    if line.type == 'polygon':
                        line = line.boundary()

                    if add_end_points:
                        out_fc_is_empty = False
                        insert_values = [line.firstPoint, row[1]]
                        if add_chainage:
                            insert_values += [0, i]
                        insert_cursor.insertRow(insert_values)
                        i += 1

                    increment = (percent or dist)
                    cur_length = increment

                    if is_percentage:
                        max_position = 1.0
                    else:
                        max_position = line.length

                    while cur_length < max_position:
                        out_fc_is_empty = False
                        new_point = line.positionAlongLine(cur_length,
                                                           is_percentage)
                        insert_values = [new_point, row[1]]
                        if add_chainage:
                            if is_percentage:
                                insert_values += [line.queryPointAndDistance(new_point)[1], i]
                            else:
                                insert_values += [cur_length, i]
                        insert_cursor.insertRow(insert_values)
                        i += 1
                        cur_length += increment

                    if add_end_points:
                        end_point = line.positionAlongLine(1, True)
                        insert_values = [end_point, row[1]]
                        if add_chainage:
                            insert_values += [line.length, i]
                        insert_cursor.insertRow(insert_values)

        try:
            oid_name = get_OID_name(input_fc)
            arcpy.management.JoinField(out_fc, fid_name, input_fc, oid_name)
        except arcpy.ExecuteError:
            # In unlikely event that JoinField fails, proceed regardless,
            # as spatial and join field are already complete
            pass

    if out_fc_is_empty:
        arcpy.AddIDMessage('WARNING', 117)

    return


def create_feature_class(input_fc, output_fc, spatial_ref):
    """Create the initial feature class

    :param input_fc: Input line features
    :param output_fc: Output point feature class
    :param spatial_ref: The spatial reference of the input
    """

    desc = arcpy.Describe(input_fc)

    # Take flag environment over Describe property unless set to default
    support_m = arcpy.env.outputMFlag.upper() if arcpy.env.outputMFlag in ['Enabled', 'Disabled'] \
        else "ENABLED" if desc.hasM else "DISABLED"
    support_z = arcpy.env.outputZFlag.upper() if arcpy.env.outputZFlag in ['Enabled', 'Disabled'] \
        else "ENABLED" if desc.hasZ else "DISABLED"

    # Create output feature class
    arcpy.management.CreateFeatureclass(
        os.path.dirname(output_fc),
        os.path.basename(output_fc),
        geometry_type="POINT",
        has_m=support_m,
        has_z=support_z,
        spatial_reference=spatial_ref)

    return


def get_OID_name(in_data):
    """Get the OIDFieldName of the data source

    :param in_data: Input data source
    :return: The OID field name of the data source
    """

    d = arcpy.Describe(in_data)
    oid = getattr(d, 'OIDFieldName') \
          if hasattr(d, 'OIDFieldName') \
          else getattr(arcpy.Describe(d.catalogPath), 'OIDFieldName')
    return oid


def convert_units(dist, param_units, spatial_info):
    """Base unit conversion

    :param dist: Distance
    :param param_units: The units as supplied from tool parameter
    :param spatial_info: arcpy.SpatialReference object
    :return: Distance after converted to new units
    """

    param_units = param_units.upper()

    if param_units in ['', None, 'UNKNOWN']:
        return dist
    else:
        if param_units != 'DECIMALDEGREES':
            p_conversion = arcpy.LinearUnitConversionFactor(param_units, 'meters')
        else:
            p_conversion = 111319.8

        try:
            sr_conversion = spatial_info.spatialReference.metersPerUnit
        except AttributeError:
            try:
                input_extent = spatial_info.extent

                centroid = input_extent.polygon.centroid
                point1 = centroid.Y, centroid.X - 0.5
                point2 = centroid.Y, centroid.X + 0.5
                sr_conversion = haversine(point1, point2) * 1000
            except Exception as err:
                # Fallback
                sr_conversion = 111319.8

        return dist * (p_conversion / sr_conversion)


def get_distance_and_units(dist):
    """ Pull distance and units from a linear unit. If units are not
    specified, return UNKNOWN.

    :param dist: Linear units
    :return: Tuple of distance (float) and units (string)
    """
    try:
        dist, units = dist.split(' ', 1)
    except ValueError:
        # ValueError occurs if units are not specified, use 'UNKNOWN'
        units = 'UNKNOWN'

    dist = dist.replace(',', '.')

    return float(dist), units


def haversine(point1, point2):
    """ Calculate the distance between two points on the Earth surface around its curvature.
    Does not account for changes in elevation (datum)

    :param point1 Tuple - Tuple of (Lat, Long) for the first point
    :param point2 Tuple - Tuple of (Lat, Long) for the second point
    :return Float - The distance between the two points about the surface of the globe in kilometers.
    """
    from math import radians, sin, cos, asin, sqrt
    radius_of_earth_km = 6371
    lat1, lng1, lat2, lng2 = list(map(radians, list(point1 + point2)))
    d = sin((lat2 - lat1) / 2) ** 2 + cos(lat1) * cos(lat2) * sin((lng2 - lng1) / 2) ** 2
    return 2 * radius_of_earth_km * asin(sqrt(d))


if __name__ == '__main__':
    in_features = arcpy.GetParameterAsText(0)  # String
    out_fc = arcpy.GetParameterAsText(1)  # String
    use_percent = point_placement[arcpy.GetParameter(2)]  # Str -> Bool
    end_points = arcpy.GetParameter(5)  # Boolean
    chainage = arcpy.GetParameter(6)  # Boolean

    describe = arcpy.Describe(in_features)
    spatial_info = namedtuple('spatial_info', 'spatialReference extent')
    sp_info = spatial_info(spatialReference=describe.spatialReference,
                           extent=describe.extent)

    if use_percent:
        percentage = arcpy.GetParameter(4) / 100  # Float
        create_points_from_lines(in_features, out_fc, sp_info.spatialReference,
                                 percent=percentage, add_end_points=end_points,
                                 add_chainage=chainage)
    else:
        distance = arcpy.GetParameterAsText(3)  # String
        distance, param_linear_units = get_distance_and_units(distance)
        distance = convert_units(distance, param_linear_units,
                                 sp_info)

        create_points_from_lines(in_features, out_fc, sp_info.spatialReference,
                                 dist=distance, add_end_points=end_points,
                                 add_chainage=chainage)

    try:
        arcpy.management.AddSpatialIndex(out_fc)
    except arcpy.ExecuteError:
        pass
