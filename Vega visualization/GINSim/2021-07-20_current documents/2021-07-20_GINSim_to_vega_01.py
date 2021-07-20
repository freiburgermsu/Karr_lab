""" Utilities for converting metabolic maps from ginsim to Vega format
:Author: Andrew Freiburger <afreiburger@uvic.ca>
:Date: 2021-05-16
"""
import xml.etree.ElementTree as et
import datetime
import json
import math
import numpy
import re
import os

__all__ = ['ginsim_to_vega']

def ginsim_to_vega(ginsim_output_filename, template_name, vega_output_filename):
    """ Convert a metabolic pathway map from GINsim format to Vega format.
    Args:
        ginsim_output_filename (:obj:`str`): path to the ginsim simulation data
        template_name (:obj:`str`): path to the Vega template for the ginsim simulation
        vega_output_filename (:obj:`str`): path to save the Vega map
    """
      
    # parse the GINsim output file into pertinent data
    with open(ginsim_output_filename, 'r') as file:
        tree = et.parse(file)
        root = tree.getroot()

    graph = root[0]
    nodes = []
    regulations = {}
    for child in graph:
        if child.tag == 'node':
            try:
                comment = child.find('annotation/comment').text
            except:
                comment = None
                
            node_coordinates = child.find('nodevisualsetting/ellipse')      
            if node_coordinates == None:
                node_coordinates = child.find('nodevisualsetting/rect')
                
            y_coordinate = node_coordinates.attrib['y']
            x_coordinate = node_coordinates.attrib['x']
            node_width = node_coordinates.attrib['width']
            node_height = node_coordinates.attrib['height']
            node_fill_color = node_coordinates.attrib['backgroundColor']
            text_color = node_coordinates.attrib['foregroundColor']
            nodes.append({'id': child.attrib['id'],
                          'x_value':x_coordinate,
                          'y_value':y_coordinate,
                          'node_width': node_width,
                          'node_height': node_height,
                          'text_color': text_color,
                          'fill_color': node_fill_color,
                          'comment': comment})
                
        elif child.tag == 'edge':
            line_points = child.find('edgevisualsetting/polyline')    
            edge_color = line_points.attrib['line_color']
            edge_width = line_points.attrib['line_width']
            coordinates = line_points.attrib['points']
            regulations[child.attrib['id']] = {'from_node':child.attrib['from'],
                                               'to_node':child.attrib['to'],
                                               'sign': child.attrib['sign'],
                                               'line_points': coordinates,
                                               'line_color': edge_color,
                                               'line_width': edge_width}                
       
        
    # the associated nodes for each node are identified and stored 
    self_loops = []
    for regulation, value in regulations.items():  
        first_node, second_node = regulation.split(':', 1)      
        value['associated_nodes'] = [first_node, second_node]
        if first_node == second_node:
            self_loops.append(regulation)
            
            
    # determine and parameterize the map coordinates
    min_x = math.inf
    min_y = math.inf
    max_x = -math.inf
    max_y = -math.inf
    for node in nodes:
        min_x = min(min_x, float(node['x_value']))
        max_x = max(max_x, float(node['x_value']))
        min_y = min(min_y, float(node['y_value']))
        max_y = max(max_y, float(node['y_value']))     

    ginsim_width = max_x - min_x
    ginsim_height = max_y - min_y
    
    quantity_nodes = len(nodes)
    print('Quantity of nodes', quantity_nodes)
    optimum_nodes_quantity = 18
    nodes_per_map_scale = (quantity_nodes / optimum_nodes_quantity) ** 0.55
    max_width_height = 700 * nodes_per_map_scale

    if ginsim_width > ginsim_height:
        width = max_width_height
        coordinate_scale = float(width / ginsim_width)
        height = ginsim_height * coordinate_scale
    else:
        height = max_width_height
        coordinate_scale = float(height / ginsim_height)
        width = ginsim_width * coordinate_scale
        
    width += 200 
    height += 50
    
    if width > height:
        print('Maximum figure dimension', width)
    else:
        print('Maximum figure dimension', height)
    
    # convert the map nodes to Vega
    edge_y_value_adjustment = 70
    edge_x_value_adjustment = 40
    net_x_adjustment = - min_x + edge_x_value_adjustment    
    net_y_adjustment = - min_y + edge_y_value_adjustment

    '''                            'x': node['x_value'] + net_x_adjustment + label_x_value_adjustment,
                                'y': node['y_value'] + net_y_adjustment + label_y_value_adjustment})
    '''               
    node_width_sum = []
    node_height_sum = []
    node_labels = []    
    for node in nodes:
        adjusted_x_node = (float(node['x_value']) + net_x_adjustment) * coordinate_scale
        adjusted_y_node = (float(node['y_value']) + net_y_adjustment) * coordinate_scale
        node["node_height"] = float(node["node_height"]) *coordinate_scale
        node["node_width"] = float(node["node_width"]) *coordinate_scale
        node["x_value"] = adjusted_x_node
        node["y_value"] = adjusted_y_node        
        
        node_labels.append({'id': node['id'],
                            'text_color': node['text_color'],
                            'x': adjusted_x_node + node["node_width"] / 2,
                            'y': adjusted_y_node + node["node_height"] / 2})
        
        node_width_sum.append(node["node_width"])
        node_height_sum.append(node["node_height"])
        
    empirical_width_adjustment = 6 * nodes_per_map_scale ** 0.8
    empirical_height_adjustment = 5 * nodes_per_map_scale ** 0.8
    average_node_width = sum(node_width_sum) / (len(node_width_sum) * 2) + empirical_width_adjustment
    average_node_height = sum(node_height_sum) / (len(node_height_sum) * 2) + empirical_height_adjustment
        
    # convert regulation paths to Vega    
    self_loop_label_y_offset = -33
    edges_label_y_offset = -10 * nodes_per_map_scale
    arc_y_offset = -20 + nodes_per_map_scale ** 2.7
    arc_x_offset = -10
    edge_coordinates = []
    edge_end_coordinates = []
    regulations_labels = []
    points = []

    for name, value in regulations.items():
        edge_svg_path = ''
        points = value['line_points'].split(' ')
        for point in range(len(points) - 1):
            x1, y1 = points[point].split(',')      
            x1 = (float(x1) + net_x_adjustment) * coordinate_scale
            y1 = (float(y1) + net_y_adjustment) * coordinate_scale
            if name in self_loops:
                node_arc_width = 20 
                node_arc_radius = 20
                edge_end_x = x1 + arc_x_offset
                edge_end_y = y1 + arc_y_offset
                edge_svg_path += 'M {} {} a {} {} 0 1 1 {} 0'.format(edge_end_x, edge_end_y, node_arc_radius, node_arc_radius, node_arc_width)
                break
            
            else:               
                if point == 0:
                    edge_svg_path += 'M {} {}'.format(x1, y1)
                    x_value = x1
                    y_value = y1
                        
                if point + 1 <= len(points) - 1:
                    x2, y2 = points[point + 1].split(',')
                    x2 = (float(x2) + net_x_adjustment) * coordinate_scale
                    y2 = (float(y2) + net_y_adjustment) * coordinate_scale
                    
                    if len(points) == 3:
                        delta_x = x2 - x1
                        delta_y = y2 - y1
                    else:
                        delta_x = x2 - x_value
                        delta_y = y2 - y_value
                    
                    if abs(delta_x) < abs(delta_y):
                        if point + 1 == len(points) - 1:
                            if delta_y < 0:
                                end_offset = (nodes_per_map_scale - 1) * (average_node_height * 0.7)
                                edge_svg_path += ' v {}'.format(delta_y + end_offset)
                            elif delta_y > 0:
                                end_offset = (0.6 - nodes_per_map_scale) * average_node_height * 0.7
                                edge_svg_path += ' v {}'.format(delta_y + end_offset)
                        else:
                            edge_svg_path += ' v {}'.format(delta_y)
                            
                        delta_x = 0
                        edge_end_x = x_value
                        edge_end_y = y_value = y2

                    elif abs(delta_x) > abs(delta_y):
                        if point + 1 == len(points) - 1:
                            if delta_x < 0:
                                end_offset = (nodes_per_map_scale - 1) * (average_node_width / 4)
                                edge_svg_path += ' h {}'.format(delta_x + end_offset)
                            elif delta_x > 0:
                                end_offset = (1 - nodes_per_map_scale) * (average_node_width / 4)
                                edge_svg_path += ' h {}'.format(delta_x + end_offset)
                        else:
                            edge_svg_path += ' h {}'.format(delta_x)
                            
                        delta_y = 0
                        edge_end_y = y_value
                        edge_end_x = x_value = x2
                        
                    else:
                        print('ERROR: The {} edge is inconsistent with a cardinal direction.'.format(name))
                    
                else:
                    print('ERROR: The {} contains an undefined point'.format(name))
                        
                '''if point + 1 == len(points) - 1:
                    if delta_x == 0:
                        edge_end_x = x1
                    elif delta_y == 0:
                        edge_end_y = y1'''
                    
        # establishing the edge ends on the figure
        if value['sign'] == 'positive':
            edge_end_shape = 'triangle'
            edge_end_angle = math.degrees(math.atan2(delta_y, delta_x))
            stroke_width = 1
        elif value['sign'] == 'negative':
            edge_end_shape = 'stroke'
            edge_end_angle = math.degrees(math.atan2(delta_y, delta_x)) + 90
            stroke_width = 5
        elif value['sign'] == 'unknown':
            edge_end_shape = 'circle'
            edge_end_angle = math.degrees(math.atan2(delta_y, delta_x))
            stroke_width = 2
                   
        if delta_x > 0:
            edge_end_x_buffer = -average_node_width
            edge_end_y_buffer = 0
            if edge_end_shape == 'triangle':
                edge_end_angle -= 150
        elif delta_x < 0:
            edge_end_x_buffer = average_node_width
            edge_end_y_buffer = 0
            if edge_end_shape == 'triangle':
                edge_end_angle += 210
        if delta_y > 0:
            edge_end_x_buffer = 0
            edge_end_y_buffer = -average_node_height
            if edge_end_shape == 'triangle':
                edge_end_angle += 90
        elif delta_y < 0:
            edge_end_x_buffer = 0
            edge_end_y_buffer = average_node_height
            if edge_end_shape == 'triangle':
                edge_end_angle += 90
            
                    
        if name in self_loops:
            labels_y_offset = self_loop_label_y_offset
            edge_ends_y_offset = -5
            edge_ends_x_offset = node_arc_width + 5
            edge_end_angle = -135
        elif name not in self_loops:
            labels_y_offset = edges_label_y_offset        
            edge_ends_y_offset = edge_end_y_buffer
            edge_ends_x_offset = edge_end_x_buffer
                        
        regulations_labels.append({'id': name,
                                   'line_color': value['line_color'],
                                   'x': x1,
                                   'y': y1 + labels_y_offset})
            
            
        edge_end_coordinates.append({
            'name': name,
            'x': edge_end_x + edge_ends_x_offset,
            'y': edge_end_y + edge_ends_y_offset,
            'edge_end_angle': edge_end_angle,
            'edge_end_shape': edge_end_shape,
            'stroke_width': stroke_width,
            'line_color': value['line_color']
        })
        
        edge_coordinates.append({
            "name": name,
            'from_node': value['from_node'],
            'to_node': value['to_node'],
            'edge_svg_path': edge_svg_path,
            'line_color': value['line_color'],
            'line_width': value['line_width'],
            'related_nodes': value['associated_nodes']
        })
            
    # incorporate the calculated details into the Vega template
    'https://run.api.biosimulations.org/results/60a1a1754c88d864b0aa4dd2/simulation.sedml%2Freport_wt?sparse=false'
    
    vega = json.load(open(template_name))
        
    signal_height = 100
    vega['width'] = width
    vega['height'] = height + signal_height
    
    url = 'https://run.api.biosimulations.org/results/60ca299e1bdfd54d622a6297/simulation.sedml%2Freport_wt?sparse=false'

    for entry in vega['data']:
        if entry['name'] == 'nodesData':
            entry['values'] = nodes
    
        elif entry['name'] == 'edgesData':
            entry['values'] = edge_coordinates
        
        elif entry['name'] == 'edgeEndCoordinatesData':
            entry['values'] = edge_end_coordinates
            
        elif entry['name'] == 'edgesLabelsData':
            entry['values'] = regulations_labels
            
        elif entry['name'] == 'nodesLabelsData':
            entry['values'] = node_labels
            
        elif entry['name'] == 'nodesValues':
            entry['url'] = url
    
    for signal in vega['signals']:
        if signal['name'] == 'edgeEndStrokeWidth':
            signal['value'] = 1.5 * coordinate_scale
    
        elif signal['name'] == 'signalHeight':
            signal['value'] = signal_height
    
        elif signal['name'] == 'signalPadding':
            signal['value'] = 0

        elif signal['name'] == 'nodeStrokeWidthData':
            signal['value'] = 1 * coordinate_scale

        elif signal['name'] == 'edgeStrokeWidthData':
            signal['value'] = 18 * coordinate_scale
    
        elif signal['name'] == 'mapMaxX':
            signal['value'] = width

        elif signal['name'] == 'mapMaxY':
            signal['value'] = height
            

    # save Vega-formatted map
    with open(vega_output_filename, 'w') as file:
        json.dump(vega, file, indent = 4)

# visualize GINSim Model 35
'''ginsim_to_vega('regulatoryGraph.ginml',
               'ginsim_to_vega.template_static_5.json',
               '{}_APF_Irons ginsim to vega, static_10.json'.format(datetime.date.today()))
'''
        
# visualize GINSim Model 79
ginsim_to_vega('regulatoryGraph_79_expanded comments.ginml',
               'ginsim_to_vega.template_static_5.json',
               '{}_APF_Ginsim to vega, static,79_1.json'.format(datetime.date.today()))