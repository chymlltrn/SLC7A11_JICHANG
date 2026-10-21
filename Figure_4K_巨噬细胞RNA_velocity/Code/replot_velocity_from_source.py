"""Redraw the historical velocity figure from exported vector plotting source.
Uses Python standard library only. Does not require the original PDF, Seurat,
velocyto.R, gene data, or a server connection.
Usage:
 python replot_velocity_from_source.py SOURCE_DIR/pdf_vector_source.json.gz output.svg
The source stores vector paths, clip paths, point circles, hulls and glyph
outlines. Coordinates are page points, x right and y down. To alter arrows in
UMAP coordinates, use the accompanying calibrated grid-arrow CSV.
"""
import argparse,gzip,json,xml.etree.ElementTree as ET
from pathlib import Path

def obj_to_node(d):
 e=ET.Element(d['tag'],d['attrs']);e.text=d['text'];e.tail=d['tail']
 for c in d['children']:e.append(obj_to_node(c))
 return e

def main():
 p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('output');a=p.parse_args()
 with gzip.open(a.source,'rt') as f:z=json.load(f)
 if z['schema']!='historical-pdf-vector-source-v1':raise ValueError('Unknown source schema')
 ET.register_namespace('','http://www.w3.org/2000/svg')
 ET.register_namespace('xlink','http://www.w3.org/1999/xlink')
 root=obj_to_node(z['vector_tree'])
 Path(a.output).write_bytes(ET.tostring(root,encoding='utf-8',xml_declaration=True))
 print(f"Redrawn: {a.output}; page {z['page_width_pt']} x {z['page_height_pt']} pt")
if __name__=='__main__':main()
