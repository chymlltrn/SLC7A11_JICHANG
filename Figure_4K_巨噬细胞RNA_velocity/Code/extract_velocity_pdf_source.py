"""Extract historical velocity plotting source, with existing-UMAP calibration.
Requires PyMuPDF, NumPy, pandas. No fitted velocity model or gene matrix is used.
"""
import argparse,csv,gzip,json,shutil,xml.etree.ElementTree as ET
from pathlib import Path
import fitz,numpy as np,pandas as pd

def node_to_obj(e):
 return {'tag':e.tag,'attrs':dict(e.attrib),'text':e.text,'tail':e.tail,'children':[node_to_obj(c) for c in e]}
def writecsv(path,rows):
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--inventory',required=True);ap.add_argument('--metadata',required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
 base=Path(a.inventory);meta=Path(a.metadata);out=Path(a.out);summary=[]
 for fig,sub,kind in [('Figure_4K','velocity_macrophages_251104','Macrophages'),('Figure_5I','velocity_Tcells_251104','T_cells')]:
  mm=pd.read_csv(meta/f'{kind}_existing_umap_metadata.csv')
  for group in ['Control','1A4_Antibody']:
   dst=out/fig/f'source_{group}';dst.mkdir(parents=True,exist_ok=True)
   p=base/'1A4_scRNA'/sub/f'velocity_umap_group_{group}_PLOT_only.pdf'
   doc=fitz.open(p);page=doc[-1];draws=page.get_drawings()
   # Initial R pch=21 points are unfilled 4-cubic circle paths, in original cell order.
   circles=[d for d in draws if d['type']=='s' and len(d['items'])==4 and all(q[0]=='c' for q in d['items'])]
   cells=mm[mm.group==group].copy().reset_index(drop=True)
   if len(circles)!=len(cells):raise ValueError(f'{fig} {group}: point/cell count mismatch')
   xy=np.array([[(d['rect'].x0+d['rect'].x1)/2,(d['rect'].y0+d['rect'].y1)/2] for d in circles])
   umap=cells[['UMAP1','UMAP2']].to_numpy();X=np.c_[umap,np.ones(len(umap))]
   affine=np.linalg.lstsq(X,xy,rcond=None)[0];pred=X@affine;res=np.linalg.norm(xy-pred,axis=1)
   if res.max()>.02:raise ValueError(f'{fig} {group}: affine residual too large: {res.max()} pt')
   inverse=np.linalg.inv(affine[:2,:])
   def to_umap(v):return (np.array(v)-affine[2,:])@inverse
   cells['pdf_x_pt']=xy[:,0];cells['pdf_y_pt']=xy[:,1]
   cells.to_csv(dst/'cells_umap_labels.csv',index=False)
   rows=[]
   for d in draws:
    # R lwd=1 -> PDF width .75 pt. Axis arrows use lwd=1.1 -> .825 pt.
    if d['type']=='s' and len(d['items'])==1 and d['items'][0][0]=='l' and abs(d['width']-.75)<1e-5:
     _,p0,p1=d['items'][0];u0=to_umap(p0);u1=to_umap(p1)
     rows.append(dict(arrow_id=len(rows)+1,x0=u0[0],y0=u0[1],x1=u1[0],y1=u1[1],pdf_x0_pt=p0.x,pdf_y0_pt=p0.y,pdf_x1_pt=p1.x,pdf_y1_pt=p1.y,pdf_draw_seqno=d['seqno']))
   writecsv(dst/'grid_arrows_from_original_pdf.csv',rows)
   # Preserve every original vector, clip path and outlined glyph. This is a
   # plotting primitive tree, not a raster, gene matrix, or analysis object.
   svg=page.get_svg_image(text_as_path=True)
   tree=node_to_obj(ET.fromstring(svg))
   payload={'schema':'historical-pdf-vector-source-v1','coordinate_system':'PDF/SVG points; origin top-left; x right; y down',
            'page_width_pt':page.rect.width,'page_height_pt':page.rect.height,'vector_tree':tree}
   with gzip.open(dst/'pdf_vector_source.json.gz','wt') as f:json.dump(payload,f,separators=(',',':'))
   # Human-editable original vector source, useful in Illustrator/Inkscape.
   (dst/'historical_plot_vectors.svg').write_text(svg)
   spans=[]
   for b in page.get_text('dict')['blocks']:
    for line in b.get('lines',[]):
     for s in line['spans']:
      u=to_umap(s['origin'])
      spans.append(dict(text=s['text'],pdf_baseline_x_pt=s['origin'][0],pdf_baseline_y_pt=s['origin'][1],UMAP_baseline_x=u[0],UMAP_baseline_y=u[1],font=s['font'],font_size_pt=s['size'],color=f"#{s['color']:06x}"))
   writecsv(dst/'pdf_text_baselines.csv',spans)
   provenance={'figure':fig,'group':group,'manuscript_group':'IgG2b' if group=='Control' else 'SLC7A11-Ab',
    'server_source_pdf':f'/data/public_database/Home_data/scRNAseq/1A4_anti/1A4_scRNA/{sub}/{p.name}',
    'server_umap_object':f'/data/public_database/Home_data/scRNAseq/1A4_anti/1A4_scRNA/{kind}_final_annotated.rds',
    'source_page':len(doc),'source_page_note':'Original PDFs have a blank page 1; plotted result is page 2.',
    'original_method':'velocyto.R; original scripts Velocity_subgroup_Macro.R / Velocity_subgroup_Tcell.R',
    'server_observed_velocyto_R_version':'0.6',
    'historical_package_version':'not recorded; the current installed package version does not establish the historical version',
    'extraction_method':'Historical PDF vector extraction; least-squares affine calibration against all existing Seurat UMAP cell coordinates in exact original cell order.',
    'cell_count':len(cells),'grid_arrow_count':len(rows),'max_affine_residual_pt':float(res.max()),'rms_affine_residual_pt':float(np.sqrt(np.mean(res**2))),
    'affine_equation':'[UMAP1, UMAP2, 1] @ affine_3x2 = [PDF_x_pt, PDF_y_pt]',
    'affine_3x2':affine.tolist(),
    'precision_note':'Arrow coordinates retain original PDF coordinate precision (roughly 0.01 point); UMAP arrow endpoints are inverse-affine reconstructions, not native saved numerical projection.',
    'scope':'Displayed grid flow, points, hulls and labels recovered. Individual cell velocity vectors and transition probabilities are not present in the historical PDF and are not asserted as recovered.',
    'excluded':'No raw single-cell data, counts/expression matrix, loom, full Seurat, fitted gene-level rvel, new embedding or new velocity fit.'}
   (dst/'provenance_and_calibration.json').write_text(json.dumps(provenance,indent=2))
   summary.append({k:provenance[k] for k in ['figure','group','cell_count','grid_arrow_count','max_affine_residual_pt','rms_affine_residual_pt']})
   # Historical PDF included for comparison, unchanged.
   shutil.copy2(p,dst/'original_PLOT_only.pdf')
 out.mkdir(parents=True,exist_ok=True);(out/'validation_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
