"""Replay the seven-dimensional lower connection from the proved mass identity.

This constructs a differential system. It does not certify a nonzero Hodge
extension at a specified fibre, nor does it include the top extension.
"""
from pathlib import Path
import json,hashlib

ROOT=Path(__file__).resolve().parents[1]
path=ROOT/'data/banana_pure_ray/system.json'
data=json.loads(path.read_text(encoding='utf-8'))
weights=[3,2,1,0]
squared_masses=[m*m for m in data['masses']]
coupling=['0']+[f'{r*m}/ss' if r else '0' for r,m in zip(weights,squared_masses)]+['0']
matrix=[row+['0'] for row in data['matrix']]+[coupling+['0']]
out={'status':'derived_from_geometric_mass_lift_identity_and_published_pure_connection',
     'system_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
     'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'basis':'(H_0,H_1,H_2,H_3,H_4,H_5,beta_c)',
     'matrix':matrix,'activity_vector':weights[:3],
     'reference':'DE_RHAM_MINIMALITY.md: nabla_s beta_i = mu_i nabla_mu_i omega',
     'fixed_fibre_nonvanishing_certified':False,'top_extension_included':False,
     'absolute_comparison_certified':False}
(path.parent/'active_lower_system.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print('saved active lower system; no fixed-fibre nonvanishing claim')
