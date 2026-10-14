from pathlib import Path
import urllib.request,hashlib,json
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'external/fixed-pure-group';out.mkdir(exist_ok=True)
sources=[
 ('charles_picard.pdf','https://www.imo.universite-paris-saclay.fr/~francois.charles/Picard.pdf',[2,5,7]),
 ('elsenhans_jahnel_components.pdf','https://wwwuser.gwdguser.de/~jjahnel/Arbeiten/Komponentengruppe_v4f.pdf',[1,20,21])]
records=[]
for name,url,pages in sources:
    p=out/name
    if not p.exists():
        data=urllib.request.urlopen(url,timeout=35).read()
        assert data[:4]==b'%PDF'
        p.write_bytes(data)
    import fitz
    doc=fitz.open(p)
    for page in pages:
        dest=out/(p.stem+'_p'+str(page+1)+'.png')
        if not dest.exists():doc[page].get_pixmap(matrix=fitz.Matrix(1.3,1.3)).save(dest)
    (out/(p.stem+'.txt')).write_text('\n'.join(page.get_text() for page in doc),encoding='utf-8')
    records.append({'name':name,'url':url,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'rendered_pages_one_based':[k+1 for k in pages]})
    print(name,'saved',flush=True)
(out/'manifest.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf-8')

