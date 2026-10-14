from pathlib import Path
import urllib.request,json,hashlib,subprocess
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'external/graph_algorithms';out.mkdir(exist_ok=True)
pdf=out/'dual_critical.pdf'
if not pdf.exists():pdf.write_bytes(urllib.request.urlopen('https://arxiv.org/pdf/1410.1653',timeout=35).read())
src=Path('E:/document/feynman_period/energy_character_classification.pdf')
pop='E:/CodexData/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftoppm.exe'
for path,pages,label in [(src,[4,5,8,11],'energy'),(pdf,[1,2,3,4],'dual_critical')]:
 for page in pages:
  dest=out/f'{label}_p{page}'
  if not dest.with_suffix('.png').exists():
   subprocess.run([pop,'-f',str(page),'-l',str(page),'-singlefile','-scale-to','1400','-png',str(path),str(dest)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
manifest={'sources':[{'path':str(src),'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'pages':[4,5,8,11]},
 {'path':str(pdf.relative_to(ROOT)),'url':'https://arxiv.org/pdf/1410.1653','sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'pages':[1,2,3,4]}]}
(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('source pages rendered')

