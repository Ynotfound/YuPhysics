from pathlib import Path
import urllib.request,hashlib,json,subprocess
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'external/phi4_benchmark';out.mkdir(exist_ok=True)
url='https://pdg.lbl.gov/2025/reviews/rpp2025-rev-kinematics.pdf';p=out/'pdg_kinematics.pdf'
if not p.exists():p.write_bytes(urllib.request.urlopen(url,timeout=30).read())
original=Path('E:/document/feynman_period/feynman_periods_expanded.pdf')
pop='E:/CodexData/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftoppm.exe'
for path,pages,prefix in [(p,[2,3,6,7],'pdg'),(original,[39],'original')]:
    for page in pages:
        dest=out/f'{prefix}_p{page}'
        if not dest.with_suffix('.png').exists():
            subprocess.run([pop,'-f',str(page),'-l',str(page),'-singlefile','-scale-to','1400','-png',str(path),str(dest)],check=True,capture_output=True)
manifest={'sources':[{'path':str(p.relative_to(ROOT)),'url':url,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'pages':[2,3,6,7]},
 {'path':str(original),'sha256':hashlib.sha256(original.read_bytes()).hexdigest(),'pages':[39]}]}
(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Original benchmark and PDG convention pages saved')

